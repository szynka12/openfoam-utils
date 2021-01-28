
#include <map>

#include "fvCFD.H"
#include "meshTools.h"
#include "polyMesh.H"

label new_cell_number(const labelList&              decrements,
                      const std::map<label, label>& deleted_cells,
                      const label                   cell_i)
{
  if (deleted_cells.find(cell_i) != deleted_cells.end())
  {
    return deleted_cells.at(cell_i) - decrements[deleted_cells.at(cell_i)];
  }
  else
  {
    return cell_i + decrements[cell_i];
  }
};

int main(int argc, char* argv[])
{
  Foam::argList::noParallel();

  // clang-format off
  #include "addOverwriteOption.H"
  #include "setRootCase.H"
  #include "createTime.H"
  #include "createMesh.H"
  // clang-format on

  const bool overwrite = args.found("overwrite");

  IOdictionary settings(IOobject("mirrorMeshDict", runTime.system(), mesh,
                                 IOobject::MUST_READ, IOobject::NO_WRITE));

  vector point;
  vector normal;

  if (settings.get<word>("planeType") == "pointAndNormal")
  {
    point  = settings.subDict("pointAndNormalDict").get<vector>("point");
    normal = settings.subDict("pointAndNormalDict").get<vector>("normal");
  }
  else
  {
    FatalError << "Plane type must be 'pointAndNormal' for this tool." << endl;
  }

  // Find all internal faces that are:
  //  a) coplanar to the mirror plane
  //  b) lie on the plane

  auto eps_compare = [](double lhs, double rhs) {
    return (lhs > rhs - SMALL) && (lhs < rhs + SMALL);
  };

  faceList faces = mesh.faces();

  auto face_centers = mesh.faceCentres();
  auto face_normals = mesh.faceAreas() / mag(mesh.faceAreas());

  auto marked_faces = boolList(faces.size(), false);

  for (label i = 0; i < mesh.nInternalFaces(); i++)
  {
    marked_faces[i]
        = eps_compare(mag(dot(face_normals()[i], normal)), 1.0)
          && eps_compare(mag(dot(face_centers[i] - point, normal)), 0.0);
  }

  labelList              new_cell_labels(mesh.nCells(), 0);
  std::map<label, label> deleted_cells;

  for (label i = 0; i < mesh.nInternalFaces(); i++)
  {
    if (marked_faces[i])
    {
      const auto owner = mesh.faceOwner()[i];
      const auto nei   = mesh.faceNeighbour()[i];

      // mark the the neighbour for deletion
      deleted_cells[nei] = owner;

      // decrease the number in each cell above the deleted one by 1
      for (label i = nei + 1; i < mesh.nCells(); i++) { new_cell_labels[i]--; }
    }
  }

  auto owner     = mesh.faceOwner();
  auto neighbour = mesh.faceNeighbour();

  for (label i = 0; i < mesh.nFaces(); i++)
  {
    owner[i] = new_cell_number(new_cell_labels, deleted_cells, owner[i]);

    if (i < mesh.nInternalFaces())
    {
      neighbour[i]
          = new_cell_number(new_cell_labels, deleted_cells, neighbour[i]);
      // if (owner[i] > neighbour[i]) {faces[i].flip();}
    }
  }

  inplaceSubset(marked_faces, faces, true);
  inplaceSubset(marked_faces, owner, true);
  inplaceSubset(marked_faces, neighbour, true);

  auto cell_info = mt::cell_neighbours(owner, neighbour);
  auto bad_cells = mt::find_multiply_connected_cells(cell_info.second);
  auto bad_faces = mt::repair_multiply_connected_cells(cell_info, bad_cells,
                                                       owner, neighbour, faces);

  label n_bad_faces = 0;
  forAll(bad_faces, i)
  {
    if (!bad_faces[i]) n_bad_faces++;
  }

  // delete all the bad faces
  Foam::Info << "\tRemoving " << n_bad_faces << " faces..." << Foam::endl;
  Foam::inplaceSubset(bad_faces, faces);
  Foam::inplaceSubset(bad_faces, owner);
  Foam::inplaceSubset(bad_faces, neighbour);
  Foam::Info << "\tDone!" << Foam::endl;

  auto points = mesh.points();

  Foam::Info << "Creating polyMesh ..." << Foam::endl;
  // Create mesh form components, patches will be added later
  Foam::polyMesh new_mesh(
      Foam::IOobject(Foam::polyMesh::defaultRegion, runTime.constant(), runTime,
                     Foam::IOobject::NO_READ,
                     Foam::IOobject::AUTO_WRITE  // this must be here! (NO_WRITE
                                                 // is a default)
                     ),
      std::move(points), std::move(faces), std::move(owner),
      std::move(neighbour));

  Foam::Info << "Done!" << Foam::endl;

  // count the marked faces (we will need to decrease bc definitions by this
  // amount)
  label n_marked_faces = 0;
  forAll(marked_faces, i)
  {
    if (marked_faces[i]) { n_marked_faces++; }
  }

  Foam::PtrList<Foam::polyPatch> patch_list;

  forAll(mesh.boundaryMesh(), i)
  {
    const auto& bc_patch = mesh.boundaryMesh()[i];
    auto        patch
        = bc_patch.clone(mesh.boundaryMesh(), i, bc_patch.size(),
                         bc_patch.start() - n_marked_faces - n_bad_faces);
    patch_list.append(patch);
  }

  new_mesh.addPatches(patch_list);

  if (!overwrite)
  {
    ++runTime;
    new_mesh.setInstance(runTime.timeName());
  }

  new_mesh.write();
}
