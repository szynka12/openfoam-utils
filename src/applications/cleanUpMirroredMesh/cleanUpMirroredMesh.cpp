
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
  Foam::argList::addBoolOption("fast");

  // clang-format off
  #include "addOverwriteOption.H"
  #include "setRootCase.H"
  #include "createTime.H"
  #include "createMesh.H"
  // clang-format on

  const bool overwrite = args.found("overwrite");
  const bool cell_subset = args.found("fast");

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
  
  Foam::label n_del_faces = 0; 
  Foam::Info << "Marking faces for deletion..." << Foam::endl;
  
  for (label i = 0; i < mesh.nInternalFaces(); i++)
  {
    marked_faces[i]
        = eps_compare(mag(dot(face_normals()[i], normal)), 1.0)
          && eps_compare(mag(dot(face_centers[i] - point, normal)), 0.0);
    if (marked_faces[i]) n_del_faces++;
  }

  Foam::Info << "Done!" << Foam::endl;
  Foam::Info << "Marking cells for deletion..." << Foam::endl;

  labelList              new_cell_labels(mesh.nCells(), 0);
  labelList              cells_near_plane;
  std::map<label, label> deleted_cells;

  for (label i = 0; i < mesh.nInternalFaces(); i++)
  {
    if (marked_faces[i])
    {
      const auto owner = mesh.faceOwner()[i];
      const auto nei   = mesh.faceNeighbour()[i];

      // mark the the neighbour for deletion
      deleted_cells[nei] = owner;

      // Decrease the number in each cell above the deleted one by 1
      //
      //for (label i = nei + 1; i < mesh.nCells(); i++) { new_cell_labels[i]--; }
      //
      // This loop costs a lot on big meshes: for every marked face in (hundred 
      // thousand) traverse a cell list (in milion). The same can be achieved by
      // decrementing only the found cell labels and than updating the
      // information in the cell list, all at once!

      if (nei + 1 < mesh.nCells()) { new_cell_labels[nei + 1]--; }
      if (cell_subset) 
      {
        cells_near_plane.append(nei);
        cells_near_plane.append(owner);
      }
    }
  }
  Foam::Info << "Done!" << Foam::endl;
  
  Foam::Info << "Updating the cell decrements list..." << Foam::endl;
  
  Foam::label cumulated_decrement = 0;
  for (label i = 0; i < mesh.nCells(); i++)
  {
    if (new_cell_labels[i] == -1) 
    {
      new_cell_labels[i] += cumulated_decrement;
      cumulated_decrement--;
    }
    else
    {
      new_cell_labels[i] += cumulated_decrement;
    }
  }
  Foam::Info << "Done!" << Foam::endl;
  
  Foam::Info << "Updating cell numbers ..." << Foam::endl;
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
  if (cell_subset)
  {
    forAll(cells_near_plane, celli)
    {
      cells_near_plane[celli] 
          = new_cell_number(
              new_cell_labels, deleted_cells, cells_near_plane[celli]);
    }
  }

  Foam::Info << "Done!" << Foam::endl;
  


  Foam::Info << "Removing " << n_del_faces << " faces..." << Foam::endl;
  inplaceSubset(marked_faces, faces, true);
  inplaceSubset(marked_faces, owner, true);
  inplaceSubset(marked_faces, neighbour, true);
  Foam::Info << "Done!" << Foam::endl;


  auto cell_info = mt::cell_neighbours(owner, neighbour);
  Foam::Info << cell_info.first.size() << Foam::endl;
  if (cell_subset)
  {
    Foam::boolList subset_list(cell_info.first.size(), false);
    for (const auto& c : cells_near_plane)
    {
      subset_list[c] = true;
    }
    //inplaceSubset(subset_list, cell_info.first);
    inplaceSubset(subset_list, cell_info.second);
  }
  Foam::Info << cell_info.first.size() << Foam::endl;

  auto bad_cells = mt::find_multiply_connected_cells(cell_info.second);
  auto bad_faces = mt::repair_multiply_connected_cells(cell_info, bad_cells,
                                                       owner, neighbour, faces);

  label n_bad_faces = 0;
  forAll(bad_faces, i)
  {
    if (!bad_faces[i]) n_bad_faces++;
  }

  // delete all the bad faces
  Foam::Info << "Removing " << n_bad_faces 
             << " faces (belonging to multiply connected cells)..." 
             << Foam::endl;
  Foam::inplaceSubset(bad_faces, faces);
  Foam::inplaceSubset(bad_faces, owner);
  Foam::inplaceSubset(bad_faces, neighbour);
  Foam::Info << "Done!" << Foam::endl;

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
