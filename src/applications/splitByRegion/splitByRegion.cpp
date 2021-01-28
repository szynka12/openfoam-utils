#include "argList.H"
#include "fvCFD.H"
#include "polyMesh.H"
#include "wallPolyPatch.H"

// TODO clean the name value pairs and make the idices increment correctly

bool any_valid_neighbour(const labelListList& face_nei,
                         const boolList& marked_faces, const label face)
{
  bool is_valid = false;
  forAll(face_nei[face], f)
  {
    is_valid = is_valid | marked_faces[face_nei[face][f]];
  }
  return is_valid;
}

void traverse_neighbours(const labelListList& face_nei, boolList& marked_faces,
                         const label& face_c, const label& face_p)
{
  forAll(face_nei[face_c], f)
  {
    const label neighbour = face_nei[face_c][f];
    if (!marked_faces[neighbour])
    {
      marked_faces[neighbour] = true;
      traverse_neighbours(face_nei, marked_faces, neighbour, face_c);
    }
  }
}

void fill_renumeration_list(List<label>& list, label& index,
                            const label last_index)
{
  while (index <= last_index)
  {
    list[index] = index;
    index++;
  }
}

int main(int argc, char* argv[])
{
  Foam::argList::noParallel();

  Foam::argList::addBoolOption(
      "noDict",
      "Ignore the dictionary file, provide a patch in '-patch' option.");
  Foam::argList::addOption("patch", "name",
                           "Use provided patch name when -noDict is used.");

  // clang-format off
  #include "setRootCase.H"
  #include "createTime.H"
  #include "createMesh.H"
  // clang-format on

  auto&      b_mesh = mesh.boundaryMesh();
  List<word> boundaries;

  if (!args.found("noDict"))
  {
    IOdictionary settings(IOobject("splitByRegionDict", runTime.system(), mesh,
                                   IOobject::MUST_READ, IOobject::NO_WRITE));
    boundaries = settings.getOrDefault("boundaries", b_mesh.names());
  }
  else
  {
    boundaries.append(args.get("patch"));
  }

  auto n_split_boundaries = boundaries.size();

  Foam::PtrList<Foam::polyPatch>
      patch_list;  //(b_mesh.size() + n_split_boundaries);

  pointField points = mesh.points();
  faceList   faces  = mesh.faces();
  labelList  own    = mesh.faceOwner();
  labelList  nei    = mesh.faceNeighbour();

  List<label> renumeration_list(faces.size(), -1);

  label renumeration_label = 0;
  fill_renumeration_list(renumeration_list, renumeration_label,
                         mesh.nInternalFaces() - 1);

  List<std::pair<string, labelRange>> name_start_list;

  forAll(b_mesh, bi)
  {
    const auto& boundary_mesh = b_mesh[bi];
    auto        range         = boundary_mesh.range();
    if (boundaries.found(boundary_mesh.name()))
    {
      // get the face neighbours
      const auto& face_conn         = boundary_mesh.faceFaces();
      label       new_patches_label = 0;
      label       first_face        = range.first();

      // Start iterating over all the faces
      label    fi = 0;
      boolList marked_faces(boundary_mesh.size(), false);

      if (!marked_faces[fi])
      {
        marked_faces[fi]  = true;
        label current_fi  = fi;
        label previous_fi = -1;

        // go into the tree as deep as you can.
        traverse_neighbours(face_conn, marked_faces, current_fi, previous_fi);
        // now all faces that are connected are marked, that means that we
        // should sort them (let's not modify the order only cherry pick the
        // faces that we need)

        // count the marked faces
        label n_marked_faces = 0;
        forAll(marked_faces, i)
        {
          if (marked_faces[i]) { n_marked_faces++; }
        }

        // create a new range
        name_start_list.append(std::make_pair(
            boundary_mesh.name() + std::to_string(new_patches_label++),
            labelRange(first_face, n_marked_faces)));

        label lower_range_index = first_face;

        Foam::autoPtr<Foam::polyPatch> patch_ptr;
        patch_ptr = Foam::autoPtr<Foam::polyPatch>(new Foam::wallPolyPatch(
            boundary_mesh.name() + std::to_string(new_patches_label++),
            n_marked_faces, first_face, 0, mesh.boundaryMesh(), "wall"));

        if (patch_ptr) { patch_list.append(patch_ptr); }

        first_face += n_marked_faces;

        // new range for the second half of the patch
        name_start_list.append(std::make_pair(
            boundary_mesh.name() + std::to_string(new_patches_label++),
            labelRange(first_face, range.size() - n_marked_faces)));

        patch_ptr = Foam::autoPtr<Foam::polyPatch>(new Foam::wallPolyPatch(
            boundary_mesh.name() + std::to_string(new_patches_label++),
            range.size() - n_marked_faces, first_face, 0, mesh.boundaryMesh(),
            "wall"));

        if (patch_ptr) { patch_list.append(patch_ptr); }
        // fill renumeration list:
        label upper_range_index = first_face;

        forAll(marked_faces, i)
        {
          marked_faces[i]
              ? renumeration_list[renumeration_label] = lower_range_index++
              : renumeration_list[renumeration_label] = upper_range_index++;
          renumeration_label++;
        }
      }
    }
    else
    {
      fill_renumeration_list(renumeration_list, renumeration_label,
                             range.last());
      name_start_list.append(std::make_pair(boundary_mesh.name(), range));

      auto patch = boundary_mesh.clone(b_mesh);
      patch_list.append(patch);
    }
  }

  // Info << name_start_list << endl;
  // Info << patch_list << endl;

  // renumeration information is now present, now we must create a new mesh and
  // overwrite the old once
  inplaceReorder(renumeration_list, faces);
  inplaceReorder(renumeration_list, own);

  Foam::Info << "Creating polyMesh ..." << Foam::endl;
  // Create mesh form components, patches will be added later
  Foam::polyMesh new_mesh(
      Foam::IOobject(Foam::polyMesh::defaultRegion, runTime.constant(), runTime,
                     Foam::IOobject::NO_READ,
                     Foam::IOobject::AUTO_WRITE  // this must be here! (NO_WRITE
                                                 // is a default)
                     ),
      std::move(points), std::move(faces), std::move(own), std::move(nei));

  Foam::Info << "Done!" << Foam::endl;

  Foam::Info << "Attaching patches ..." << Foam::endl;

  new_mesh.addPatches(patch_list);
  Foam::Info << "Done!" << Foam::endl;

  Foam::Info << "Writing mesh!" << Foam::endl;
  new_mesh.removeFiles();
  new_mesh.write();

  Foam::Info << "End!" << Foam::endl;
}
