#include "OFstream.H"
#include "argList.H"
#include "fvCFD.H"
#include "polyMesh.H"

// TODO refactor that later
bool check_faces(const label& f1, const label& f2,
                 const vector& first_to_second, const vectorField& face_centers,
                 const scalar& tol)
{
  const auto diff
      = mag(face_centers[f2] - (face_centers[f1] + first_to_second));

  return diff < tol;
}

bool comp_vec(const vector& v1, const vector& v2, const vector& first_to_second,
              const scalar& tol)
{
  const auto diff = mag(v2 - (v1 + first_to_second));

  return diff < tol;
}

void renumber_neighbours(labelList&           renumeration_list,
                         const vector&        slave_to_master,
                         const labelListList& master_connectivity,
                         const labelListList& slave_connectivity,
                         const label& master_root, const label& slave_root,
                         const labelRange&  master_range,
                         const labelRange&  slave_range,
                         const vectorField& face_centers,
                         const scalar&      tolerance)
{
  const auto& master_neighbours
      = master_connectivity[master_root - master_range.first()];
  const auto& slave_neighbours
      = slave_connectivity[slave_root - master_range.first()];

  forAll(master_neighbours, mi)
  {
    const auto target_face = master_neighbours[mi] + master_range.first();

    const label target_index = slave_range.first() + master_neighbours[mi];

    forAll(slave_neighbours, si)
    {
      const auto current_slave_face
          = slave_neighbours[si] + slave_range.first();
      if (renumeration_list[current_slave_face] == -1)
      {
        if (check_faces(current_slave_face, target_face, slave_to_master,
                        face_centers, tolerance))
        {
          renumeration_list[current_slave_face] = target_index;

          renumber_neighbours(renumeration_list, slave_to_master,
                              master_connectivity, slave_connectivity,
                              target_face, current_slave_face, master_range,
                              slave_range, face_centers, tolerance);

          break;
        }
      }
    }
  }
}

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

  const auto& b_mesh = mesh.boundaryMesh();

  // those are copies
  pointField points            = mesh.points();
  faceList   faces             = mesh.faces();
  auto       renumeration_list = labelList(faces.size(), -1);

  const vectorField& face_centers = mesh.faceCentres();

  IOdictionary settings(IOobject("renumberCoupledPatchesDict", runTime.system(),
                                 mesh, IOobject::MUST_READ,
                                 IOobject::NO_WRITE));

  // Read dictionary
  PtrList<dictionary> patches(settings.lookup("patches"));
  const scalar        tolerance = settings.lookupOrDefault("tolerance", SMALL);

  for (const dictionary& dict : patches)
  {
    word  master_name = dict.get<word>("master");
    label master_id   = mesh.boundaryMesh().findPatchID(master_name);

    word  slave_name = dict.get<word>("slave");
    label slave_id   = mesh.boundaryMesh().findPatchID(slave_name);

    vector slave_to_master = dict.get<vector>("slaveToMaster");

    auto master_connectivity = b_mesh[master_id].faceFaces();
    auto slave_connectivity  = b_mesh[slave_id].faceFaces();

    auto master_range = b_mesh[master_id].range();
    auto slave_range  = b_mesh[slave_id].range();

    // we do not use this for now
    /*
    // find the frist matching cell, than traverse connectivity while checking
    // we only renumber the slave

    label first_master = master_range.first();
    label first_slave  = -1;

    for (label i = slave_range.first(); i <= slave_range.last(); i++)
    {
      if (check_faces(i, first_master, slave_to_master, face_centers,
                      tolerance))
      {
        first_slave = i;
      }
    }

    renumeration_list[first_slave] = slave_range.first();

    // recursive find matching faces
    renumber_neighbours(renumeration_list, slave_to_master, master_connectivity,
                        slave_connectivity, first_master, first_slave,
                        master_range, slave_range, face_centers, tolerance);
    */

    for (label mi = master_range.first(); mi <= master_range.last(); mi++)
    {
      auto target_slave_label
          = slave_range.first() + (mi - master_range.first());
      for (label si = slave_range.first(); si <= slave_range.last(); si++)
      {
        if (renumeration_list[si] == -1)
        {
          if (check_faces(si, mi, slave_to_master, face_centers, tolerance))
          {
            renumeration_list[si] = target_slave_label;

            const auto& p1_m = points[faces[mi][0]];
            // cycle the slave face until they have same points ordering
            while (!comp_vec(points[faces[si][0]], p1_m, slave_to_master,
                             tolerance))
            {
              inplaceRotateList<List, label>(faces[si], -1);
            }
          }
        }
      }
    }

    // check if all the faces in the range recieved the number
    for (label i = slave_range.first(); i <= slave_range.last(); i++)
    {
      if (renumeration_list[i] == -1)
      {
        FatalError
            << "One (or more) of the slave faces is not renumbered on patch '"
            << slave_name << "'!" << endl;
      }
    }
  }

  // fil the rest of the list
  forAll(renumeration_list, i)
  {
    if (renumeration_list[i] == -1)
    {
      renumeration_list[i] = i;
      i++;
    }
  }

  // renumber faces and owner
  inplaceReorder(renumeration_list, faces);
  labelList owner = reorder(renumeration_list, mesh.faceOwner());

  Foam::Info << "Creating polyMesh ..." << Foam::endl;
  // Create mesh form components, patches will be added later
  Foam::polyMesh new_mesh(
      Foam::IOobject(Foam::polyMesh::defaultRegion, runTime.constant(), runTime,
                     Foam::IOobject::NO_READ,
                     Foam::IOobject::AUTO_WRITE  // this must be here! (NO_WRITE
                                                 // is a default)
                     ),
      std::move(points), std::move(faces), std::move(owner),
      std::move(mesh.faceNeighbour().clone().ref()));

  Foam::Info << "Done!" << Foam::endl;

  if (!overwrite)
  {
    ++runTime;
    new_mesh.setInstance(runTime.timeName());
  }

  // Insert boundary patches
  Foam::PtrList<Foam::polyPatch> patch_list;

  forAll(mesh.boundaryMesh(), i)
  {
    const auto& bc_patch = mesh.boundaryMesh()[i];
    auto        patch    = bc_patch.clone(mesh.boundaryMesh());
    patch_list.append(patch);
  }

  new_mesh.addPatches(patch_list);

  Info << "Writing mesh..." << endl;
  new_mesh.write();

  Info << "End!" << endl;
  return 0;
}
