// TODO splitting into patches
// TODO refactor of the header (might put all the functions here)
// TODO check other element types (how fluent exports them)
// TODO add ability to disable checks
// TODO mesh scaling
// TODO dictionary

// Cgns headers
#include "cgnslib.h"

// Foam headers
#include "IOobject.H"
#include "ListOps.H"
#include "Time.H"
#include "argList.H"
#include "autoPtr.H"
#include "cellList.H"
#include "faceList.H"
#include "pointField.H"
#include "polyMesh.H"
#include "polyPatch.H"
#include "wallPolyPatch.H"

// helper functions
#include "cgnsToFoam.h"

int main(int argc, char *argv[])
{
  Foam::argList::noParallel();
  Foam::argList::addArgument(".cgns file");
  Foam::argList::addOption("scale", "scalar", "scale the mesh");

// clang-format off
  #include "setRootCase.H"
  #include "createTime.H"
  // clang-format on

  int file = cgh::open_cgns(args[1].c_str());

  bool split = args.found("regSplit");
  Foam::scalar scale = args.getOrDefault<Foam::scalar>("scale", 1.0);

  cgh::Base b(file, 1);
  b.info();

  cgh::Zone z(file, b.n, 1);
  z.info();

  Foam::Info << "Reading coordinates ..." << Foam::endl;

  std::vector<double> cx = cgh::read_coord<double>(file, b, z, 1);
  std::vector<double> cy = cgh::read_coord<double>(file, b, z, 2);
  std::vector<double> cz = cgh::read_coord<double>(file, b, z, 3);

  Foam::pointField points(z.n_nodes);
  forAll(points, i) { points[i] = Foam::vector(cx[i], cy[i], cz[i]) * scale; }

  Foam::Info << "Done!" << Foam::endl;

  // First read information about boundary conditions

  Foam::Info << "Reading boundary conditions ..." << Foam::endl;
  int bc_number;
  cgh::cgns_check_error(cg_nbocos(file, b.n, z.n, &bc_number));

  std::vector<cgh::Boundary> bcs;
  for (int bc_id = 1; bc_id <= bc_number; bc_id++)
  {
    bcs.push_back(cgh::Boundary(file, b.n, z.n, bc_id));
    bcs.back().info();
  }
  Foam::Info << "Done!" << Foam::endl;

  int n_sections;
  cgh::cgns_check_error(cg_nsections(file, b.n, z.n, &n_sections));

  std::pair<int, int> N       = cgh::count_faces_and_cells(file, b.n, z.n);
  int                 n_faces = N.first;
  int                 n_cells = N.second;

  Foam::Info << "In " << n_sections << " sections there is:" << Foam::endl;
  Foam::Info << "\tFaces: " << n_faces << Foam::endl;
  Foam::Info << "\tCells: " << n_cells << Foam::endl;

  Foam::Info << "Reading sections ..." << Foam::endl;

  // Preallocate the containers
  Foam::faceList           face_list(n_faces);
  Foam::faceList::iterator f_it = face_list.begin();

  Foam::cellList           cell_list(n_cells);
  Foam::cellList::iterator c_it = cell_list.begin();

  for (int section_i = 1; section_i <= n_sections; section_i++)
  {
    char          sectionname[33];
    ElementType_t type;
    cgsize_t      start, end;
    int           nbndry, parent_flag;
    cgh::cgns_check_error(cg_section_read(file, b.n, z.n, section_i,
                                          sectionname, &type, &start, &end,
                                          &nbndry, &parent_flag));

    Foam::Info << "\tSection '" << sectionname << "'..." << Foam::endl;

    int n_elements = end - start + 1;

    if (type == CGNS_ENUMV(NFACE_n))
    {
      cgh::read_ngon<Foam::cell, 0>(file, b.n, z.n, n_elements, section_i,
                                    c_it);
    }
    else if (type == CGNS_ENUMV(NGON_n))
    {
      cgh::read_ngon<Foam::face>(file, b.n, z.n, n_elements, section_i, f_it);
    }
    else
    {
      Foam::FatalError << "Unsupported element type: " << type << Foam::endl;
    }
  }
  Foam::Info << "Done!" << Foam::endl;

  Foam::Info << "Computing the neighbour/owner ..." << Foam::endl;

  // Create owner neighbour tables
  Foam::labelList owner(face_list.size(), -1);
  Foam::labelList neighbour(face_list.size(), -1);
  Foam::boolList  marked_faces(face_list.size(), false);

  Foam::label n_internal_faces = 0;
  bool        should_flip      = false;
  Foam::label face_            = -1;

  // Check internal faces for such situation that two cells are sharing two
  // faces For example:
  //
  // 0 -- 1 --- 2
  // |    |     |
  // |    6     |
  // |    |     |
  // 3 -- 4 --- 5
  //
  // This can probably happen using the poly-hexacore approach from fluent. So
  // let's declare the container for the cell neighbours.
  Foam::List<Foam::List<Foam::label>> cell_nei(n_cells);

  // Iterate over the container in reverse so to generate correct addressing.
  // The mesh form fluent is written 'backwards' in terms of FOAM mesh. So we
  // will have to reverse every single list. Here the change of index is also
  // necessary to make the reversal happen (cell_list.size() - 1 - celli).

  forAllReverse(cell_list, celli)
  {
    // we will reverse the ordering of the cells so this is important
    const Foam::label rev_celli = cell_list.size() - 1 - celli;

    // get reference to face labels for current cell
    Foam::labelList &cellfaces = cell_list[celli];

    forAll(cellfaces, facei)
    {
      // If cellfaces[facei] < 0 is true then the normal vector points inward
      // the cell. We will use that information to check if faces needs flipping
      if (cellfaces[facei] < 0)
      {
        // cells still have the CGNS numbering in them
        cellfaces[facei] = -cellfaces[facei] - 1;
        should_flip      = true;
      }
      else
      {
        cellfaces[facei] = cellfaces[facei] - 1;
        should_flip      = false;
      }
      face_ = cellfaces[facei];

      if (!marked_faces[face_])
      {
        // First visit: owner
        owner[face_] = rev_celli;

        marked_faces[face_] = true;

        // Normal points to the cell with bigger number (no neighbour). If the
        // condition is true we need to flip
        if (should_flip) { face_list[face_].flip(); }
      }
      else
      {
        // Second visit: neighbour
        neighbour[face_] = rev_celli;
        n_internal_faces++;

        // We also know that the current cell is a neighbour to the owner of
        // face 'face_'
        cell_nei[owner[face_]].append(rev_celli);
        cell_nei[rev_celli].append(owner[face_]);
      }
    }
  }

  Foam::Info << "\t# Internal faces: " << n_internal_faces << Foam::endl;

  // Sort every list, so that finding duplicates would be easy
  std::for_each(cell_nei.begin(), cell_nei.end(),
                [](Foam::labelList &l) { std::sort(l.begin(), l.end()); });

  // Find multiply connected cells (as in the drawing). Lists are sorted so now
  // we simple compare the consecutive elements. We store here ordered pairs
  // (lower cell is first).
  Foam::List<std::pair<Foam::label, Foam::label>> bad_cells;
  forAll(cell_nei, celli)
  {
    const auto &curr_nei = cell_nei[celli];
    for (Foam::label neii = 1; neii < curr_nei.size(); neii++)
    {
      if (curr_nei[neii - 1] == curr_nei[neii])
      {
        auto bcell = celli < curr_nei[neii]
                         ? std::make_pair(celli, curr_nei[neii])
                         : std::make_pair(curr_nei[neii], celli);

        bad_cells.append(bcell);
      }
    }
  }

  Foam::inplaceUniqueSort(bad_cells);

  Foam::label n_bad_faces = 0;
  if (bad_cells.size())
  {
    Foam::Info << "\t\t"
               << "Found " << bad_cells.size()
               << " pair(s) of multiply "
                  "connected cells:"
               << Foam::endl;
    Foam::Info << "\t\t" << bad_cells << Foam::endl;
    Foam::Info << "\t\t"
               << "Attempting correction..." << Foam::endl;

    // We will use this lambda to reverse the reordering of cell indices that
    // was done previously to ensure the upper triangular ordering of faces. We
    // need it because we want to get the faces from the (already) renumbered
    // cells.
    auto reverse_mapping = [n_cells](Foam::label rev_celli) {
      return -(rev_celli - n_cells + 1);
    };

    Foam::boolList bad_faces(face_list.size(), true);
    forAll(bad_cells, bi)
    {
      // grab the cells in question
      const auto &c1 = cell_list[reverse_mapping(bad_cells[bi].first)];
      const auto &c2 = cell_list[reverse_mapping(bad_cells[bi].second)];

      Foam::List<Foam::label> multiple_faces;

      forAll(c1, facei)
      {
        if (owner[c1[facei]] == bad_cells[bi].first
            && neighbour[c1[facei]] == bad_cells[bi].second)
        {
          multiple_faces.append(c1[facei]);
        }
      }
      forAll(c2, facei)
      {
        if (owner[c2[facei]] == bad_cells[bi].second
            && neighbour[c2[facei]] == bad_cells[bi].first)
        {
          multiple_faces.append(c2[facei]);
        }
      }
      Foam::Info << "\t\t"
                 << "Correcting faces: " << multiple_faces << Foam::endl;

      // First cell will remain the rest has to be merged. They have the same
      // (and correct) orientation so we look for repeated vertices, and stitch
      // a face using that information:
      // 1 -- 2 -- 3    1 -- 2 -- 3
      // |    |    | -> |         |
      // 4 -- 5 -- 6    4 -- 5 -- 6
      // so we must identify the indices and now when to switch
      Foam::face &recipient_face = face_list[multiple_faces[0]];
      for (Foam::label facei = 1; facei < multiple_faces.size(); facei++)
      {
        const Foam::face &donor_face     = face_list[multiple_faces[facei]];
        bad_faces[multiple_faces[facei]] = false;
        n_bad_faces++;

        // This is dirty hack used to exit two nested loops
        bool loop_switch = false;

        forAll(recipient_face, rpi)
        {
          if (loop_switch) { break; }
          forAll(donor_face, dpi)
          {
            if (recipient_face[rpi] == donor_face[dpi])
            {
              // we have found a same point, that means that the other one for
              // recipient it is the next index and for a donor it is previous.
              // The new cell will be created from indices
              // [rpi + 1, rpi] + [dpi + 1 , dpi - 2 ]
              Foam::face merged_face(recipient_face.size() + donor_face.size()
                                     - 2);

              for (Foam::label mi = 0; mi < recipient_face.size(); mi++)
              {
                rpi++;
                if (rpi == recipient_face.size()) { rpi = 0; }
                merged_face[mi] = recipient_face[rpi];
              }

              for (Foam::label mi = recipient_face.size();
                   mi < merged_face.size(); mi++)
              {
                dpi++;
                if (dpi == donor_face.size()) { dpi = 0; }
                merged_face[mi] = donor_face[dpi];
              }

              recipient_face = merged_face;

              loop_switch = true;
              break;
            }
          }
        }
      }
    }
    // delete all the bad faces
    Foam::Info << "\t\tRemoving " << n_bad_faces << " faces..." << Foam::endl;
    Foam::inplaceSubset(bad_faces, face_list);
    Foam::inplaceSubset(bad_faces, owner);
    Foam::inplaceSubset(bad_faces, neighbour);
    Foam::Info << "\t\tDone!" << Foam::endl;
  }

  // All bad faces must have been internal
  n_internal_faces -= n_bad_faces;
  n_faces -= n_bad_faces;

  // Reverse and shift numbering in boundary lists
  for (auto &b : bcs)
  {
    std::for_each(std::begin(b.faces), std::end(b.faces),
                  [n_faces, n_bad_faces](cgsize_t &x) { x = n_faces - 1 - x; });
    std::reverse(std::begin(b.faces), std::end(b.faces));
  }

  // Sort the boundary patches in terms of the faces indices, we need to do that
  // since OpenFOAM expects that the patches are in order (I think)
  std::sort(bcs.begin(), bcs.end());

  // Check if all the boundary arrays are contiguous together
  int first_b_face = bcs.front().faces.front();
  int last_b_face  = bcs.back().faces.back();
  int n_b_faces    = n_faces - n_internal_faces;

  bool contaigous_boundary
      = (first_b_face == n_internal_faces) && (last_b_face == n_faces - 1);

  if (!contaigous_boundary)
  {
    Foam::Warning << "Something is wrong with cell numbering in the original "
                     "mesh.\n Run mesh checks in *both* Fluent Mesher and "
                     "Solver to make sure that the mesh is correct."
                  << Foam::endl;
  }

  // Reverse all the stuff. Now neighbours are in a correct order for a resize
  std::reverse(face_list.begin(), face_list.end());
  std::reverse(owner.begin(), owner.end());
  std::reverse(neighbour.begin(), neighbour.end());
  neighbour.resize(n_internal_faces);

  Foam::Info << "Done!" << Foam::endl;

  Foam::Info << "Creating polyMesh ..." << Foam::endl;
  // Create mesh form components, patches will be added later
  Foam::polyMesh mesh(
      Foam::IOobject(Foam::polyMesh::defaultRegion, runTime.constant(), runTime,
                     Foam::IOobject::NO_READ,
                     Foam::IOobject::AUTO_WRITE  // this must be here! (NO_WRITE
                                                 // is a default)
                     ),
      std::move(points), std::move(face_list), std::move(owner),
      std::move(neighbour));

  Foam::Info << "Done!" << Foam::endl;

  Foam::Info << "Attaching patches ..." << Foam::endl;

  Foam::PtrList<Foam::polyPatch> patch_list(bcs.size());

  for (int i = 0; i < bcs.size(); i++)
  {
    Foam::autoPtr<Foam::polyPatch> patch_ptr;
    patch_ptr = Foam::autoPtr<Foam::polyPatch>(
        new Foam::wallPolyPatch(bcs[i].name, bcs[i].size(), bcs[i].faces[0], i,
                                mesh.boundaryMesh(), "wall"));

    if (patch_ptr) { patch_list.set(i, patch_ptr); }
  }
  mesh.addPatches(patch_list);
  Foam::Info << "Done!" << Foam::endl;

  Foam::Info << "Writing mesh!" << Foam::endl;
  mesh.removeFiles();
  mesh.write();

  Foam::Info << "End!" << Foam::endl;

  return 0;
}
