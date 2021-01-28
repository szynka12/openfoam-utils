

#include "fvCFD.H"
#include "scalarList.H"

int main(int argc, char* argv[])
{
  Foam::argList::noParallel();
// Foam::argList::addArgument("distort magnitude [0-1]");

// clang-format off
  #include "setRootCase.H"
  #include "createTime.H"
  #include "createMesh.H"
  // clang-format on

  const cellList cells     = mesh.cells();
  const auto     c_centres = mesh.cellCentres();
  const auto     f_centres = mesh.faceCentres();
  const auto     faces     = mesh.faces();
  const auto     own       = mesh.faceOwner();
  const auto     nei       = mesh.faceNeighbour();

  List<scalar> oq(cells.size(), 0.0);

  forAll(oq, i)
  {
    const auto cur_faces = cells[i];
    const auto c_c       = c_centres[i];

    List<scalar> oq_faces(cur_faces.size(), 0.0);

    // normal is in the directionm of a neighbour
    forAll(cur_faces, fi)
    {
      face face = faces[cur_faces[fi]];
      if (own[cur_faces[fi]] != i) { face.flip(); }

      const auto c_f = f_centres[cur_faces[fi]];

      oq_faces[fi] = (c_f - c_c)
                     & face.areaNormal(mesh.points()) / mag(c_f - c_c)
                           / mag(face.areaNormal(mesh.points()));

      if (cur_faces[fi] < mesh.nInternalFaces())
      {
        label cell_neighbour;
        if (own[cur_faces[fi]] == i) { cell_neighbour = nei[cur_faces[fi]]; }
        else if (nei[cur_faces[fi]] == i)
        {
          cell_neighbour = own[cur_faces[fi]];
        }
        else
        {
          Info << "Dupa" << endl;
        }
        auto c_nei = c_centres[cell_neighbour];

        scalar oq_cell = (c_nei - c_c)
                         & face.areaNormal(mesh.points()) / mag(c_nei - c_c)
                               / mag(face.areaNormal(mesh.points()));

        oq_faces[fi] = min(oq_faces[fi], oq_cell);
      }
    }
    oq[i] = min(oq_faces);
  }

  scalar min_oq = min(oq);

  Info << "Min. orthogonal quality = " << min_oq << endl;
  Info << "Corresponding angle = "
       << Foam::acos(min_oq) * 180.0 / Foam::constant::mathematical::pi << endl;

  return 0;
}
