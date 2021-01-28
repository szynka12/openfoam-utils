
#include "fvCFD.H"

int main(int argc, char* argv[])
{
  Foam::argList::noParallel();
// Foam::argList::addArgument("distort magnitude [0-1]");

// clang-format off
  #include "setRootCase.H"
  #include "createTime.H"
  #include "createMesh.H"
  // clang-format on

  IOdictionary settings(IOobject("distortBoundaryMeshDict", runTime.system(),
                                 mesh, IOobject::MUST_READ,
                                 IOobject::NO_WRITE));

  // compute the mean cell size and use it to set the distortion magnitude
  auto distort_mag = 0.5 * settings.get<double>("distortMagnitude")
                     * Foam::pow(Foam::average(mesh.cellVolumes()), 1. / 3.);

  auto& b_mesh     = mesh.boundaryMesh();
  auto  boundaries = settings.getOrDefault("boundaries", b_mesh.names());

  // We only deal with quads here so, if we count a point 4 times it means it is
  // in the "interior" of the boundary
  Foam::List<Foam::label> points_count(mesh.nPoints(), 0);

  auto             r_gen = Foam::Random();
  Foam::pointField points(mesh.points());

  Foam::Info << "Modifing points..." << Foam::endl;

  forAll(b_mesh, bi)
  {
    if (boundaries.found(b_mesh[bi].name()))
    {
      const auto boundary_normal = b_mesh[bi][0].unitNormal(mesh.points());

      Foam::List<Foam::label> points_count(mesh.nPoints(), 0);

      forAll(b_mesh[bi], fi)
      {
        const auto& face = b_mesh[bi][fi];
        forAll(face, pi) { points_count[face[pi]]++; }
      }

      forAll(points, pi)
      {
        if (points_count[pi] == 4)
        {
          const auto r_v = r_gen.sample01<Foam::vector>();
          const auto mv  = r_v - boundary_normal * (boundary_normal & r_v);
          points[pi]     = mesh.points()[pi] + distort_mag * mv;
        }
        points_count[pi] = 0;
      }
    }
  }

  mesh.movePoints(points);

  runTime.setTime(instant(runTime.constant()), 0);
  mesh.write();

  Foam::Info << "End!" << Foam::endl;
  return 0;
}
