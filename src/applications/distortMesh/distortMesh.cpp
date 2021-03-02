
#include "fvCFD.H"

int main(int argc, char* argv[])
{
  Foam::argList::noParallel();
// Foam::argList::addArgument("distort magnitude [0-1]");

// clang-format off
  #include "setRootCase.H"
  #include "createTime.H"
  #include "createMesh.H"
  #include "addOverwriteOption.H"
  // clang-format on
  
  const bool overwrite = args.found("overwrite");

  IOdictionary settings(IOobject("distortMeshDict", runTime.system(),
                                 mesh, IOobject::MUST_READ,
                                 IOobject::NO_WRITE));

  const auto internal_dist_f = settings.get<double>("internal");

  const auto boundary_dist_f = settings.get<double>("boundary");

  const auto ignore_features = 
    settings.getOrDefault<bool>("ignoreFeatures", false);



  const auto n_points = mesh.nPoints();
  labelList boundary_indices(n_points, -1);
  boolList is_feature_point(n_points, false);


  auto& b_mesh     = mesh.boundaryMesh();

  forAll(b_mesh, bi)
  {
    const auto point_labels = b_mesh[bi].meshPoints();
    forAll(point_labels, pi)
    {
      const auto p = point_labels[pi];
      if (boundary_indices[p] == -1) // was internal
      {
        boundary_indices[p] = bi;
      }
      else // was not internal, meaning it must be a feature
      {
        is_feature_point[p] = true;
      }
    }
  }

  const auto point_cells = mesh.pointCells();
  const scalarField h = Foam::pow( mesh.cellVolumes(), 1. / 3.)();
  scalarField dual_cell_dimension(n_points, 0);
  
  forAll(dual_cell_dimension, i)
  {
    const auto& cells = point_cells[i];
    forAll(cells, ci)
    {
      dual_cell_dimension[i] += h[cells[ci]];
    }
    dual_cell_dimension[i] /= cells.size();
  }

  auto r_gen = Foam::Random();

  pointField points = mesh.points();

  forAll( points , pi)
  {
    if(!is_feature_point[pi] || ignore_features)
    {
      const auto r_v = r_gen.sample01<Foam::vector>();
      
      if (boundary_indices[pi] == -1) //internal
      {
        points[pi] += internal_dist_f * dual_cell_dimension[pi] * r_v;
      }
    }
  }

  mesh.movePoints(points);

  //if (!overwrite)
  //{
    //++runTime;
    //mesh.setInstance(runTime.timeName());
  //}

  runTime.setTime(instant(runTime.constant()), 0);
  
  mesh.write();

  return 0;

  // compute the mean cell size and use it to set the distortion magnitude
  //auto distort_mag = 0.5 * settings.get<double>("distortmagnitude")
                     //* Foam::pow(Foam::average(mesh.cellVolumes()), 1. / 3.);

  //auto& b_mesh     = mesh.boundaryMesh();
  //auto  boundaries = settings.getOrDefault("boundaries", b_mesh.names());

  //// We only deal with quads here so, if we count a point 4 times it means it is
  //// in the "interior" of the boundary
  //Foam::List<Foam::label> points_count(mesh.nPoints(), 0);

  //auto             r_gen = Foam::Random();
  //Foam::pointField points(mesh.points());

  //Foam::Info << "Modifing points..." << Foam::endl;

  //forAll(b_mesh, bi)
  //{
    //if (boundaries.found(b_mesh[bi].name()))
    //{
      //const auto boundary_normal = b_mesh[bi][0].unitNormal(mesh.points());

      //Foam::List<Foam::label> points_count(mesh.nPoints(), 0);

      //forAll(b_mesh[bi], fi)
      //{
        //const auto& face = b_mesh[bi][fi];
        //forAll(face, pi) { points_count[face[pi]]++; }
      //}

      //forAll(points, pi)
      //{
        //if (points_count[pi] == 4)
        //{
          //const auto r_v = r_gen.sample01<Foam::vector>();
          //const auto mv  = r_v - boundary_normal * (boundary_normal & r_v);
          //points[pi]     = mesh.points()[pi] + distort_mag * mv;
        //}
        //points_count[pi] = 0;
      //}
    //}
  //}

  //mesh.movePoints(points);

  //runTime.setTime(instant(runTime.constant()), 0);
  //mesh.write();

  //Foam::Info << "End!" << Foam::endl;
  //return 0;
}
