add_executable(distortBoundaryMesh distortBoundaryMesh.cpp)

target_include_directories(distortBoundaryMesh PUBLIC
  $ENV{FOAM_SRC}/finiteVolume/lnInclude
  $ENV{FOAM_SRC}/meshTools/lnInclude
  $ENV{FOAM_SRC}/OpenFOAM/lnInclude
  $ENV{FOAM_SRC}/OSspecific/POSIX/lnInclude
  )

target_link_directories(distortBoundaryMesh PUBLIC
  $ENV{FOAM_LIBBIN}
  $ENV{FOAM_LIBBIN}/$ENV{FOAM_MPI}
  )

target_link_libraries(distortBoundaryMesh PUBLIC 
  finiteVolume
  meshTools
  OpenFOAM
  fileFormats
  surfMesh
  Pstream
  )
