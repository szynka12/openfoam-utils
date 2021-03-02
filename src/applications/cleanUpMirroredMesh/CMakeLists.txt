add_executable(cleanUpMirroredMesh cleanUpMirroredMesh.cpp)

target_include_directories(cleanUpMirroredMesh PUBLIC
  $ENV{FOAM_SRC}/finiteVolume/lnInclude
  $ENV{FOAM_SRC}/meshTools/lnInclude
  $ENV{FOAM_SRC}/OpenFOAM/lnInclude
  $ENV{FOAM_SRC}/OSspecific/POSIX/lnInclude
  )

target_link_directories(cleanUpMirroredMesh PUBLIC
  $ENV{FOAM_LIBBIN}
  $ENV{FOAM_LIBBIN}/$ENV{FOAM_MPI}
  )


target_link_libraries(cleanUpMirroredMesh PUBLIC 
  finiteVolume
  meshTools
  OpenFOAM
  fileFormats
  surfMesh
  Pstream
  )
