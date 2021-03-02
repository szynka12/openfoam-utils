add_executable(renumberCoupledPatches renumberCoupledPatches.cpp)

target_include_directories(renumberCoupledPatches PUBLIC
  $ENV{FOAM_SRC}/finiteVolume/lnInclude
  $ENV{FOAM_SRC}/meshTools/lnInclude
  $ENV{FOAM_SRC}/OpenFOAM/lnInclude
  $ENV{FOAM_SRC}/OSspecific/POSIX/lnInclude
  )

target_link_directories(renumberCoupledPatches PUBLIC
  $ENV{FOAM_LIBBIN}
  $ENV{FOAM_LIBBIN}/$ENV{FOAM_MPI}
  )


target_link_libraries(renumberCoupledPatches PUBLIC 
  finiteVolume
  meshTools
  OpenFOAM
  fileFormats
  surfMesh
  Pstream
  )
