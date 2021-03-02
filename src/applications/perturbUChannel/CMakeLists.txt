add_executable(perturbUChannel perturbUChannel.C)

target_include_directories(perturbUChannel PUBLIC
  $ENV{FOAM_SRC}/meshTools/lnInclude
  $ENV{FOAM_SRC}/finiteVolume/lnInclude
  $ENV{FOAM_SRC}/OpenFOAM/lnInclude
  $ENV{FOAM_SRC}/OSspecific/POSIX/lnInclude
  )

target_link_directories(perturbUChannel PUBLIC
  $ENV{FOAM_LIBBIN}
  $ENV{FOAM_LIBBIN}/$ENV{FOAM_MPI}
  )


target_link_libraries(perturbUChannel PUBLIC 
  finiteVolume
  meshTools
  OpenFOAM
  Pstream
  )

