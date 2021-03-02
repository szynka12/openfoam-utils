add_executable(deltaSolver   deltaSolver.C)

target_include_directories(deltaSolver PUBLIC
  $ENV{FOAM_SRC}/OpenFOAM/lnInclude
  $ENV{FOAM_SRC}/OSspecific/POSIX/lnInclude
  $ENV{FOAM_SRC}/finiteVolume/lnInclude
  $ENV{FOAM_SRC}/meshTools/lnInclude
  $ENV{FOAM_SRC}/sampling/lnInclude
  $ENV{FOAM_SRC}/TurbulenceModels/turbulenceModels/lnInclude
  $ENV{FOAM_SRC}/TurbulenceModels/incompressible/lnInclude
  $ENV{FOAM_SRC}/transportModels
  $ENV{FOAM_SRC}/transportModels/incompressible/singlePhaseTransportModel
  $ENV{FOAM_SRC}/dynamicMesh/lnInclude
  $ENV{FOAM_SRC}/dynamicFvMesh/lnInclude
  ${CMAKE_CURRENT_LIST_DIR}
  )

target_link_directories(deltaSolver  PUBLIC
  $ENV{FOAM_LIBBIN}
  $ENV{FOAM_LIBBIN}/$ENV{FOAM_MPI}
  )

target_link_libraries(deltaSolver PUBLIC
  OpenFOAM
  Pstream
  finiteVolume
  fvOptions
  meshTools
  sampling
  turbulenceModels
  incompressibleTurbulenceModels
  incompressibleTransportModels
  dynamicMesh
  dynamicFvMesh
  topoChangerFvMesh
  atmosphericModels
  )
