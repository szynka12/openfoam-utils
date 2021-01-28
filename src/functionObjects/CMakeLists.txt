set(SOURCES 
  turbulentKineticEnergy/turbulentKineticEnergy.C
  )

add_library(userFunctionObjects SHARED ${SOURCES})

target_include_directories(userFunctionObjects PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  $ENV{FOAM_SRC}/OpenFOAM/lnInclude
  $ENV{FOAM_SRC}/OSspecific/POSIX/lnInclude
  $ENV{FOAM_SRC}/finiteVolume/lnInclude
  $ENV{FOAM_SRC}/functionObjects/field/lnInclude

  )

target_link_directories(userFunctionObjects PUBLIC
  # OpenFOAM libraries
  $ENV{FOAM_LIBBIN}
  $ENV{FOAM_LIBBIN}/openmpi-system
  )

target_link_libraries(userFunctionObjects PUBLIC 
  OpenFOAM
  )
  
