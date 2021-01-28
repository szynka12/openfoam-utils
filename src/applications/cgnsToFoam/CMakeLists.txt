add_executable(cgnsToFoam cgnsToFoam.cpp)



target_include_directories(cgnsToFoam PUBLIC
  ${CGNS_INC_DIR}
  ${CGNS_LIB_DIR}
  lnInclude
  $ENV{FOAM_SRC}/finiteVolume/lnInclude
  $ENV{FOAM_SRC}/meshTools/lnInclude
  $ENV{FOAM_SRC}/OpenFOAM/lnInclude
  $ENV{FOAM_SRC}/OSspecific/POSIX/lnInclude
  )

target_link_directories(cgnsToFoam PUBLIC
  ${CGNS_LIB_DIR}
  $ENV{FOAM_LIBBIN}
  )

target_link_libraries(cgnsToFoam PUBLIC 
  cgns
  finiteVolume
  meshTools
  OpenFOAM
  )
