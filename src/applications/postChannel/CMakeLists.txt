add_executable(userPostChannel userPostChannel.C channelIndex.C)

target_include_directories(userPostChannel PUBLIC
  $ENV{FOAM_SRC}/finiteVolume/lnInclude
  $ENV{FOAM_SRC}/meshTools/lnInclude
  $ENV{FOAM_SRC}/OpenFOAM/lnInclude
  $ENV{FOAM_SRC}/OSspecific/POSIX/lnInclude
  $ENV{FOAM_SRC}/sampling/lnInclude
  )

target_link_directories(userPostChannel PUBLIC
  $ENV{FOAM_LIBBIN}
  $ENV{FOAM_LIBBIN}/openmpi-system
  )


target_link_libraries(userPostChannel PUBLIC 
  finiteVolume
  meshTools
  OpenFOAM
  fileFormats
  genericPatchFields
  sampling
  Pstream
  )

