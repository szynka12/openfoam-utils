cmake_minimum_required (VERSION 3.5)
project (my_project)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_definitions( 
  -DWM_LABEL_SIZE=$ENV{WM_LABEL_SIZE}
  -DNoRepository
  -DOPENFOAM=2006
  -DWM_DP
  -DWM_LABEL_SIZE=32
  )

if(NOT DEFINED ENV{WM_PROJECT_DIR})
  message( FATAL_ERROR "OpenFOAM is not sourced (no $WM_PROJECT_DIR env. \
    variable has been found)!" )
endif()

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ENV{FOAM_USER_APPBIN})
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ENV{FOAM_USER_LIBBIN})

add_subdirectory(functionObjects)


add_custom_target(
    copy-compile-commands ALL
    ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_BINARY_DIR}/compile_commands.json
        ${CMAKE_CURRENT_SOURCE_DIR}
    )
