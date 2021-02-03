#!/usr/bin/env bash

PROJECT_DIR=$(pwd)

# CGNS
CGNS_BUILD_DIR=${PROJECT_DIR}/dependencies/build/CGNS
CGNS_DIR=${PROJECT_DIR}/dependencies/CGNS

mkdir -p $CGNS_BUILD_DIR
cd $CGNS_BUILD_DIR
cmake $CGNS_DIR
make -j
cd $PROJECT_DIR
