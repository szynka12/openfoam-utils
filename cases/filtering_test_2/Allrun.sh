#!/usr/bin/env bash

blockMesh -case constant/targetMesh
gmshToFoam cylinders.msh
bash setBoundaryTypes.sh
pisoFoam

