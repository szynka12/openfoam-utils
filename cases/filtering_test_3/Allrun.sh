#!/usr/bin/env bash

blockMesh -case constant/targetMesh
gmshToFoam cylinders.msh
bash setBoundaryTypes.sh

decomposePar
broadcastTargetMesh
mpirun -n 2 pisoFoam -parallel
gatherTargetMesh


