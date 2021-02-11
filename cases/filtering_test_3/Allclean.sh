#!/usr/bin/env bash

rm -r processor*

foamListTimes -rm
foamListTimes -withZero -rm -case constant/targetMesh

foamCleanPolyMesh 
foamCleanPolyMesh -case constant/targetMesh
