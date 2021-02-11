#!/usr/bin/env bash

foamListTimes -rm
foamListTimes -withZero -rm -case constant/targetMesh

foamCleanPolyMesh 
foamCleanPolyMesh -case constant/targetMesh
