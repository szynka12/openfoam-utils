#!/usr/bin/env bash

mesh_name="half_cube8k5.cgns"



entry="foamDictionary constant/polyMesh/boundary -entry entry0"

rename_patch () {
  tmp=$(foamDictionary constant/polyMesh/boundary -entry entry0.${1} -value)
  foamDictionary constant/polyMesh/boundary -entry entry0.${1} -remove
  foamDictionary constant/polyMesh/boundary -entry entry0.${2} -set "${tmp}"
}


foamCleanPolyMesh
cgnsToFoam $mesh_name
splitByRegion -noDict -patch "xm-periodic"
mirrorMesh -overwrite
cleanUpMirroredMesh -overwrite
combinePatchFaces 90 -overwrite
foamDictionary constant/polyMesh/boundary -entry entry0.ym-wall -remove
splitByRegion -noDict -patch "yp-wall"


# rename all boundaries
rename_patch "zp-wall" "zp"
rename_patch "zm-wall" "zm"
rename_patch "xm-periodic1" "xm"
rename_patch "xm-periodic3" "xp"
rename_patch "yp-wall1" "yp"
rename_patch "yp-wall3" "ym"

renumberCoupledPatches -overwrite

# set up cyclic boundaries
${entry}.xm.type -set "cyclic"
${entry}.xm.inGroups -set "1(cycylic)"
${entry}.xm.matchTolerance -set "0.0001"
${entry}.xm.transform -set "translational"
${entry}.xm.separationVector -set "(1 0 0)"
${entry}.xm.neighbourPatch -set "xp"

${entry}.xp.type -set "cyclic"
${entry}.xp.inGroups -set "1(cycylic)"
${entry}.xp.matchTolerance -set "0.0001"
${entry}.xp.transform -set "translational"
${entry}.xp.separationVector -set "(-1 0 0)"
${entry}.xp.neighbourPatch -set "xm"

${entry}.ym.type -set "cyclic"
${entry}.ym.inGroups -set "1(cycylic)"
${entry}.ym.matchTolerance -set "0.0001"
${entry}.ym.transform -set "translational"
${entry}.ym.separationVector -set "(0 1 0)"
${entry}.ym.neighbourPatch -set "yp"

${entry}.yp.type -set "cyclic"
${entry}.yp.inGroups -set "1(cycylic)"
${entry}.yp.matchTolerance -set "0.0001"
${entry}.yp.transform -set "translational"
${entry}.yp.separationVector -set "(0 -1 0)"
${entry}.yp.neighbourPatch -set "ym"

renumberMesh -overwrite

checkMesh >> "checkMesh_${mesh_name}.log"
