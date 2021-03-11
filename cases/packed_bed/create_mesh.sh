#!/usr/bin/env zsh

nx=1
ny=2
nz=1

set_block () {
  x=$1
  y=$2
  z=$3
  
  # normals
  i="(1 0 0)"
  j="(0 1 0)"
  k="(0 0 1)"

  case_name="${x}_${y}_${z}"

  if [ -d "$case_name" ]; then rm -Rf $case_name; fi
  cp -r "quarter_template" $case_name
  
  foamGetDict -case "$case_name" -f mirrorMeshDict 

  blm="$case_name"/system/blockMeshDict
  mm=${case_name}/system/mirrorMeshDict

  # add i, j, k to the blockMeshDict, cant use set because it makes macro
  # expansion of our variables fucked
  sed -i "s/index1/${x}/g" $blm
  sed -i "s/index2/${y}/g" $blm
  sed -i "s/index3/${z}/g" $blm
  
  or=$(foamDictionary $blm -entry geometry.sphere.origin -value)

  foamDictionary ${mm} -entry pointAndNormalDict.point -set "$or"

  blockMesh -case "$case_name"

  foamDictionary ${mm} -entry pointAndNormalDict.normal -set "${i}"
  mirrorMesh -overwrite -case "$case_name"

  foamDictionary ${mm} -entry pointAndNormalDict.normal -set "${j}"
  mirrorMesh -overwrite -case "$case_name"

  foamDictionary ${mm} -entry pointAndNormalDict.normal -set "${k}"
  mirrorMesh -overwrite -case "$case_name"

  splitByRegion -case "$case_name" -noDict -patch x
  sed -i "s/x1/xp_${case_name}/g" "$case_name"/constant/polyMesh/boundary
  sed -i "s/x3/xm_${case_name}/g" "$case_name"/constant/polyMesh/boundary
  
  splitByRegion -case "$case_name" -noDict -patch y
  sed -i "s/y1/yp_${case_name}/g" "$case_name"/constant/polyMesh/boundary
  sed -i "s/y3/ym_${case_name}/g" "$case_name"/constant/polyMesh/boundary
  
  splitByRegion -case "$case_name" -noDict -patch z
  sed -i "s/z1/zp_${case_name}/g" "$case_name"/constant/polyMesh/boundary
  sed -i "s/z3/zm_${case_name}/g" "$case_name"/constant/polyMesh/boundary
}

for a in {1..$nx}
do
  for b in {1..$ny}
  do
    for c in {1..$nz}
    do
      set_block $a $b $c
      
      if [[ "$a$b$c" != "111" ]]; then
        case_name="${a}_${b}_${c}"
        mergeMeshes -overwrite "1_1_1" $case_name
        stitchMesh -case "1_1_1" "yp_1_1_1" "ym_${case_name}"
      fi
    
    done
  done
done
