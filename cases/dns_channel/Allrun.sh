mv 0.org 0
blockMesh
renumberMesh -overwrite
perturbUChannel
decomposePar
mpirun -n 24 pimpleFoam -parallel > log.pimpleFoam &


