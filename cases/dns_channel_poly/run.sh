rm -r 0
cp -r 0.org 0
perturbUChannel
#decomposePar
#mpirun -n 24 pimpleFoam -parallel > log.pimpleFoam 2>&1 &


