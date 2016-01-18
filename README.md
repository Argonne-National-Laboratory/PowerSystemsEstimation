# PowerSystemsEstimation
Source files for the 9-bus power system estimation.

To run these codes one needs to:

1. If there is no local petsc git checkout, go to:

http://www.mcs.anl.gov/petsc/developers/index.html

and follow the instructions there to obtain the development version of
PETSc.

2. Configure: (adapt the following to include your path)

./configure --PETSC_ARCH=arch-mac --PETSC_DIR=/Users/noemi/research/workgit/ps-test/petsc --download-exodusii --download-hdf5 --download-mpich=1 --download-netcdf --download-saws --download-sundials=yes --with-cc=gcc --with-fc=gfortran --with-mpi4py

3. Follow instructions on the PETSc tutorial website to do make and
test your local petsc distribution

Note: if you want to make only parts of the code (if running on mac) do:
./arch-mac/lib/petsc/conf/reconfigure-arch-mac.py

4. In your .bashrc or .profile add the following:

export PETSC_DIR=path-to-petsc
export PETSC_ARCH=arch-mac

5. To run the adjoint-based gradient code go to your local dir and do:

make ex9businertiaest_adj
./ex9businertiaest_adj -tao_monitor

or if would like to load a file containing the observations, do:
./ex9businertiaest_adj -tao_monitor -loadObservations obs-perturbed.bin

or if you'd like to see the covariance parameters do:
./ex9businertiaest_adj -outputCov

6. For postprocessing you can use matlab. To start matlab with PETSc
do:

MATLABPATH=$PETSC_DIR/share/petsc/matlab matlab -nodesktop

Then in your matlab code use 'PetscBinaryRead' function to read bin
files, for example A = PetscBinaryRead('out_pert.bin');

7. To run the estimation code with parameters set at the command line,
type for example:

./ex9businertiaest_adj -tao_monitor -data_noise 0.01 -load_disturb 7.0 -tfinal 2.0
