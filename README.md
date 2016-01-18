# PowerSystemsEstimation
Source files for the 9-bus power system estimation.

To run these codes one needs to:

1. If there is no local petsc git checkout, go to:

http://www.mcs.anl.gov/petsc/developers/index.html

and follow the instructions there to obtain the development version
of PETSc.

****IMPORTANT****

The time integrator in the lates dev files has a bug. The gradient and
adjoints work with Hong's branch:

hongzh/ts-adjoint-over-split-intervals

****IMPORTANT****

2. Configure: (adapt the following to include your path)

./configure --PETSC_ARCH=arch-mac --PETSC_DIR=/Users/noemi/research/workgit/ps-test/petsc --download-exodusii --download-hdf5 --download-mpich=1 --download-netcdf --download-saws --download-sundials=yes --with-cc=gcc --with-fc=gfortran --with-mpi4py

3. Follow instructions to do make and test your local petsc distribution

Note: if you want to make only parts of the code (if running on mac) do:
./arch-mac/lib/petsc/conf/reconfigure-arch-mac.py

4. In your .bashrc or .profile add the following:

export PETSC_DIR= path-to-petsc
export PETSC_ARCH=arch-mac

5. To run the adjoint-based gradient code go to your local dir (to
state-estimation/codes/paramest_9bus) and do:

make ex9busopt_adj_noemi
./ex9busopt_adj_noemi -tao_monitor

or if would like to load a file containing the observations, do:
./ex9busopt_adj_noemi -tao_monitor -loadObservations obs-perturbed.bin

or if you'd like to see the covariance parameters do:
./ex9busopt_adj_noemi -outputCov

6. Start matlab with:

MATLABPATH=$PETSC_DIR/share/petsc/matlab matlab -nodesktop

Then in your matlab code use 'PetscBinaryRead' function, for example
A = PetscBinaryRead('out_pert.bin');

7. Run scripts to produce surface plots

python inertia_var_response/generate_experiments.py
python inertia_var_response/run_experiments.py params_in.txt 100 test2_large > test2_large.output &
python inertia_var_response/run_experiments_table1.py inertia_var_response/params_table1.txt 1 table1_var_eps1e-5 > table1_var_eps1e-5.output &

8. run code from command line:

./ex9busopt_adj_noemi -tao_monitor -data_noise 0.01 -load_disturb 7.0 -tfinal 2.0
