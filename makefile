CFLAGS	        = -g
FFLAGS	        =
CPPFLAGS        =
FPPFLAGS        =
LOCDIR          = ./paramest_9bus
EXAMPLESC       = 
EXAMPLESF       =
EXAMPLESFH      =
MANSEC          = TS
DIRS            =0000

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

ex9businertiaest_adj: ex9businertiaest_adj.o common-inpest.o chkopts
	-${CLINKER} -o ex9businertiaest_adj common-inpest.o ex9businertiaest_adj.o  ${PETSC_LIB}

driver_pce.exe: driver_pce.o common-inpest.o chkopts
	-${CLINKER} -o driver_pce.exe common-inpest.o driver_pce.o  ${PETSC_LIB}

clean_files:
	${RM} ex9businertiaest_adj.o driver_pce.o common-inpest.o

include ${PETSC_DIR}/lib/petsc/conf/test
