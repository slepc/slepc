#
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  SLEPc - Scalable Library for Eigenvalue Problem Computations
#  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
#
#  This file is part of SLEPc.
#  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#

MANSEC     = EPS
TESTS      = test1 test2 test3 test4 test5 test6 test7f test8 test9 test10 test11 test12 test13 test14 test14f test15f test16 test17 test17f test18 test19 test20 test21 test22 test23 test24 test25 test26 test27 test28 test29 test30 test31 test32 test34 test35 test36 test37 test38 test39 test40 test41 test42 test43 test44 test45

include ${SLEPC_DIR}/lib/slepc/conf/slepc_common

#------------------------------------------------------------------------------------

testtest5_blopex: test5.PETSc
	@test=test5_2_blopex; check=test5_2; \
        if [ "${PETSC_SCALAR}" = "complex" ]; then check=$${check}_complex; fi; \
        ${MPIEXEC} -n 1 ./test5 -symm -eps_type blopex -eps_nev 4 > $${test}.tmp 2>&1; \
        if (${DIFF} output/$${check}.out $${test}.tmp > /dev/null 2>&1) then \
          echo "BLOPEX example src/eps/tests/test5 run successfully with 1 MPI process"; \
        else \
          echo "Possible error running BLOPEX src/eps/tests/test5 with 1 MPI process"; \
          cat $${test}.tmp; \
          touch ../../../check_error; \
        fi; \
        ${RM} $${test}.tmp

testtest7f: test7f.PETSc
	@test=test7f_1; check=test7f_1; \
	GFORTRAN_UNBUFFERED_ALL=1 ${MPIEXEC} -n 1 ./test7f -eps_nev 4 | ${SED} -e 's/83791/83792/'> $${test}.tmp 2>&1; \
        if (${DIFF} output/$${check}.out $${test}.tmp > /dev/null 2>&1) then \
          echo "Fortran example src/eps/tests/test7f run successfully with 1 MPI process"; \
        else \
          echo "Possible error running Fortran src/eps/tests/test7f with 1 MPI process"; \
          cat $${test}.tmp; \
          touch ../../../check_error; \
        fi; \
        ${RM} $${test}.tmp

testtest10: test10.PETSc
	@test=test10_1; check=test10_1_ks; \
        options="-eps_nev 4 -eps_ncv 14 -m 11 -eps_largest_magnitude"; \
        ${MPIEXEC} -n 1 ./test10 $${options} | ${SED} -e 's/82109/82110/' | ${SED} -e "s/58463/58462/" > $${test}.tmp 2>&1; \
        if (${DIFF} output/$${check}.out $${test}.tmp > /dev/null 2>&1) then \
          echo "C/C++ example src/eps/tests/test10 run successfully with 1 MPI process"; \
        else \
          echo "Possible error running C/C++ src/eps/tests/test10 with 1 MPI process"; \
          cat $${test}.tmp; \
          touch ../../../check_error; \
        fi; \
        ${RM} $${test}.tmp

testtest10_mpi: test10.PETSc
	@test=test10_1; check=test10_1_ks; \
        options="-eps_nev 4 -eps_ncv 14 -m 11 -eps_largest_magnitude"; \
        ${MPIEXEC} -n 2 ./test10 $${options} | ${SED} -e 's/82109/82110/' | ${SED} -e "s/58463/58462/" > $${test}.tmp 2>&1; \
        if (${DIFF} output/$${check}.out $${test}.tmp > /dev/null 2>&1) then \
          echo "C/C++ example src/eps/tests/test10 run successfully with 2 MPI processes"; \
        else \
          echo "Possible error running C/C++ src/eps/tests/test10 with 2 MPI processes"; \
          cat $${test}.tmp; \
          touch ../../../check_error; \
        fi; \
        ${RM} $${test}.tmp

testtest10_cuda: test10.PETSc
	@test=test10_1; check=test10_1_ks_cuda; \
        options="-eps_nev 4 -eps_ncv 14 -m 12 -eps_largest_magnitude -mat_type aijcusparse"; \
        ${MPIEXEC} -n 1 ./test10 $${options} | ${SED} -e 's/54988/54989/' | ${SED} -e "s/35009/35008/" > $${test}.tmp 2>&1; \
        if (${DIFF} output/$${check}.out $${test}.tmp > /dev/null 2>&1) then \
          echo "C/C++ example src/eps/tests/test10 run successfully with CUDA"; \
        else \
          echo "Possible error running C/C++ src/eps/tests/test10 with CUDA"; \
          cat $${test}.tmp; \
          touch ../../../check_error; \
        fi; \
        ${RM} $${test}.tmp;

testtest10_hip: test10.PETSc
	@test=test10_1; check=test10_1_ks_cuda; \
        options="-eps_nev 4 -eps_ncv 14 -m 12 -eps_largest_magnitude -mat_type aijhipsparse"; \
        ${MPIEXEC} -n 1 ./test10 $${options} | ${SED} -e 's/54988/54989/' | ${SED} -e "s/35009/35008/" > $${test}.tmp 2>&1; \
        if (${DIFF} output/$${check}.out $${test}.tmp > /dev/null 2>&1) then \
          echo "C/C++ example src/eps/tests/test10 run successfully with HIP"; \
        else \
          echo "Possible error running C/C++ src/eps/tests/test10 with HIP"; \
          cat $${test}.tmp; \
          touch ../../../check_error; \
        fi; \
        ${RM} $${test}.tmp;
