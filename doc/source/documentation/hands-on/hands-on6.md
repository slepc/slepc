Tutorials

# Parallel Execution

The objective of this exercise is to run an example program with different number of processors to see how execution time is reduced. This time, we are going to solve a standard eigensystem {math}`Ax= \lambda x` with the matrix loaded from a file. In particular, the matrix we are going to use is QC2534. It is a complex matrix of order 2534 arising from a quantum chemistry application (more details can be found at [Matrix Market](https://math.nist.gov/MatrixMarket/data/NEP/h2plus/h2plus)).

## Compiling
c
Copy the file {{'[ex4.c](https://slepc.upv.es/{}/src/eps/tutorials/ex4.c.html)'.format(branch)}} to your directory and add these lines to the makefile

```{code} make
ex4: ex4.o
	-${CLINKER} -o ex4 ex4.o ${SLEPC_EPS_LIB}
	${RM} ex4.o
```

:::{note}
In the above text, the blank space in the 2nd and 3rd lines is a tab.
:::

Build the executable with the command (optimized complex version)

```{code} console
$ make PETSC_ARCH=arch-linux-gnu-c-opt-complex ex4
```

## Source Code Details

This example program is very similar to that of exercise 3. It uses the PETSc function {{'[MatLoad](https://petsc.org/{}/manualpages/Mat/MatLoad)'.format(branch)}} to load a matrix from a file. The matrix file is specified in the command line.

## Running the Program

In order to run this example, you will need the file `qc2534.petsc`. Locate it in the file system and then run the program with the command

```{code} console
$ ./ex4 -file qc2534.petsc
```

For execution with more than one process:

```{code} console
$ mpiexec -n 2 ./ex4 -file qc2534.petsc
```

Check the output of the program. It should be the same as with one process.

Try using the `-log_view` option to have a look at the profiling information collected by PETSc and SLEPc. For instance, check the size and number of MPI messages.

```{code} console
$ mpiexec -n 2 ./ex4 -file qc2534.petsc -log_view
```

Try to find out how much time was spent is solving the eigenvalue problem. Is there a significant reduction when we increase the number of processors?

## Instrumenting the Source Code

If we are just interested in knowing the time used by the eigensolver, then it may be better to let our example program inform us. With the function {{'[PetscTime](https://petsc.org/{}/manualpages/Sys/PetscTime)'.format(branch)}}, it is possible to obtain the current time of day (wall-clock time) in seconds. Edit the source code and add two calls to this function just before and after the `EPSSolve` call, as in the following fragment of code

```{code} c
PetscCall(PetscTime(&t1;));
PetscCall(EPSSolve(eps));
PetscCall(PetscTime(&t2;));
PetscCall(PetscPrintf(PETSC_COMM_WORLD," Elapsed Time: %f\n",t2-t1));
```

Also you must add the definition of the two new variables

```{code} c
PetscLogDouble t1,t2;
```

as well as the header

```{code} c
#include <petsctime.h>
```

Recompile the program with

```{code} console
$ make PETSC_ARCH=arch-linux-gnu-c-opt-complex ex4
```

Run it with one, two and four processors checking the time spent by the solver

```{code} console
$ mpiexec -n 1 ./ex4 -file qc2534.petsc
$ mpiexec -n 2 ./ex4 -file qc2534.petsc
$ mpiexec -n 4 ./ex4 -file qc2534.petsc
```
