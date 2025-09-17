Tutorials

# Standard Symmetric Eigenvalue Problem

This example solves a standard symmetric eigenproblem {math}`Ax= \lambda x`, where A is the matrix resulting from the discretization of the Laplacian operator in 1 dimension by centered finite differences.

```{code}
     | 2 -1  0  0  0  0 |
     |-1  2 -1  0  0  0 |
A =  | 0 -1  2 -1  0  0 |
     | 0  0 -1  2 -1  0 |
     | 0  0  0 -1  2 -1 |
     | 0  0  0  0 -1  2 |
```

## Compiling

Copy the file [ex1.c](https://slepc.upv.es/documentation/current/src/eps/tutorials/ex1.c) [[plain text]](https://slepc.upv.es/documentation/current/src/eps/tutorials/ex1.c) to the working directory and add these lines to the makefile

```{code} make
ex1: ex1.o
	-${CLINKER} -o ex1 ex1.o ${SLEPC_EPS_LIB}
	${RM} ex1.o
```

:::{note}
In the above text, the blank space in the 2nd and 3rd lines is a tab.
:::

Build the executable with the command

```{code} console
$ make ex1
```

:::{note}
Example ex1 is also available in Fortran [ex1f.F90](https://slepc.upv.es/documentation/current/src/eps/tutorials/ex1f.F90) [[plain text]](https://slepc.upv.es/documentation/current/src/eps/tutorials/ex1f.F90).
:::

## Running the Program

In order to run the program for a problem of order 50, type the following

```{code} console
$ ./ex1 -n 50
```

You should get an output similar to this

```{code}
1-D Laplacian Eigenproblem, n=50

 Number of iterations of the method: 8
 Solution method: krylovschur

 Number of requested eigenvalues: 1
 Stopping condition: tol=1e-08, maxit=100
 Number of converged eigenpairs: 3

           k          ||Ax-kx||/||kx||
   ----------------- ------------------
       3.996207        4.30363e-10
       3.984841        2.08631e-09
       3.965946        9.98404e-09
```

## Source Code Details

Examine the source code of the sample program and locate the function calls mentioned in the following comments.

**The Options Database** : All the PETSc functionality related to the options database is available in SLEPc. This allows the user to input control data at run time very easily. In this example, the function [PetscOptionsGetInt](https://petsc.org/release/manualpages/Sys/PetscOptionsGetInt) is used to check whether the user has provided a command line option to set the value of n, the problem dimension. If so, the variable n is set accordingly; otherwise, n remains unchanged.

**Vectors and Matrices** : Usage of matrices and vectors in SLEPc is exactly the same as in PETSc. The user can create a new parallel or sequential matrix, A, with subroutine [MatCreate](https://petsc.org/release/manualpages/Mat/MatCreate), where the matrix format can be specified at runtime. The example creates a matrix, sets the nonzero values with [MatSetValues](https://petsc.org/release/manualpages/Mat/MatSetValues) and then assembles it.

**Solving the Eigenvalue Problem** : Usage of eigensolvers is very similar to other kinds of solvers provided by PETSc. After creating the matrix, the problem is solved by means of an EPS object (Eigenvalue Problem Solver) via the following sequence of function calls:

* [EPSCreate](https://slepc.upv.es/documentation/current/docs/manualpages/EPS/EPSCreate)`(MPI_Comm comm,EPS *eps);`
* [EPSSetOperators](https://slepc.upv.es/documentation/current/docs/manualpages/EPS/EPSSetOperators)`(EPS eps,Mat A,Mat B);`
* [EPSSetProblemType](https://slepc.upv.es/documentation/current/docs/manualpages/EPS/EPSSetProblemType)`(EPS eps,EPSProblemType type);`
* [EPSSetFromOptions](https://slepc.upv.es/documentation/current/docs/manualpages/EPS/EPSSetFromOptions)`(EPS eps);`
* [EPSSolve](https://slepc.upv.es/documentation/current/docs/manualpages/EPS/EPSSolve)`(EPS eps);`
* [EPSGetConverged](https://slepc.upv.es/documentation/current/docs/manualpages/EPS/EPSGetConverged)`(EPS eps, int *nconv);`
* [EPSGetEigenpair](https://slepc.upv.es/documentation/current/docs/manualpages/EPS/EPSGetEigenpair)`(EPS eps,int i,PetscScalar *kr,PetscScalar *ki,Vec xr,Vec xi);`
* [EPSDestroy](https://slepc.upv.es/documentation/current/docs/manualpages/EPS/EPSDestroy)`(EPS eps)`;

First, the eigenproblem solver (EPS) context is created and the operator(s) associated with the eigensystem are set, as well as the problem type. Then various options are set for customized solution. After that, the program solves the problem, retrieves the solution, and finally destroys the EPS context.

The above function calls are very important and will be present in most SLEPc programs. In the example source code [ex1.c](https://slepc.upv.es/documentation/current/src/eps/tutorials/ex1.c) you will find other functions apart from these. What do they do?

## Playing with EPS Options

Now we are going to experiment with different options of the EPS object. A full list of command-line options can be obtained by running the example with the option `-help`.

To show information about the solver object:

```{code} console
$ ./ex1 -eps_view
```

:::{note}
This option internally calls the function [EPSView](https://slepc.upv.es/documentation/current/docs/manualpages/EPS/EPSView).  Alternatively, we could include a direct call to this function in the source code. Almost all command-line options have a related function call.
:::
:::{note}
All the command-line options related to the EPS object have the `-eps_` prefix.
:::

This time, your output will include something like this

```{code}
EPS Object: 1 MPI processes
  type: krylovschur
    50% of basis vectors kept after restart
    using the locking variant
  problem type: symmetric eigenvalue problem
  selected portion of the spectrum: largest eigenvalues in magnitude
  number of eigenvalues (nev): 1
  number of column vectors (ncv): 16
  maximum dimension of projected problem (mpd): 16
  maximum number of iterations: 100
  tolerance: 1e-08
  convergence test: relative to the eigenvalue
BV Object: 1 MPI processes
  type: svec
  17 columns of global length 30
  vector orthogonalization method: classical Gram-Schmidt
  orthogonalization refinement: if needed (eta: 0.7071)
  block orthogonalization method: GS
  doing matmult as a single matrix-matrix product
DS Object: 1 MPI processes
  type: hep
  solving the problem with: Implicit QR method (_steqr)
ST Object: 1 MPI processes
  type: shift
  shift: 0
  number of matrices: 1
```

This option is very useful to see which solver and options the program is using.

Try solving a much larger problem, for instance with n=400. Note that in that case the program does not return a solution. This means that the solver has reached the maximum number of allowed iterations and the convergence criterion was not satisfied. What we can do is either increase the number of iterations

```{code} console
$ ./ex1 -n 400 -eps_max_it 400
```

or relax the convergence criterion.

```{code} console
$ ./ex1 -n 400 -eps_tol 1e-3
```

Note that in the latter case the relative error displayed by the program is significantly larger, meaning that the solution has only 3 correct decimal digits, as expected.

It is possible to change the number of requested eigenvalues. Try the following execution

```{code} console
$ ./ex1 -n 400 -eps_nev 3 -eps_tol 1e-7
```

In this case, the program did not succeed to compute two of the requested eigenpairs. This is again due to the convergence criterion, which is satisfied by some eigenpairs but not for all. As in the previous case, we could increase further the number of iterations or relax the convergence criterion. Another alternative is to increase the number of column vectors (i.e., the dimension of the subspace with which the eigensolver works). This usually improves the convergence behavior at the expense of larger memory requirements.

```{code} console
$ ./ex1 -n 400 -eps_nev 3 -eps_ncv 24
```

Note that the default value of `ncv` depends on the value of `nev`.

Try to set some of the above options directly in the source code by calling the related functions [EPSSetTolerances](https://slepc.upv.es/documentation/current/docs/manualpages/EPS/EPSSetTolerances) and [EPSSetDimensions](https://slepc.upv.es/documentation/current/docs/manualpages/EPS/EPSSetDimensions).  Modify and recompile the program. Use `-eps_view` to check that the values are correctly set. Is it now possible to change these options from the command- line? Does this change whether you place the calls before or after the call to [EPSSetFromOptions](https://slepc.upv.es/documentation/current/docs/manualpages/EPS/EPSSetFromOptions)?

Convergence is usually bad when eigenvalues are close to each other, which is the case in this example. In order to see what is happening while the eigensolver iterates, we can use a monitor to display information associated to the convergence of eigenpairs at each iteration:

```{code} console
$ ./ex1 -eps_monitor
```
or

```{code} console
$ ./ex1 -eps_monitor_all
```
Also, in some SLEPc installations, it is possible to monitor convergence graphically with the `draw` viewer with `draw_lg` format. For example, try this:

```{code} console
$ ./ex1 -n 700 -eps_nev 5 -eps_ncv 35 -eps_monitor_all draw::draw_lg -draw_pause .1
```

:::{note}
The plot is drawn in an X11 pop-up window. So this requires that the display is correctly exported.
:::

## Changing the Eigensolver

The convergence behavior for a particular problem also depends on the properties of the eigensolver being used. SLEPc provides several eigensolvers which can be selected in the source code with the function [EPSSetType](https://slepc.upv.es/documentation/current/docs/manualpages/EPS/EPSSetType), or at run time:

```{code} console
$ ./ex1 -eps_nev 4 -eps_type lobpcg
```

The following table shows some of the solvers available in SLEPc.

Solver                                  |  Command-line Name  |  Parameter
---                                     |  ---                |  ---
Krylov-Schur                            |  krylovschur        |  EPSKRYLOVSCHUR
Generalized Davidson                    |  gd                 |  EPSGD
Jacobi-Davidson                         |  jd                 |  EPSJD
Rayleigh-quotient conjugate gradient    |  rqcg               |  EPSRQCG
Locally optimal block preconditioned CG |  lobpcg             |  EPSLOBPCG
Contour integral spectrum slice         |  ciss               |  EPSCISS
Lanczos with Explicit Restart           |  lanczos            |  EPSLANCZOS
Arnoldi with Explicit Restart           |  arnoldi            |  EPSARNOLDI
Subspace Iteration                      |  subspace           |  EPSSUBSPACE
Power / RQI                             |  power              |  EPSPOWER
Lapack                                  |  lapack             |  EPSLAPACK
ARPACK                                  |  arpack             |  EPSARPACK

:::{note}
The Lapack solver is not really a full-featured eigensolver but simply an interface to some LAPACK routines. These routines operate sequentially in dense mode and therefore are suitable only for small size problems. This solver should be used only for debugging purposes.
:::

:::{note}
The last one (ARPACK) may or may not be available on your system, depending on whether it was enabled during installation of SLEPc. It consists in an interface to the external [ARPACK library](https://github.com/opencollab/arpack-ng). Interfaces to other external libraries may be available as well. These can be used as any other SLEPc native eigensolver.
:::

:::{note}
The default solver is `krylovschur` for both symmetric and non-symmetric problems.
:::
