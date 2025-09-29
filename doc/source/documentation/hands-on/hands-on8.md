# Exercise 8: Quadratic Eigenvalue Problem

Now we are going to focus on the solution of quadratic eigenvalue problems with `PEP` solvers. In this case, the problem to be solved is formulated as {math}`(\lambda^{2}M+ \lambda C+K)x=0`. In our simple example, {math}`M` is a diagonal matrix, {math}`C` is tridiagonal, and {math}`K` is the 2-D Laplacian.

## Compiling

Copy the file {{'[ex16.c](https://slepc.upv.es/{}/src/pep/tutorials/ex16.c.html)'.format(branch)}} to your directory and add these lines to the makefile

```{code} make
ex16: ex16.o
	-${CLINKER} -o ex16 ex16.o ${SLEPC_PEP_LIB}
	${RM} ex16.o
```

Build the executable with the command

```{code} console
$ make ex16
```

## Running the Program

Run the program without arguments and see the output:

```{code}
Quadratic Eigenproblem, N=100 (10x10 grid)

 Number of requested eigenvalues: 1

           k          ||P(k)x||/||kx||
   ----------------- ------------------
 -1.164037+1.653625i    2.02021e-12
 -1.164037-1.653625i    2.02021e-12
```

## Source Code Details

The `PEP` object is used very much like `EPS` or `SVD`, as can be seen in the source code. Here is a summary of the main function calls:

* `PEPCreate(MPI_Comm comm,PEP *pep);`
* `PEPSetOperators`(`PEP` pep,{external:doc}`PetscInt` `nmat`,{external:doc}`Mat` `A[]`);
* `PEPSetProblemType(PEP pep,PEPProblemType type);`
* `PEPSetFromOptions(PEP pep);`
* `PEPSolve(PEP pep);`
* `PEPGetConverged`(`PEP` `pep`, {external:doc}`PetscInt` `*nconv`);
* `PEPGetEigenpair`(`PEP` `pep`,{external:doc}`PetscInt` `i`,{external:doc}`PetscScalar` `*kr`,{external:doc}`PetscScalar` `*ki`,{external:doc}`Vec` `xr`,{external:doc}`Vec` `xi`);
* `PEPDestroy(PEP pep)`;

First, the solver context (`PEP`) is created and the three problem matrices are specified. Then various options are set for customized solution. After that, the program solves the problem, retrieves the solution, and finally destroys the context.

## PEP Options

Most of the options available in the `EPS` object have their equivalent in `PEP`.  A full list of command-line options can be obtained by running the example
with the option `-help`.

To show information about the `PEP` solver, add the `-pep_view` option:

```{code}
PEP Object: 1 MPI processes
  type: toar
    50% of basis vectors kept after restart
    using the locking variant
  problem type: symmetric polynomial eigenvalue problem
  polynomial represented in MONOMIAL basis
  selected portion of the spectrum: largest eigenvalues in magnitude
  number of eigenvalues (nev): 1
  number of column vectors (ncv): 16
  maximum dimension of projected problem (mpd): 16
  maximum number of iterations: 100
  tolerance: 1e-08
  convergence test: relative to the eigenvalue
  extraction type: NORM
BV Object: 1 MPI processes
  type: svec
  18 columns of global length 100
  vector orthogonalization method: classical Gram-Schmidt
  orthogonalization refinement: if needed (eta: 0.7071)
  block orthogonalization method: GS
  doing matmult as a single matrix-matrix product
DS Object: 1 MPI processes
  type: nhep
ST Object: 1 MPI processes
  type: shift
  shift: 0.
  number of matrices: 3
  all matrices have different nonzero pattern
  KSP Object: (st_) 1 MPI processes
    type: preonly
    maximum iterations=10000, initial guess is zero
    tolerances:  relative=1e-08, absolute=1e-50, divergence=10000.
    left preconditioning
    using NONE norm type for convergence test
  PC Object: (st_) 1 MPI processes
    type: lu
      out-of-place factorization
      tolerance for zero pivot 2.22045e-14
      matrix ordering: nd
      factor fill ratio given 5., needed 1.
        Factored matrix follows:
          Mat Object: 1 MPI processes
            type: seqaij
            rows=100, cols=100
            package used to perform factorization: petsc
            total: nonzeros=100, allocated nonzeros=100
            total number of mallocs used during MatSetValues calls =0
              not using I-node routines
    linear system matrix = precond matrix:
    Mat Object: 1 MPI processes
      type: seqaij
      rows=100, cols=100
      total: nonzeros=100, allocated nonzeros=500
      total number of mallocs used during MatSetValues calls =0
        not using I-node routines
```

:::{note}
All the command-line options related to the `PEP` object have the `-pep_` prefix.
:::

Try changing some of the values, for example:

```{code} console
$ ./ex16 -pep_nev 4 -pep_ncv 24 -pep_smallest_magnitude -pep_tol 1e-5
```

## Choosing the Solver Method

Several polynomial eigensolvers are available, which can be selected in the source code with the function `PEPSetType`, or at run time:

```{code} console
$ ./ex16 -pep_type qarnoldi
```

The following table shows the list of `PEP` solvers available in SLEPc.

Solver                        |  Command-line Name  |  Parameter
---                           |---                  |---
Linearization                 |  linear             |  PEPLINEAR
Quadratic Arnoldi             |  qarnoldi           |  PEPQARNOLDI
Two-level orthogonal Arnoldi  |  toar               |  PEPTOAR
Symmetric TOAR                |  stoar              |  PEPSTOAR
Jacobi-Davidson               |  jd                 |  PEPJD
Contour integral              |  ciss               |  PEPCISS

:::{note}
The default solver is `toar`.
:::

The `linear` solver performs an explicit linearization of the quadratic eigenproblem, resulting in a generalized eigenproblem. This linearization can be customized, see section [](#sec:linearization) in the Users Manual for details. For instance:

```{code} console
$ ./ex16 -pep_type linear -pep_linear_linearization 1,0 -pep_linear_explicitmatrix
```

Since in this problem all matrices are symmetric, the problem type is set to `PEP_HERMITIAN` in the source code with `PEPSetProblemType`, and this obliges us to set the explicit matrix flag, see `PEPLinearSetExplicitMatrix`.  It is also possible to use a non-symmetric linearization by choosing the corresponding problem type:

```{code} console
$ ./ex16 -pep_type linear -pep_general
```

In the `linear` solver it is also possible to tune any of the `EPS` options, including those corresponding to `ST` and the linear solvers. For instance:

```{code} console
$ ./ex16 -pep_type linear -pep_general -pep_linear_st_ksp_type bcgs -pep_linear_st_pc_type jacobi
```
