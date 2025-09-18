Tutorials

# Generalized Eigenvalue Problem Stored in a File

In this exercise we are going to work with a generalized eigenvalue problem, {math}`Ax=\lambda Bx`. The example program loads two matrices A and B from a file and then solves the associated eigensystem.

The matrices we are going to work with are BFW62A and BFW62B, which are available at [Matrix Market](https://math.nist.gov/MatrixMarket/data/NEP/bfwave/bfwave). This particular problem is non-symmetric. Some of the eigenvalues (those of largest magnitude) come in complex conjugate pairs while the rest are real.

## Compiling

Copy the file {{'[ex7.c](https://slepc.upv.es/{}/src/eps/tutorials/ex7.c.html)'.format(branch)}} to your directory and add these lines to the makefile

```{code} make
ex7: ex7.o
	-${CLINKER} -o ex7 ex7.o ${SLEPC_EPS_LIB}
	${RM} ex7.o
```

:::{note}
In the above text, the blank space in the 2nd and 3rd lines is a tab.
:::

Build the executable with the command

```{code} console
$ make ex7
```

## Source Code Details

This example uses the PETSc function [MatLoad](https://petsc.org/release/manualpages/Mat/MatLoad) to load a matrix from a file. The two matrix files are specified in the command line.  Note that these files have been converted from Matrix Market format to PETSc binary format.

Compare the source code of the example program with the previous ones. Note that, in this case, two matrix objects are passed in the [EPSSetOperators](../../manualpages/EPS/EPSSetOperators) function call:

```{code} c
PetscCall(EPSSetOperators(eps,A,B));
```

## Running the Program

Run the program with the following command line:

```{code} console
$ ./ex7 -f1 ${SLEPC_DIR}/share/slepc/datafiles/matrices/bfw62a.petsc -f2 ${SLEPC_DIR}/share/slepc/datafiles/matrices/bfw62b.petsc
```

Run the program to compute more than one eigenpair. Use the following option to plot the computed eigenvalues:

```{code} console
$ ./ex7 -f1 ${SLEPC_DIR}/share/slepc/datafiles/matrices/bfw62a.petsc -f2 ${SLEPC_DIR}/share/slepc/datafiles/matrices/bfw62b.petsc -eps_type subspace -eps_nev 6 -eps_view_values draw -draw_pause -1
```

:::{note}
The plot is drawn in an X11 pop-up window. So this requires that the display is correctly exported.
:::

## Spectral Transformations in Generalized Problems

The following table shows the expressions of the operator in each of the available spectral transformations in the case of generalized problems. Note that both matrices A and B are involved.

Spectral Transformation  |  Operator
---                      |  ---
Shift of Origin          |  {math}`B^{-1} A + \sigma I`
Shift-and-invert         |  {math}`(A - \sigma B)^{-1} B`
Cayley                   |  {math}`(A- \sigma B)^{-1} (A+ \nu B)`
Preconditioner           |  {math}`K^{-1} \approx (A - \sigma B)^{-1}`

In the case of generalized problems, the shift-and-invert transformation does not represent a cost penalty with respect to the simpler shift of origin, since in both cases the inverse of a matrix is required.

```{code} console
$ ./ex7 -f1 ${SLEPC_DIR}/share/slepc/datafiles/matrices/bfw62a.petsc -f2 ${SLEPC_DIR}/share/slepc/datafiles/matrices/bfw62b.petsc -eps_target 0 -st_type sinvert
```

The above execution computes the eigenvalues closest to the origin. Use a target near the left end of the spectrum to compute the largest magnitude eigenvalues

```{code} console
$ ./ex7 -f1 ${SLEPC_DIR}/share/slepc/datafiles/matrices/bfw62a.petsc -f2 ${SLEPC_DIR}/share/slepc/datafiles/matrices/bfw62b.petsc -eps_target -250000 -st_type sinvert
```

## Preconditioned Eigensolvers

As hinted above, generalized eigenproblems have the drawback that in the default mode (shift of origin) one of the matrices have to be (implicitly) inverted. However, preconditioned eigensolvers do not have this limitation, and may be able to solve the problem with just a preconditioner or a few iterations of an iterative linear solver.

Here is an example with Generalized Davidson:

```{code} console
$ ./ex7 -f1 ${SLEPC_DIR}/share/slepc/datafiles/matrices/bfw62a.petsc -f2 ${SLEPC_DIR}/share/slepc/datafiles/matrices/bfw62b.petsc -eps_type gd -eps_nev 6
```

Try the above example also with `-eps_target` and `-eps_harmonic`.
