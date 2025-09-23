Tutorials

# Use of Deflation Subspaces

The term deflation refers to the use of the knowledge of one or more eigenpairs to find other eigenpairs. For instance, most eigensolvers try to approximate a number of eigenpairs and, as soon as one of them has converged, they deflate it for better approximating the other ones. Another case is when one eigenpair is known a priori and one wants to use this knowledge to compute other eigenpairs. SLEPc supports this by means of _deflation subspaces_.

This example illustrates the use of deflation subspaces to compute the smallest nonzero eigenvalue of the Laplacian of a graph corresponding to a 2-D regular mesh. The problem is a standard symmetric eigenproblem {math}`Ax= \lambda x`, where {math}`A = L(G)` is the Laplacian of graph G, defined as follows: Aii = degree of node i, Aij = -1 if edge (i,j) exists in G, zero otherwise. This matrix is symmetric positive semidefinite and singular, and [1 1 ... 1]T is the eigenvector associated with the zero eigenvalue. In graph theory, one is usually interested in computing the eigenvector associated with the next eigenvalue (the so-called Fiedler vector).

## Compiling

Copy the file {{'[ex11.c](https://slepc.upv.es/{}/src/eps/tutorials/ex11.c.html)'.format(branch)}} to your directory and add these lines to the makefile

```{code} make
ex11: ex11.o
	-${CLINKER} -o ex11 ex11.o ${SLEPC_EPS_LIB}
	${RM} ex11.o
```

:::{note}
In the above text, the blank space in the 2nd and 3rd lines is a tab.
:::

Build the executable with the command

```{code} console
$ make ex11
```

## Source Code Details

This example computes the smallest eigenvalue by setting `EPS_SMALLEST_REAL` in `EPSSetWhichEigenpairs`. An alternative would be to use a shift-and-invert spectral transformation with a zero target to compute the eigenvalues closest to the origin, or to use harmonic extraction with a zero target.

By specifying a deflation subspace (the one associated to the eigenvector [1 1 ... 1]T) with the function `EPSSetDeflationSpace`, the convergence to the zero eigenvalue is avoided. Thus, the program should compute the smallest nonzero eigenvalues.

## Running the Program

Run the program simply with

```{code} console
$ ./ex11
```

For the case of using an inexact spectral transformation, the command line would be:

```{code} console
$ ./ex11 -eps_target 0.0 -eps_target_real -st_type sinvert -st_ksp_rtol 1e-10 -st_ksp_type gmres -st_pc_type jacobi
```

Note that a shift-and-invert spectral transformation should always be used in combination with `EPS_TARGET_MAGNITUDE` or `EPS_TARGET_REAL`.

And for the case of harmonic extraction:

```{code} console
$ ./ex11 -eps_target 0.0 -eps_target_real -eps_harmonic
```
