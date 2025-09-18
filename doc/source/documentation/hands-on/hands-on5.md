Tutorials

# Problem without Explicit Matrix Storage

In many applications, it may be better to keep the matrix (or matrices) that define the eigenvalue problem implicit, that is, without storing its nonzero entries explicitly. An example is when we have a matrix-vector routine available. SLEPc allows easy management of this case. This exercise tries to illustrate it by solving a standard symmetric eigenproblem corresponding to the Laplacian operator in 2 dimensions in which the matrix is not built explicitly.

## Compiling

Copy the file {{'[ex3.c](https://slepc.upv.es/{}/src/eps/tutorials/ex3.c.html)'.format(branch)}} to your directory and add these lines to the makefile

```{code} make
ex3: ex3.o
	-${CLINKER} -o ex3 ex3.o ${SLEPC_EPS_LIB}
	${RM} ex3.o
```

:::{note}
In the above text, the blank space in the 2nd and 3rd lines is a tab.
:::

Build the executable with the command

```{code} console
$ make ex3
```

## Source Code Details

PETSc provides support for matrix-free problems via the _shell_ matrix type.
This kind of matrices is created with a call to [MatCreateShell](https://petsc.org/release/manualpages/Mat/MatCreateShell), and their operations are specified with [MatShellSetOperation](https://petsc.org/release/manualpages/Mat/MatShellSetOperation).
For basic use of these matrices with EPS solvers only the matrix-vector product operation is required. In the example, this operation is performed by a separate function `MatMult_Laplacian2D`.

## Running the Program

Run the program without any spectral transformation options. For instance:

```{code} console
$ ./ex3 -eps_type subspace -eps_tol 1e-9 -eps_nev 8
```

Now try running the program with shift-and-invert to get the eigenvalues closest to the origin

```{code} console
$ ./ex3 -eps_target 0.0 -st_type sinvert
```

Note that the above command yields a run-time error. Observe the information printed on the screen and try to deduce the reason of the error. In this case, the error is due to the fact that SLEPc tries to use a direct linear solver within the ST object, and this is not possible unless the matrix has been created explicitly as in previous examples.

There are more chances to have success if an inexact shift-and-invert scheme is used. Try using an interative linear solver without preconditioning:

```{code} console
$ ./ex3 -eps_target 0.0 -st_type sinvert -st_ksp_rtol 1e-10 -st_ksp_type gmres -st_pc_type none
```

You can also try with a nonzero target.

This example is slow, because the iterative linear solver takes a lot of iterations to reach the required precision (add `-st_ksp_converged_reason` to monitor the convergence of the linear solver). In order to alleviate this problem, a preconditioner should be used. However, a powerful preconditioner such as ILU cannot be used in this case, for the same reason a direct solver is not available. The only possibility is to use a simple preconditioner such as Jacobi. Try running the last example again with `-st_pc_type jacobi`. Note that this works because the `MATOP_GET_DIAGONAL` has been defined in our program. In this particular example the Jacobi preconditioner does not help reduce the linear iterations since the matrix is not strongly diagonally dominant.

Some of the above difficulties can be avoided by using a preconditioned eigensolver, as in the examples shown in previous exercises.
