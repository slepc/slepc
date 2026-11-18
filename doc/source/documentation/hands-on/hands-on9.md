# Exercise 9: Structured Eigenvalue Problem

The goal of this exercise is to illustrate how to solve structured eigenvalue problems. In particular, we are going to solve a problem with [](#sec:bse) structure (BSE). Due to its algebraic properties, all eigenvalues are real and come in pairs {math}`\{\lambda,-\lambda\}`. The particular matrix generated in the example is fictitious, but is helpful to mimic problems that one can find in materials science applications.

We recall that in this problem the matrix has the following block structure,

```{math}
H = \begin{bmatrix}
        R & C \\
        -C^* & -R^T
    \end{bmatrix},
```

where $R$ is Hermitian and $C$ is complex symmetric.

:::{warning}
The example is intended to be used with a PETSc configuration with complex scalars. If you use it with real scalars, the problem to be solved will be different, and instead of the BSE structure it will be a [](#sec:lrep) eigenvalue problem (LREP).
:::

## Compiling

Copy the file {{'[ex55.c](https://slepc.upv.es/{}/src/eps/tutorials/ex55.c.html)'.format(branch)}} to your directory and add these lines to the makefile

```{code} make
ex55: ex55.o
	-${CLINKER} -o ex55 ex55.o ${SLEPC_EPS_LIB}
	${RM} ex55.o
```

Build the executable with the command

```{code} console
$ make ex55
```

## Source Code Details

The organization of the example program is similar to that of previous exercises. Note that the program is prepared to read two matrices {math}`R` and {math}`C` from a file with the `-R` and `-C` command-line options. If not provided, the program will generate simple {math}`R` and {math}`C` that are pentadiagonal and tridiagonal, respectively, and satisfy the definition of the BSE problem.

The main difference with respect to other examples is that we have to use a special matrix construction function `MatCreateBSE()` to obtain a Hamiltonian matrix {math}`H` with the wanted structure:

```{code} c
PetscCall(MatCreateBSE(R,C,&H));
```

The matrix {math}`H` is a 2-by-2 block matrix in PETSc nested format (`MATNEST`) with a special mark so that SLEPc knows it has a BSE structure.

Note that we have to configure the `EPS` object with `EPSSetProblemType`(`eps`,`EPS_BSE`) so that the solver takes into account the BSE structure, otherwise it will solve the problem as non-Hermitian.

## Running the Program

Run the program with problem size 100 requesting two eigenvalues.

```{code} console
$ ./ex55 -n 100 -eps_nev 2
```

The output will look like this:

```{code}
Bethe-Salpeter eigenproblem, n=100

 Linear eigensolve converged (2 eigenpairs) due to CONVERGED_TOL; iterations 22
 ---------------------- --------------------
            k             ||Ax-kx||/||kx||
 ---------------------- --------------------
        2.151876            1.84634e-09
       -2.151876            1.84634e-09
 ---------------------- --------------------
```

You can see that the solver returns both positive and negative eigenvalues. These eigenvalues are those of smallest magnitude. You can check this by increasing the number of wanted eigenvalues, e.g., `-eps_nev 12`. The number of iterations is smaller in this case, the reason is that the size of the subspace is larger. You can see this with the `-eps_view` option. Use `-eps_ncv` to change the subspace size if you wish.

One thing you can try is solving the problem as non-Hermitian:

```{code} console
$ ./ex55 -n 100 -eps_nev 12 -eps_non_hermitian -eps_smallest_magnitude
```

The output should be similar. However, the computation takes longer, and furthermore you can see that the printed eigenvalues have a zero imaginary part, which indicates that the solver is producing a small error because it is disregarding the BSE structure.

The structured solver can also compute the largest magnitude eigenvalues with `-eps_largest_magnitude`, although this is generally not needed in applications (that is why the default in BSE is `-eps_smallest_magnitude`).

The example has an additional command-line option `-checkorthog` that will obtain both the right and left eigenvectors for each computed eigenvalue and check that they are bi-orthogonal.

## Using a Threshold instead of Specifying the Eigenvalue Count

The previous runs explicitly indicate how many eigenvalues we want, with `-eps_nev`. However, sometimes we want to compute eigenvalues that are below a given threshold, but we do not know how many satisfy this restriction. SLEPc provides an option for this, that can be used in BSE problems as well as in other contexts. In BSE problems, the meaning of the threshold parameter is compute all eigenvalues whose absolute value is below the given threshold. See `EPSSetThreshold()` for additional details.

```{code} console
$ ./ex55 -n 100 -eps_threshold_absolute 2.3
```

The above execution will return 18 eigenvalues. Note that we can also specify a relative threshold instead of absolute.
