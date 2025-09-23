Tutorials

# Standard Non-Symmetric Eigenvalue Problem

In this exercise we are going to work with a non-symmetric problem. The example solves the eigenvalue problem associated with a Markov model of a random walk on a triangular grid. Although the matrix is non-symmetric, all eigenvalues are real. Eigenvalues come in pairs with the same magnitude and different signs. The values 1 and -1 are eigenvalues for any matrix size. More details about this problem can be found at [Matrix Market](https://math.nist.gov/MatrixMarket/data/NEP/mvmrwk/mvmrwk).

## Compiling

Copy the file {{'[ex5.c](https://slepc.upv.es/{}/src/eps/tutorials/ex5.c.html)'.format(branch)}} to your directory and add these lines to the makefile

```{code} make
ex5: ex5.o
	-${CLINKER} -o ex5 ex5.o ${SLEPC_EPS_LIB}
	${RM} ex5.o
```

:::{note}
In the above text, the blank space in the 2nd and 3rd lines is a tab.
:::

Build the executable with the command

```{code} console
$ make ex5
```

## Source Code Details

The example program is very similar to that in Exercise 1. The main difference is that the problem is set to be non-symmetric with `EPSSetProblemType`:

```{code} c
PetscCall(EPSSetProblemType(eps,EPS_NHEP));
```

In this example we also illustrate the use of `EPSSetInitialSpace`.

## Running the Program

Run the program requesting four eigenvalues.

```{code} console
$ ./ex5 -eps_nev 4
```

The output will look like this:

```{code}
Markov Model, N=120 (m=15)

 Solution method: krylovschur

 Number of requested eigenvalues: 4
 Linear eigensolve converged (4 eigenpairs) due to CONVERGED_TOL; iterations 5
 ---------------------- --------------------
            k             ||Ax-kx||/||kx||
 ---------------------- --------------------
       -1.000000            7.05086e-10
        1.000000            1.47517e-09
       -0.971367            1.04527e-10
        0.971367            2.23781e-10
 ---------------------- --------------------
```

You can see that the solver returns both positive and negative eigenvalues.  This is because largest magnitude eigenvalues are computed by default, that is, internally the solver sorts the eigenvalue approximations according to {math}`|\lambda|`, and the same criterion is used for sorting the finally computed eigenvalues.

Other criteria can be used, see `EPSSetWhichEigenpairs` for details. For instance, for computing only the rightmost eigenvalues, try the following.

```{code} console
$ ./ex5 -eps_nev 4 -eps_largest_real
```

Similarly, it is possible to request the smallest magnitude eigenvalues with `-eps_smallest_magnitude`. The difference in that case is that the solver needs much more iterations to converge. The justification is that in this problem the smallest magnitude eigenvalues are located in the interior of the spectrum, and computing interior eigenvalues is always harder as explained next.

## Computing Interior Eigenvalues

It is well known that computing eigenvalues located at the interior of the spectrum is much more difficult than those in the periphery. We are going to discuss different strategies available in SLEPc.  The general way of computing interior eigenvalues is to specify a target value, around which the eigenvalues must be sought.

```{code} console
$ ./ex5 -eps_nev 4 -eps_target 0.75
```

Note that apart from the target value, one should specify a sorting criterion relative to the target (`-eps_target_magnitude`). However, this option can be omitted because it is the default when a target is specified. The output is in this case:

```{code}
Markov Model, N=120 (m=15)

 Solution method: krylovschur

 Number of requested eigenvalues: 4
 Linear eigensolve converged (6 eigenpairs) due to CONVERGED_TOL; iterations 14
 ---------------------- --------------------
            k             ||Ax-kx||/||kx||
 ---------------------- --------------------
        0.771298            1.13746e-09
        0.714286            9.26568e-10
        0.704038            2.49204e-09
        0.702317            8.51587e-10
        0.800066            2.80221e-14
        0.660191            7.61205e-09
        0.840634            1.94667e-14
        0.857143            1.67828e-10
 ---------------------- --------------------
```

We have obtained eigenvalues both on the left and on the right of the target value {math}`\tau=0.75`, and they are sorted according to the distance to {math}`\tau`.

The number of iterations is higher than in the default case. The theory says that Krylov methods (and other methods as well) approximate eigenvalues from the periphery to the interior, meaning that before getting eigenvalues closest to 0.75 the solver has to find out the eigenvalues from 0.75 to the rightmost extreme. If we choose a target close to the extreme the number of iterations will be small, and they will increase as {math}`\tau` is moved inside of the spectrum. Therefore, this is not a good strategy because it will not be viable for difficult problems.

Sometimes, an improvement may come from changing the way in which the method extracts the spectral information from the built subspace; see `EPSSetExtraction` for details. One such technique is called harmonic extraction. Try the following:

```{code} console
$ ./ex5 -eps_nev 4 -eps_target 0.75 -eps_harmonic
```

In this simple problem, harmonic extraction gives no benefit but in difficult problems it may be a significant improvement, especially in combination with preconditioned solvers (discussed later below).
A better solution may be to use a spectral transformation, but with several considerations to take into account regarding cost.

## Getting Started with Spectral Transformations

The general idea of the spectral transformation is to substitute the original problem, {math}`Ax=\lambda x`, by another one, {math}`Tx= \theta x`, in which the eigenvalues are mapped to a different position but eigenvectors remain unchanged. With this strategy, one can move interior eigenvalues to the periphery.

Each EPS object uses an ST object internally to manage the spectral transformation. The following table shows the available spectral transformations, which can be selected with the function `STSetType` or at run time.

Spectral Transformation  |  Operator                                   |  Command-line Name  |  Parameter
---                      |  ---                                        |  ---                |  ---
Shift of Origin          |  {math}`A - \sigma I`                       |  shift              |  STSHIFT
Shift-and-invert         |  {math}`(A- \sigma I)^{-1}`                 |  sinvert            |  STSINVERT
Cayley                   |  {math}`(A- \sigma I)^{-1} (A + \nu I)`     |  cayley             |  STCAYLEY
Preconditioner           |  {math}`K^{-1} \approx (A- \sigma I)^{-1}`  |  precond            |  STPRECOND
Polynomial filter        |  {math}`p(A)`                               |  filter             |  STFILTER

:::{note}
The default is to do shift of origin with a value {math}`\sigma =0`. This was reported by `-eps_view` in the previous example.
:::
:::{note}
The preconditioner is not really a spectral transformation like the rest. It will be discussed later below.
:::

The shift-and-invert spectral transformation can be used for computing interior eigenvalues:

```{code} console
$ ./ex5 -eps_nev 4 -eps_target 0.75 -st_type sinvert
```

:::{note}
All the command-line options related to the ST object have the `-st_` prefix.
:::

With the above execution, the number of iterations is very small, but each iteration is much more costly than in the previous cases because linear systems must be solved to handle the inverted operator (the issue of how to solve linear systems is discussed below). The value of the parameter {math}`\sigma` (the shift) is taken to be equal to {math}`\tau` (the target). Run with `-eps_view` to check that it is indeed the case.

Try also with `cayley`, which is nearly equivalent.

## Handling the Inverses

In the table of spectral transformations shown above, there are some operators that include the inverse of a certain matrix. These operators are not computed explicitly in order to preserve sparsity. Instead, in the ST object the multiplication by these inverses is replaced by a linear equation solve via a KSP object from PETSc.

SLEPc allows us to pass options to this KSP linear solver object. For instance,

```{code} console
$ ./ex5 -eps_nev 4 -eps_target 0.75 -st_type sinvert -st_ksp_type preonly -st_pc_type lu
```

:::{note}
In order to specify a command-line option related to the linear solver contained in ST, simply add the `-st_` prefix in front.
:::

The options of the above example specify a direct linear solver (LU factorization). This is what SLEPc does by default. This strategy is usually called _exact shift-and-invert_. Its main drawback is that direct solvers are more costly in terms of flops and storage and are less parallelizable.

An alternative is to do an _inexact shift-and-invert_ , that is, to use an iterative linear solver. The following line illustrates how to use an iterative solver

```{code} console
$ ./ex5 -eps_nev 4 -eps_target 0.75 -st_type sinvert -st_ksp_type gmres -st_pc_type bjacobi -st_ksp_rtol 1e-12
```

Iterative linear solvers may fail to converge if the coefficient matrix is ill-conditioned or close to singular. Also, the accuracy of the eigensolver may be compromised if the iterative linear solver provides a solution far from full working precision.

Note that in SLEPc it is extremely easy to switch between exact and inexact schemes.

## Preconditioned Eigensolvers

As mentioned above, the inexact shift-and-invert scheme is very sensitive to the accuracy with which the linear systems are solved. This usually implies using a very stringent tolerance (10-12 in the example) and makes it impractical for difficult problems.

An alternative is to use a preconditioned eigensolver, such as those of Davidson type: `EPSGD` and `EPSJD`. These solvers try to emulate the idea of shift-and-invert but they are very robust with respect to bad accuracy (i.e., large tolerance) of the iterative linear solve.

Here is an example where Jacobi-Davidson is used:

```{code} console
$ ./ex5 -eps_nev 4 -eps_target 0.75 -eps_type jd -st_type precond -st_ksp_type bcgsl -st_pc_type bjacobi -st_ksp_rtol 0.001
```

:::{note}
The `-st_type precond` key can be omitted in this case, since it is the default in all preconditioned eigensolvers.
:::

Try adding `-eps_harmonic` to the above example. As mentioned before, harmonic extraction is usually better when used in preconditioned solvers.
