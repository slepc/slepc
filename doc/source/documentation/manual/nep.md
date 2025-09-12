(cap:nep)=
# NEP: Nonlinear Eigenvalue Problems

The Nonlinear Eigenvalue Problem (`NEP`) solver object covers the general case where the eigenproblem is nonlinear with respect to the eigenvalue, but it cannot be expressed in terms of a polynomial. We will write the problem as $T(\lambda)x=0$, where $T$ is a matrix-valued function of the eigenvalue $\lambda$. Note that `NEP` does not cover the even more general case of having a nonlinear dependence on the eigenvector $x$.

In terms of the user interface, `NEP` is quite similar to previously discussed solvers. The main difference is how to represent the function $T$. We will show different alternatives for this.

The `NEP` module of SLEPc has been explained with more detail in {cite:p}`Campos:2021:NEP`, including an algorithmic description of the implemented solvers.

{#sec:nep label="sec:nep"}
## General Nonlinear Eigenproblems

As in previous chapters, we first set up the notation and briefly review basic properties of the eigenvalue problems to be addressed. In this case, we focus on general nonlinear eigenproblems, that is, those that cannot be expressed in a simpler form such as a polynomial eigenproblem. These problems arise in many applications, such as the discretization of delay differential equations. Some of the methods designed to solve such problems are based on Newton-type iterations, so in some ways `NEP` has similarities to PETSc's nonlinear solvers `SNES`. For background material on the nonlinear eigenproblem, the reader is referred to {cite:p}`Guttel:2017:NEP`, {cite:p}`Mehrmann:2004:NEP`.

We consider nonlinear eigenvalue problems of the form

```{math}
:label: eq:nep

T(\lambda)x=0,\qquad x\neq 0,
```

 where $T:\Omega\rightarrow\mathbb{C}^{n\times n}$ is a matrix-valued function that is analytic on an open set of the complex plane $\Omega\subseteq\mathbb{C}$. Assuming that the problem is regular, that is, $\det T(\lambda)$ does not vanish identically, any pair $(\lambda,x)$ satisfying equation {math:numref}`eq:nep` is an eigenpair, where $\lambda\in\Omega$ is the eigenvalue and $x\in\mathbb{C}^n$ is the eigenvector. Linear and polynomial eigenproblems are particular cases of equation {math:numref}`eq:nep`.

An example application is the rational eigenvalue problem

```{math}
:label: eq:rep

-Kx+\lambda Mx+\sum_{j=1}^k\frac{\lambda}{\sigma_j-\lambda}C_jx=0,
```

 arising in the study of free vibration of plates with elastically attached masses. Here, all matrices are symmetric, $K$ and $M$ are positive-definite and $C_j$ have small rank. Another example comes from the discretization of parabolic partial differential equations with time delay $\tau$, resulting in

```{math}
:label: eq:delay

(-\lambda I + A + e^{-\tau\lambda}B)x = 0.
```

### Split Form

Equation {math:numref}`eq:nep` can always be rewritten as

```{math}
:label: eq:split

\big(A_0f_0(\lambda)+A_1f_1(\lambda)+\cdots+A_{\ell-1}f_{\ell-1}(\lambda)\big)x=
\left(\sum_{i=0}^{\ell-1}A_if_i(\lambda)\right)x = 0,
```

 where $A_i$ are $n\times n$ matrices and $f_i:\Omega\rightarrow\mathbb{C}$ are analytic functions. We will call equation {math:numref}`eq:split` the split form of the nonlinear eigenvalue problem. Often, the formulation arising from applications already has this form, as illustrated by the examples above. Also, a polynomial eigenvalue problem fits this form, where in this case the $f_i$ functions are the polynomial bases of degree $i$, either monomial or non-monomial.

## Defining the Problem

The user interface of the `NEP` package is quite similar to `EPS` and `PEP`. As mentioned above, the main difference is the way in which the eigenproblem is defined. In equation [](#sec:nepjac), we focus on the case where the problem is defined as in PETSc's nonlinear solvers `SNES`, that is, providing user-defined callback functions to compute the nonlinear function matrix, $T(\lambda)$, and its derivative, $T'(\lambda)$. We defer the discussion of using the split form of the nonlinear eigenproblem to section [](#sec:nepsplit).

{#sec:nepjac label="sec:nepjac"}
### Using Callback Functions

A sample code for solving a nonlinear eigenproblem with `NEP` is shown in figure [](#fig:ex-nep). The usual steps are performed, starting with the creation of the solver context with `NEPCreate`. Then the problem matrices are defined, see discussion below. The call to `NEPSetFromOptions` captures relevant options specified in the command line. The actual solver is invoked with `NEPSolve`. Then, the solution is retrieved with `NEPGetConverged` and `NEPGetEigenpair`. Finally, `NEPDestroy` destroys the object.

```{code-block} c
:name: fig:ex-nep
:caption: Example code for basic solution with `NEP` using callbacks.

NEP         nep;       /*  eigensolver context */
Mat         F, J;      /*  Function and Jacobian matrices  */
Vec         xr, xi;    /*  eigenvector, x       */
PetscScalar kr, ki;    /*  eigenvalue, k        */
ApplCtx     ctx;       /*  user-defined context */
PetscInt    j, nconv;
PetscReal   error;

NEPCreate( PETSC_COMM_WORLD, &nep );
/* create and preallocate F and J matrices */
NEPSetFunction( nep, F, F, FormFunction, &ctx );
NEPSetJacobian( nep, J, FormJacobian, &ctx );
NEPSetFromOptions( nep );
NEPSolve( nep );
NEPGetConverged( nep, &nconv );
for (j=0; j<nconv; j++) {
  NEPGetEigenpair( nep, j, &kr, &ki, xr, xi );
  NEPComputeError( nep, j, NEP_ERROR_RELATIVE, &error );
}
NEPDestroy( &nep );
```

In `SNES`, the usual way to define a set of nonlinear equations $F(x)=0$ is to provide two user-defined callback functions, one to compute the residual vector, $r=F(x)$ for a given $x$, and another one to evaluate the Jacobian matrix, $J(x)=F'(x)$. In the case of `NEP` there are some differences, since the function $T$ depends on the parameter $\lambda$ only. For a given value of $\lambda$ and its associated vector $x$, the residual vector is defined as

```{math}
:label: eq:nlres

r=T(\lambda)x.
```

 We require the user to provide a callback function to evaluate $T(\lambda)$, rather than computing the residual $r$. Once $T(\lambda)$ has been built, `NEP` solvers can compute its action on any vector $x$. Regarding the derivative, in `NEP` we use $T'(\lambda)$, which will be referred to as the Jacobian matrix by analogy to `SNES`. This matrix must be computed with another callback function.

Hence, both callback functions must compute a matrix. The nonzero pattern of these matrices does not usually change, so they must be created and preallocated at the beginning of the solution process. Then, these `Mat` objects are passed to the solver, together with the pointers to the callback functions, with `NEPSetFunction``NEPSetJacobian`

```{code} c
NEPSetFunction(NEP nep,Mat F,Mat P,PetscErrorCode (*fun)(NEP,PetscScalar,
               Mat,Mat,void*),void *ctx);
NEPSetJacobian(NEP nep,Mat J,PetscErrorCode (*jac)(NEP,PetscScalar,
               Mat,void*),void *ctx)
```

The argument `ctx` is an optional user-defined context intended to contain application-specific parameters required to build $T(\lambda)$ or $T'(\lambda)$, and it is received as the last argument in the callback functions. The callback routines also get an argument containing the value of $\lambda$ at which $T$ or $T'$ must be evaluated. Note that the `NEPSetFunction` callback takes two `Mat` arguments instead of one. The rationale for this is that some `NEP` solvers require to perform linear solves with $T(\lambda)$ within the iteration (in `SNES` this is done with the Jacobian), so $T(\lambda)$ will be passed as the coefficient matrix to a `KSP` object. The second `Mat` argument `P` is the matrix from which the preconditioner is constructed (which is usually the same as `F`).

There is the possibility of solving the problem in a matrix-free fashion, that is, just implementing subroutines that compute the action of $T(\lambda)$ or $T'(\lambda)$ on a vector, instead of having to explicitly compute all nonzero entries of these two matrices. The SLEPc distribution contains an example illustrating this, using the concept of *shell* matrices (see section [](#sec:supported) for details).

#### Parameters for Problem Definition.

Once $T$ and $T'$ have been set up, the definition of the problem is completed with the number and location of the eigenvalues to compute, in a similar way as eigensolvers discussed in previous chapters.

The number of requested eigenvalues (and eigenvectors), `nev`, is established with `NEPSetDimensions`

```{code} c
NEPSetDimensions(NEP nep,PetscInt nev,PetscInt ncv,PetscInt mpd);
```

By default, `nev`=1 (and some solvers will return only one eigenpair, even if a larger `nev` is requested). The other two arguments control the dimension of the subspaces used internally (the number of column vectors, `ncv`, and the maximum projected dimension, `mpd`), although they are relevant only in eigensolvers based on subspace projection (basic algorithms ignore them). There are command-line keys for these parameters: `-nep_nev`, `-nep_ncv` and `-nep_mpd`.

:::{table} Available possibilities for selection of the eigenvalues of interest in `NEP`.
:name: tab:portionn

 |`NEPWhich`                |Command line key           |Sorting criterion
 |--------------------------|---------------------------|----------------------------------------
 |`NEP_LARGEST_MAGNITUDE`   |`-nep_largest_magnitude`   |Largest $\|\lambda\|$
 |`NEP_SMALLEST_MAGNITUDE`  |`-nep_smallest_magnitude`  |Smallest $\|\lambda\|$
 |`NEP_LARGEST_REAL`        |`-nep_largest_real`        |Largest $\mathrm{Re}(\lambda)$
 |`NEP_SMALLEST_REAL`       |`-nep_smallest_real`       |Smallest $\mathrm{Re}(\lambda)$
 |`NEP_LARGEST_IMAGINARY`   |`-nep_largest_imaginary`   |Largest $\mathrm{Im}(\lambda)$
 |`NEP_SMALLEST_IMAGINARY`  |`-nep_smallest_imaginary`  |Smallest $\mathrm{Im}(\lambda)$
 |`NEP_TARGET_MAGNITUDE`    |`-nep_target_magnitude`    |Smallest $\|\lambda-\tau\|$
 |`NEP_TARGET_REAL`         |`-nep_target_real`         |Smallest $\|\mathrm{Re}(\lambda-\tau)\|$
 |`NEP_TARGET_IMAGINARY`    |`-nep_target_imaginary`    |Smallest $\|\mathrm{Im}(\lambda-\tau)\|$
 |`NEP_ALL`                 |`-nep_all`                 |All $\lambda\in\Omega$

:::

For the selection of the portion of the spectrum of interest, there are several alternatives listed in table [](#tab:portionn), to be selected with the function `NEPSetWhichEigenpairs`

```{code} c
NEPSetWhichEigenpairs(NEP nep,NEPWhich which);
```

The default is to compute the largest magnitude eigenvalues. For the sorting criteria relative to a target value, $\tau$ must be specified with `NEPSetTarget` or in the command-line with `-nep_target`.

`NEP` solvers can also work with a region of the complex plane (`RG`), as discussed in section [](#sec:region) for linear problems. Some eigensolvers (NLEIGS) use the definition of the region to compute `nev` eigenvalues in its interior. If *all* eigenvalues inside the region are required, then a contour-integral method is required, see discussion in {cite:t}`str-11`.

#### Left Eigenvectors

As in the case of linear eigensolvers, some `NEP` solvers have two-sided variants to compute also left eigenvectors. In the case of `NEP`, left eigenvectors are defined as

```{math}
:label: eq:nepleft

y^*T(\lambda)=0^*,\qquad y\neq 0.
```

 Two-sided variants can be selected with `NEPSetTwoSided`

```{code} c
NEPSetTwoSided(NEP eps,PetscBool twosided);
```

{#sec:nepsplit label="sec:nepsplit"}
### Expressing the NEP in Split Form

```{code-block} c
:name: fig:ex-split
:caption: Example code for defining the `NEP` eigenproblem in the split form.

FNCreate(PETSC_COMM_WORLD,&f1);  /* f1 = -lambda */
FNSetType(f1,FNRATIONAL);
coeffs[0] = -1.0; coeffs[1] = 0.0;
FNRationalSetNumerator(f1,2,coeffs);

FNCreate(PETSC_COMM_WORLD,&f2);  /* f2 = 1 */
FNSetType(f2,FNRATIONAL);
coeffs[0] = 1.0;
FNRationalSetNumerator(f2,1,coeffs);

FNCreate(PETSC_COMM_WORLD,&f3);  /* f3 = exp(-tau*lambda) */
FNSetType(f3,FNEXP);
FNSetScale(f3,-tau,1.0);

mats[0] = A;  funs[0] = f2;
mats[1] = Id; funs[1] = f1;
mats[2] = B;  funs[2] = f3;
NEPSetSplitOperator(nep,3,mats,funs,SUBSET_NONZERO_PATTERN);
```

Instead of implementing callback functions for $T(\lambda)$ and $T'(\lambda)$, a usually simpler alternative is to use the split form of the nonlinear eigenproblem, equation {math:numref}`eq:split`. Note that in split form, we have $T'(\lambda)=\sum_{i=0}^{\ell-1}A_if'_i(\lambda)$, so the derivatives of $f_i(\lambda)$ are also required. As described below, we will represent each of the analytic functions $f_i$ by means of an auxiliary object `FN` that holds both the function and its derivative.

Hence, for the split form representation we must provide $\ell$ matrices $A_i$ and the corresponding functions $f_i(\lambda)$, by means of `NEPSetSplitOperator`

```{code} c
NEPSetSplitOperator(NEP nep,PetscInt l,Mat A[],FN f[],MatStructure str);
```

Here, the `MatStructure` flag is used to indicate whether all matrices have the same (or subset) nonzero pattern with respect to the first one. Figure [](#fig:ex-split) illustrates this usage with the problem of equation {math:numref}`eq:delay`, where $\ell=3$ and the matrices are $I$, $A$ and $B$ (note that in the code we have changed the order for efficiency reasons, since the nonzero pattern of $I$ and $B$ is a subset of $A$'s in this case). Two of the associated functions are polynomials ($-\lambda$ and $1$) and the other one is the exponential $e^{-\tau\lambda}$.

Note that using the split form is required in order to be able to use some eigensolvers, in particular, those that project the nonlinear eigenproblem onto a low dimensional subspace and then use a dense nonlinear solver for the projected problem.

Details of how to define the $f_i$ functions by using the `FN` class are provided in section [](#sec:sys).

:::{table} Problem types considered in `NEP`.
:name: tab:ntypeq

 |Problem Type  |`NEPProblemType`  |Command line key
 |--------------|------------------|------------------
 |General       |`NEP_GENERAL`     |`-nep_general`
 |Rational      |`NEP_RATIONAL`    |`-nep_rational`

:::

When defining the problem in split form, it may also be useful to specify a problem type. For example, if the user knows that all $f_i$ functions are rational, as in equation {math:numref}`eq:rep`, then setting the problem type to `NEP_RATIONAL` gives a hint to the solver that may simplify the solution process. The problem types currently supported for `NEP` are listed in table [](#tab:ntypeq). When in doubt, use the default problem type (`NEP_GENERAL`).

The problem type can be specified at run time with the corresponding command line key or, more usually, within the program with the function `NEPSetProblemType`

```{code} c
NEPSetProblemType(NEP nep,NEPProblemType type);
```

Currently, the problem type is ignored in most solvers and it is taken into account only in NLEIGS for determining singularities automatically.

## Selecting the Solver

The solution method can be specified procedurally with `NEPSetType`

```{code} c
NEPSetType(NEP nep,NEPType method);
```

or via the options database command `-nep_type` followed by the name of the method (see table [](#tab:solversn)). The methods currently available in `NEP` are the following:

-   Residual inverse iteration (RII), where in each iteration the eigenvector correction is computed as $T(\sigma)^{-1}$ times the residual $r$.

-   Successive linear problems (SLP), where in each iteration a linear eigenvalue problem $T(\tilde\lambda)\tilde x=\mu T'(\tilde\lambda)\tilde x$ is solved for the eigenvalue correction $\mu$.

-   Nonlinear Arnoldi, which builds an orthogonal basis $V_j$ of a subspace expanded with the vectors generated by RII, then chooses the approximate eigenpair $(\tilde\lambda,\tilde x)$ such that $\tilde x=V_jy$ and $V_j^*T(\tilde\lambda)V_jy=0$.

-   NLEIGS, which is based on a (rational) Krylov iteration operating on a companion-type linearization of a rational interpolant of the nonlinear function.

-   CISS, a contour-integral solver that allows computing all eigenvalues in a given region.

-   Polynomial interpolation, where a matrix polynomial $P(\lambda)$ is built by evaluating $T(\cdot)$ at a few points, then `PEP` is used to solve the polynomial eigenproblem.

:::{table} Nonlinear eigenvalue solvers available in the `NEP` module.
:name: tab:solversn

 |Method                      |`NEPType`      |Options Database |  Need $T'(\cdot)$ | Two-sided
 |----------------------------|---------------|-----------------|-------------------|--------------
 |Residual inverse iteration  |`NEPRII`       |`rii`            |       no          |
 |Successive linear problems  |`NEPSLP`       |`slp`            |       yes         |    yes
 |Nonlinear Arnoldi           |`NEPNARNOLDI`  |`narnoldi`       |       no          |
 |Rational Krylov (NLEIGS)    |`NEPNLEIGS`    |`nleigs`         |       no          |    yes
 |Contour integral SS         |`NEPCISS`      |`ciss`           |       yes         |
 |Polynomial interpolation    |`NEPINTERPOL`  |`interpol`       |       no          |

:::

The `NEPSLP` method performs a linearization that results in a (linear) generalized eigenvalue problem. This is handled by an `EPS` object created internally. If required, this `EPS` object can be extracted with the operation `NEPSLPGetEPS`

```{code} c
NEPSLPGetEPS(NEP nep,EPS *eps);
```

This allows the application programmer to set any of the `EPS` options directly within the code. These options can also be set through the command-line, simply by prefixing the `EPS` options with `-nep_slp_`.

Similarly, `NEPINTERPOL` works with a `PEP` object internally, that can be retrieved by `NEPInterpolGetPEP`. Another relevant option of this solver is the degree of the interpolation polynomial, that can be set with `NEPInterpolSetInterpolation`

```{code} c
NEPInterpolSetInterpolation(NEP nep,PetscReal tol,PetscInt deg);
```

The polynomial interpolation solver currently uses Chebyshev polynomials of the 1st kind and requires the user to specify an interval of the real line where the eigenvalues must be computed, e.g.

```{code} console
$ ./ex22 -nep_type interpol -rg_interval_endpoints 0.1,14.0,-0.1,0.1
         -nep_nev 2 -nep_interpol_interpolation_degree 15 -nep_target 1.0
```

For details about specifying a region, see section [](#sec:sys).

Some solvers such as `NEPRII` and `NEPNARNOLDI` need a `KSP` object to handle the solution of linear systems of equations. This `KSP` and can be retrieved with e.g. `NEPRIIGetKSP`

```{code} c
NEPRIIGetKSP(NEP nep,KSP *ksp);
```

This `KSP` object is typically used to compute the action of $T(\sigma)^{-1}$ on a given vector. In principle, $\sigma$ is an approximation of an eigenvalue, but it is usually more efficient to keep this value constant, otherwise the factorization or preconditioner must be recomputed every time since eigensolvers update eigenvalue approximations in each iteration. This behaviour can be changed with `NEPSetLagPreconditioner`

```{code} c
NEPRIISetLagPreconditioner(NEP nep,PetscInt lag);
```

Recomputing the preconditioner every 2 iterations, say, will introduce a considerable overhead, but may reduce the number of iterations significantly. Another related comment is that, when using an iterative linear solver, the requested accuracy is adapted as the outer iteration progresses, being the tolerance larger in the first solves. Again, the user can modify this behaviour with `NEPRIISetConstCorrectionTol`. Both options can also be changed at run time. As an example, consider the following command line:

```{code} console
$ ./ex22 -nep_type rii -nep_rii_lag_preconditioner 2
         -nep_rii_ksp_type bcgs -nep_rii_pc_type ilu
         -nep_rii_const_correction_tol 1 -nep_rii_ksp_rtol 1e-3
```

The example uses RII with BiCGStab plus ILU, where the preconditioner is updated every two outer iterations and linear systems are solved up to a tolerance of $10^{-3}$.

The NLEIGS solver is most appropriate for problems where $T(\cdot)$ is singular at some known parts of the complex plane, for instance the case that $T(\cdot)$ contains $\sqrt{\lambda}$. To treat this case effectively, the NLEIGS solver requires a discretization of the singularity set, which can be provided by the user in the form of a callback function: `NEPNLEIGSSetSingularitiesFunction`

```{code} c
NEPNLEIGSSetSingularitiesFunction(NEP nep,PetscErrorCode (*fun)
                          (NEP,PetscInt*,PetscScalar*,void*),void *ctx);
```

Alternatively, if the problem is known to be a rational eigenvalue problem, the user can avoid the computation of singularities by just specifying the problem type with `NEPSetProblemType`, as explained at the end of the previous section. If none of the above functions is invoked by the user, then the NLEIGS solver attempts to determine the singularities automatically.

## Retrieving the Solution

The procedure for obtaining the computed eigenpairs is similar to previously discussed eigensolvers. After the call to `NEPSolve`, the computed results are stored internally and a call to `NEPGetConverged` must be issued to obtain the number of converged solutions. Then calling `NEPGetEigenpair` repeatedly will retrieve each eigenvalue-eigenvector pair.

`NEPGetEigenpair`

```{code} c
NEPGetEigenpair(NEP nep,PetscInt j,PetscScalar *kr,PetscScalar *ki,
                Vec xr,Vec xi);
```

In two-sided solvers (see table [](#tab:solversn)), it is also possible to retrieve left eigenvectors with `NEPGetLeftEigenvector`

```{code} c
NEPGetLeftEigenvector(NEP nep,PetscInt j,Vec yr,Vec yi);
```

**Note about real/complex scalar versions**: The interface makes provision for returning a complex eigenvalue (or eigenvector) when doing the computation in a PETSc/SLEPc version built with real scalars, as is done in other eigensolvers such as `EPS`. However, in some cases this will not be possible. In particular, when callback functions are used and a complex eigenvalue approximation is hit, the solver will fail unless configured with complex scalars. The reason is that the user interface for callback functions only have a single `PetscScalar lambda` argument and hence cannot handle complex arguments in real arithmetic.

The function `NEPComputeError`

```{code} c
NEPComputeError(NEP nep,PetscInt j,NEPErrorType type,PetscReal *error);
```

can be used to assess the accuracy of the computed solutions. The error is based on the 2-norm of the residual vector $r$ defined in equation {math:numref}`eq:nlres`.

As in the case of `EPS`, in `NEP` the number of iterations carried out by the solver can be determined with `NEPGetIterationNumber`, and the tolerance and maximum number of iterations can be set with `NEPSetTolerances`. Also, convergence can be monitored with either textual monitors `-nep_monitor`, `-nep_monitor_all`, `-nep_monitor_conv`, or graphical monitors `-nep_monitor draw::draw_lg`, `-nep_monitor_all draw::draw_lg`. See section [](#sec:monitor) for additional details. Similarly, there is support for viewing the computed solution as explained in section [](#sec:epsviewers).

The `NEP` class also provides some kind of iterative refinement, similar to the one available in `PEP`, see section [](#sec:refine). The parameters can be set with `NEPSetRefine`

```{code} c
NEPSetRefine(NEP nep,NEPRefine refine,PetscInt npart,
             PetscReal tol,PetscInt its,NEPRefineScheme scheme);
```

```{rubric} Footnotes
```

```{eval-rst}
.. bibliography::
   :filter: docname in docnames
```
