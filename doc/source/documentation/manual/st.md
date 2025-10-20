(ch:st)=
# ST: Spectral Transformation

The Spectral Transformation (`ST`) is the object that encapsulates the functionality required for acceleration techniques based on the transformation of the spectrum. Most eigensolvers in `EPS` work by applying an operator to a set of vectors and this operator can adopt different forms. The `ST` object handles all the different possibilities in a uniform way, so that the solver can proceed without knowing which transformation has been selected. The spectral transformation can be specified at run time, together with related options such as which linear solver to use.

Despite being a rather unrelated concept, `ST` is also used to handle the preconditioners and correction-equation solvers used in preconditioned eigensolvers such as GD and JD.

The description in this chapter focuses on the use of `ST` in the context of `EPS`. For usage within other solver classes, we will provide further details, e.g., shift-and-invert for polynomial eigenproblems in section [](#sec:qst).

## General Description

Spectral transformations are powerful tools for adjusting the way in which eigensolvers behave when coping with a problem. The general strategy consists in transforming the original problem into a new one in which eigenvalues are mapped to a new position while eigenvectors remain unchanged. These transformations can be used with several goals in mind:

-   Compute internal eigenvalues. In some applications, the eigenpairs of interest are not the extreme ones (largest magnitude, rightmost, leftmost), but those contained in a certain interval or those closest to a certain value of the complex plane.

-   Accelerate convergence. Convergence properties typically depend on how close the eigenvalues are from each other. With some spectral transformations, difficult eigenvalue distributions can be remapped in a more favorable way in terms of convergence.

-   Handle some special situations. For instance, in generalized problems when matrix $B$ is singular, it may be necessary to use a spectral transformation.

SLEPc separates spectral transformations from solution methods so that any combination of them can be specified by the user. To achieve this, most eigensolvers contained in `EPS` are implemented in such a way that they are independent of which transformation has been selected by the user (the exception are preconditioned solvers, see below). That is, the solver algorithm has to work with a generic operator, whose actual form depends on the transformation used. After convergence, eigenvalues are transformed back appropriately.

For technical details of the transformations described in this chapter, the interested user is referred to {cite:p}`Eri80,Sco82,Nou87,Mee94`.

**Preconditioners**:
As explained in the previous chapter, `EPS` contains preconditioned eigensolvers such as GD or JD. These solvers either apply a preconditioner at a certain step of the computation, or need to solve a correction equation with a preconditioned linear solver. One of the main goals of these solvers is to achieve a similar effect as an inverse-based spectral transformation such as shift-and-invert, but with less computational cost. For this reason, a "preconditioner" spectral transformation has been included in the `ST` object. However, this is just a convenient way of organizing the functionality, since this fake spectral transform cannot be used with non-preconditioned eigensolvers, and conversely preconditioned eigensolvers cannot be used with conventional spectral transformations.

## Basic Usage

The `ST` module is the analog of some PETSc modules such as {external:doc}`PC`. The user does not usually need to create a stand-alone `ST` object explicitly. Instead, every `EPS` object internally sets up an associated `ST`. Therefore, the usual object management methods such as `STCreate()`, `STDestroy()`, `STView()`, `STSetFromOptions()`, are not usually called by the user.

Although the `ST` context is hidden inside the `EPS` object, the user still has control over all the options, by means of the command line, or also inside the program. To allow application programmers to set any of the spectral transformation options directly within the code, the following function is provided to extract the `ST` context from the `EPS` object:

```{code} c
EPSGetST(EPS eps,ST *st);
```

After this, one is able to set any options associated with the `ST` object. For example, to set the value of the shift, $\sigma$, the following function is available:

```{code} c
STSetShift(ST st,PetscScalar shift);
```

This can also be done with the command line option `-st_shift <shift>`. Note that the argument `shift` is defined as a {external:doc}`PetscScalar`, and this means that complex shifts are not allowed unless the complex version of SLEPc is used.

:::{warning}
Usually, `STSetShift()` is never called from application code, as the value of the shift is taken from the target. In all transformations except `STSHIFT`, there is a direct connection between the target $\tau$ (described in section [](#sec:which)) and the shift $\sigma$, as will be discussed below. The normal usage is that the user sets the target and then $\sigma$ is set to $\tau$ automatically (though it is still possible for the user to set a different value of the shift).
:::

Other object operations are available, which are not usually called by the user. The most important ones are `STApply()`, that applies the operator to a vector, and `STSetUp()`, that prepares all the necessary data structures before the solution process starts. The term "operator" refers to one of $A$, $B^{-1}\!A$, $A-\sigma I$, \... depending on which kind of spectral transformation is being used.

## Available Transformations

This section describes the spectral transformations that are provided in SLEPc. As in the case of eigensolvers, the spectral transformation to be used can be specified procedurally or via the command line. The application programmer can set it by means of:

```{code} c
STSetType(ST st,STType type);
```

The `ST` type can also be set with the command-line option `-st_type` followed by the name of the method (see table [](#tab:transforms)). The first five spectral transformations are described in detail in the rest of this section. The last possibility, `STSHELL`, uses a specific, application-provided spectral transformation. Section [](#sec:shell) in the final chapter describes how to implement one of these transformations.

:::{table} Spectral transformations available in the `ST` package
:name: tab:transforms

  Spectral Transformation  |`STType`     |Options Name |            Operator
  :------------------------|:------------|:------------|:-------------------------------:
  Shift of Origin          |`STSHIFT`    |`shift`      |          $B^{-1}A-\sigma I$
  Shift-and-invert         |`STSINVERT`  |`sinvert`    |         $(A-\sigma B)^{-1}B$
  Generalized Cayley       |`STCAYLEY`   |`cayley`     |     $(A-\sigma B)^{-1}(A+\nu B)$
  Preconditioner           |`STPRECOND`  |`precond`    |   $K^{-1}\approx(A-\sigma B)^{-1}$
  Polynomial Filter        |`STFILTER`   |`filter`     |                $p(A)$
  Shell Transformation     |`STSHELL`    |`shell`      |            *user-defined*

:::

The last column of table [](#tab:transforms) shows a general form of the operator used in each case. This generic operator can adopt different particular forms depending on whether the eigenproblem is standard or generalized, or whether the value of the shift ($\sigma$) and anti-shift ($\nu$) is zero or not. All the possible combinations are listed in table [](#tab:op).

:::{table} Operators used in each spectral transformation mode
:name: tab:op

| `ST`       |Choice of $\sigma,\nu$    |         Standard problem         |       Generalized problem        |
|:-----------|:-------------------------|:--------------------------------:|:--------------------------------:|
| `shift`    |$\sigma=0$                |               $A$                |            $B^{-1}A$             |
| ` `        |$\sigma\not=0$            |           $A-\sigma I$           |        $B^{-1}A-\sigma I$        |
| `sinvert`  |$\sigma=0$                |             $A^{-1}$             |            $A^{-1}B$             |
| ` `        |$\sigma\not=0$            |       $(A-\sigma I)^{-1}$        |       $(A-\sigma B)^{-1}B$       |
| `cayley`   |$\sigma\not=0,\nu=0$      |       $(A-\sigma I)^{-1}A$       |       $(A-\sigma B)^{-1}A$       |
| ` `        |$\sigma=0,\nu\not=0$      |          $I+\nu A^{-1}$          |         $I+\nu A^{-1}B$          |
| ` `        |$\sigma\not=0,\nu\not=0$  |   $(A-\sigma I)^{-1}(A+\nu I)$   |   $(A-\sigma B)^{-1}(A+\nu B)$   |
| `precond`  |$\sigma=0$                |      $K^{-1}\approx A^{-1}$      |      $K^{-1}\approx A^{-1}$      |
| ` `        |$\sigma\not=0$            | $K^{-1}\approx(A-\sigma I)^{-1}$ | $K^{-1}\approx(A-\sigma B)^{-1}$ |

:::

The expressions shown in table [](#tab:op) are not built explicitly. Instead, the appropriate operations are carried out when applying the operator to a certain vector. The inverses imply the solution of a system of linear equations that is managed by setting up an associated {external:doc}`KSP` object. The user can control the behavior of this object by adjusting the appropriate options, as will be illustrated with examples in section [](#sec:lin).

### Shift of Origin

By default, no spectral transformation is performed. This is equivalent to a shift of origin (`STSHIFT`) with $\sigma=0$, that is, the first line of table [](#tab:op). The solver works with the original expressions of the eigenvalue problems,

```{math}
:label: eq:std-eig-problem

Ax=\lambda x,
```

 for standard problems, and $Ax=\lambda Bx$ for generalized ones. Note that this last equation is actually treated internally as

```{math}
:label: eq:gen-eig-problem

B^{-1}Ax=\lambda x.
```

 When the eigensolver in `EPS` requests the application of the operator to a vector, a matrix-vector multiplication by matrix $A$ is carried out (in the standard case) or a matrix-vector multiplication by matrix $A$ followed by a linear system solve with coefficient matrix $B$ (in the generalized case). Note that in the last case, the operation will fail if matrix $B$ is singular.

When the shift, $\sigma$, is given a value different from the default, 0, the effect is to move the whole spectrum by that exact quantity, $\sigma$, which is called *shift of origin*. To achieve this, the solver works with the shifted matrix, that is, the expressions it has to cope with are

```{math}
:label: eq:std-eig-problem-shift

(A-\sigma I)x=\theta x,
```

 for standard problems, and

```{math}
:label: eq:gen-eig-problem-shift

(B^{-1}A-\sigma I) x=\theta x,
```

 for generalized ones. The important property that is used is that shifting does not alter the eigenvectors and that it does change the eigenvalues in a simple known way, it shifts them by $\sigma$. In both the standard and the generalized problems, the following relation holds

```{math}
:label: eq:eig-problem-shift

\theta=\lambda-\sigma.
```

This means that after the solution process, the value $\sigma$ has to be added[^slepc-v3.5] to the computed eigenvalues, $\theta$, in order to retrieve the solution of the original problem, $\lambda$. This is done by means of the function `STBackTransform()`, which does not need to be called directly by the user.

### Shift-and-invert

```{figure} ../../_static/images/manual/svg/fig-sinvert.svg
:alt: The shift-and-invert spectral transformation
:name: fig:sinvert

The shift-and-invert spectral transformation
```

The shift-and-invert spectral transformation (`STSINVERT`) is used to enhance convergence of eigenvalues in the neighborhood of a given value. In this case, the solver deals with the expressions

```{math}
:label: eq:std-eig-problem-shift-and-invert

(A-\sigma I)^{-1}x=\theta x,
```

```{math}
:label: eq:gen-eig-problem-shift-and-invert

(A-\sigma B)^{-1}B x=\theta x,
```

 for standard and generalized problems, respectively. This transformation is effective for finding eigenvalues near $\sigma$ since the eigenvalues $\theta$ of the operator that are largest in magnitude correspond to the eigenvalues $\lambda$ of the original problem that are closest to the shift $\sigma$ in absolute value, as illustrated in figure [](#fig:sinvert) for an example with real eigenvalues. Once the wanted eigenvalues have been found, they may be transformed back to eigenvalues of the original problem. Again, the eigenvectors remain unchanged. In this case, the relation between the eigenvalues of both problems is

```{math}
:label: eq:eig-problem-shift-and-invert

\theta=1/(\lambda-\sigma).
```

Therefore, after the solution process, the operation to be performed in function `STBackTransform()` is $\lambda=\sigma+1/\theta$ for each of the computed eigenvalues.

This spectral transformation is used in the spectrum slicing technique, see section [](#sec:slice).

{#sec:cayley}
### Cayley

The generalized Cayley transform (`STCAYLEY`) is defined from the expressions

```{math}
:label: eq:std-eig-problem-cayley

(A-\sigma I)^{-1}(A+\nu I)x=\theta x,
```

```{math}
:label: eq:gen-eig-problem-cayley

(A-\sigma B)^{-1}(A+\nu B)x=\theta x,
```

 for standard and generalized problems, respectively. Sometimes, the term Cayley transform is applied for the particular case in which $\nu=\sigma$. This is the default if $\nu$ is not given a value explicitly. The value of $\nu$ (the anti-shift) can be set with the following function:

```{code} c
STCayleySetAntishift(ST st,PetscScalar nu);
```

or in the command line with `-st_cayley_antishift`.

This transformation is mathematically equivalent to shift-and-invert and, therefore, it is effective for finding eigenvalues near $\sigma$ as well. However, in some situations it is numerically advantageous with respect to shift-and-invert (see {cite:p}`Bai00{11.2},Leh01`).

In this case, the relation between the eigenvalues of both problems is

```{math}
:label: eq:eig-problem-cayley

\theta=(\lambda+\nu)/(\lambda-\sigma).
```

 Therefore, after the solution process, the operation to be performed in function `STBackTransform()` is $\lambda=(\theta\sigma+\nu)/(\theta-1)$ for each of the computed eigenvalues.

{#sec:precond}
### Preconditioner

As mentioned in the introduction of this chapter, the special type `STPRECOND` is used for handling preconditioners or preconditioned iterative linear solvers, which are used in the context of preconditioned eigensolvers for expanding the subspace. For instance, in the GD solver the so-called correction vector $d_i$ to be added to the subspace in each iteration is computed as

```{math}
:label: eq:correction-vector

d_i=K^{-1}P_i(A-\theta_i B)x_i,
```

 where $(\theta_i,x_i)$ is the current approximation of the sought-after eigenpair, and $P_i$ is a projector involving $x_i$ and $K^{-1}x_i$. In the above expressions, $K$ is a preconditioner matrix that is built from $A-\theta_i B$. However, since $\theta_i$ changes at each iteration, which would force recomputation of the preconditioner, we opt for using

```{math}
:label: eq:precon

K^{-1}\approx (A-\sigma B)^{-1}.
```

Similarly, in the JD eigensolver the expansion of the subspace is carried out by solving a correction equation similar to

```{math}
:label: eq:jd-correction

(I-x_ix_i^*)(A-\theta_i B)(I-x_ix_i^*)d_i=-(A-\theta_i B)x_i,
```

 where the system is solved approximately with a preconditioned iterative linear solver. For building the preconditioner of this linear system, the projectors $I-x_ix_i^*$ are ignored, and again it is not recomputed in each iteration. Therefore, the preconditioner is built as in equation {math:numref}`eq:precon` as well.

It should be clear from the previous discussion, that `STPRECOND` does not work in the same way as the rest of spectral transformations. In particular, it does not rely on `STBackTransform()`. It is rather a convenient mechanism for handling the preconditioner and linear solver (see examples in section [](#sec:lin)). The expressions shown in tables [](#tab:transforms) and [](#tab:op) are just a reference to indicate from which matrix the preconditioner is built by default.

There is the possibility that the user overrides the default behavior, that is, to explicitly supply a matrix from which the preconditioner is to be built, with:

```{code} c
STSetPreconditionerMat(ST st,Mat mat);
```

The above function can also be used in other spectral transformations such as shift-and-invert in case the user has a cheap approximation $K$ of the coefficient matrix $A-\sigma B$. An alternative is to pass approximations of both $A$ and $B$ so that `ST` builds the preconditioner matrix internally, with:

```{code} c
STSetSplitPreconditioner(ST st,PetscInt n,Mat Psplit[],MatStructure strp);
```

Note that preconditioned eigensolvers in `EPS` select `STPRECOND` by default, so the user does not need to specify it explicitly.

{#sec:filter}
### Polynomial Filtering

The type `STFILTER` is also special. It is used in the case of standard symmetric (or Hermitian) eigenvalue problems when the eigenvalues of interest are interior to the spectrum and we want to avoid the high cost associated with the matrix factorization of the shift-and-invert spectral transformation. The techniques generically known as *polynomial filtering* aim at this goal.

The polynomial filtering methods address the eigenvalue problem

```{math}
:label: eq:polyfilt

p(A)x=\theta x,
```

 where $p(\cdot)$ is a suitable high-degree polynomial. Once the polynomial is built, the eigensolver relies on `STApply()` to compute approximations of the eigenvalues $\theta$ of the transformed problem. These approximations must be processed in some way in order to recover the $\lambda$ eigenvalues. Note that in this case there is no `STBackTransform()` operation. Details of the method can be found in {cite:p}`Fan12`.

Currently, SLEPc provides several types of polynomial filtering techniques, which can be selected via `STFilterSetType()`. Note that the external package EVSL also implements polynomial filters to compute all eigenvalues in an interval.

## Advanced Usage

Using the `ST` object is very straightforward. However, when using spectral transformations many things are happening behind the scenes, mainly the solution of systems of linear equations. The user must be aware of what is going on in each case, so that it is possible to guide the solution process in the most beneficial way. This section describes several advanced aspects that can have a considerable impact on efficiency.

{#sec:lin}
### Solution of Linear Systems

In many of the cases shown in table [](#tab:op), the operator contains an inverted matrix, which means that a system of linear equations must be solved whenever the application of the operator to a vector is required. These cases are handled internally by means of a {external:doc}`KSP` object.

In the simplest case, a generalized problem is to be solved with a zero shift. Suppose you run a program that solves a generalized eigenproblem, with default options:

```{code} console
$ ./program
```

In this case, the `ST` object associated with the `EPS` solver creates a {external:doc}`KSP` object whose coefficient matrix is $B$. By default, this {external:doc}`KSP` object is set to use a direct solver, in particular an LU factorization. However, default settings can be changed, as illustrated below.

The following command-line is equivalent to the previous one:

```{code} console
$ ./program -st_ksp_type preonly -st_pc_type lu
```

The two options specify the type of the linear solver and preconditioner to be used. The `-st_` prefix indicates that the option corresponds to the linear solver within `ST`. The combination `preonly`$+$`lu` instructs to use a direct solver (LU factorization, see {{'[PETSc documentation](https://petsc.org/{}/manual/ksp/)'.format(branch)}} for details), so this is the same as the default. Adding a new option changes the default behavior, for instance

```{code} console
$ ./program -st_ksp_type preonly -st_pc_type lu -st_pc_factor_mat_solver_type mumps
```

In this case, an external linear solver package is used (MUMPS, see {{'[PETSc documentation](https://petsc.org/{}/overview/linear_solve_table/)'.format(branch)}} for other available packages). Note that an external package is required for computing a matrix factorization in parallel, since PETSc itself only provides sequential direct linear solvers.

Instead of a direct linear solver, it is possible to use an iterative solver. This may be necessary in some cases, specially for very large problems. However, the user is warned that using an iterative linear solver makes the overall solution process less robust (see also the discussion of preconditioned eigensolvers below). As an example, the command-line

```{code} console
$ ./program -st_ksp_type gmres -st_pc_type bjacobi -st_ksp_rtol 1e-9
```

selects the GMRES solver with block Jacobi preconditioning. In the case of iterative solvers, it is important to use an appropriate tolerance, usually slightly more stringent for the linear solves relative to the desired accuracy of the eigenvalue calculation ($10^{-9}$ in the example, compared to $10^{-8}$ for the eigensolver).

Although the direct solver approach may seem too costly, note that the factorization is only carried out at the beginning of the eigenvalue calculation and this cost is amortized in each subsequent application of the operator. This is also the case for iterative methods with preconditioners with high-cost set-up such as ILU.

The application programmer is able to set the desired linear systems solver options also from within the code. In order to do this, first the context of the {external:doc}`KSP` object must be retrieved with the following function:

```{code} c
STGetKSP(ST st,KSP *ksp);
```

The above functionality is also applicable to the other spectral transformations. For instance, for the shift-and-invert technique with $\tau=10$ using BiCGStab+Jacobi:

```{code} console
$ ./program -st_type sinvert -eps_target 10 -st_ksp_type bcgs -st_pc_type jacobi
```

In shift-and-invert and Cayley, unless $\sigma=0$, the coefficient matrix is not a simple matrix but an expression that can be explicitly constructed or not, depending on the user's choice. This issue is examined in detail in section [](#sec:explicit) below.

:::{note}
By default the preconditioner is built from the matrix $A-\sigma I$ (or $A-\sigma B$ in generalized problems). In some special situations, the user may want to build the preconditioner from a different matrix, e.g., stemming from simplified versions of $A$ and $B$. To do this, one should use the functions mentioned in section [](#sec:precond).
:::

In many cases, especially if a shift-and-invert or Cayley transformation is being used, iterative methods may not be well suited for solving linear systems (because of the properties of the coefficient matrix that can be indefinite and ill-conditioned). When using an iterative linear solver, it may be helpful to run with the option `-st_ksp_converged_reason`, which will display the number of iterations required in each operator application. In extreme cases, the iterative solver fails, so `EPSSolve()` aborts with an error

```{code} c
[0]PETSC ERROR: KSP did not converge (reason=DIVERGED_ITS)!
```

If this happens, it is necessary to use a direct method for solving the linear systems, as explained above.

#### The Case of Preconditioned Eigensolvers

The {external:doc}`KSP` object contained internally in `ST` is also used for applying the preconditioner or solving the correction equation in preconditioned eigensolvers.

The GD eigensolver employs just a preconditioner. Therefore, by default it sets the {external:doc}`KSP` type to `preonly` (no other {external:doc}`KSP` is allowed) and the {external:doc}`PC` type to `jacobi`. The user may change the preconditioner, for example as

```{code} console
$ ./ex5 -eps_type gd -st_pc_type asm
```

The JD eigensolver uses both an iterative linear solver and a preconditioner, so both `KSP` and `PC` are meaningful in this case. It is important to note that, contrary to the ordinary spectral transformations where a direct linear solver is recommended, in JD using an iterative linear solver is usually better than a direct solver. Indeed, the best performance may be achieved with a few iterations of the linear solver (or a large tolerance). For instance, the next example uses JD with GMRES+Jacobi limiting to 10 the number of allowed iterations for the linear solver:

```{code} console
$ ./ex5 -eps_type jd -st_ksp_type gmres -st_pc_type jacobi -st_ksp_max_it 10
```

A discussion on the different options available for the Davidson solvers can be found in {cite:p}`Rom14`.

{#sec:explicit}
### Explicit Computation of the Coefficient Matrix

Three possibilities can be distinguished regarding the form of the coefficient matrix of the systems of linear equations associated with the different spectral transformations. The possible coefficient matrices are:

-   Simple: $B$.

-   Shifted: $A-\sigma I$.

-   Axpy: $A-\sigma B$.

The first case has already been described and presents no difficulty. In the other two cases, there are three possible approaches:

"`shell`"

:   To work with the corresponding expression without forming the matrix explicitly. This is achieved by internally setting a matrix-free matrix with {external:doc}`MatCreateShell`().

"`inplace`"

:   To build the coefficient matrix explicitly. This is done by means of a {external:doc}`MatShift`() or a {external:doc}`MatAXPY`() operation, which overwrites matrix $A$ with the corresponding expression. This alteration of matrix $A$ is reversed after the eigensolution process has finished.

"`copy`"

:   To build the matrix explicitly, as in the previous option, but using a working copy of the matrix, that is, without modifying the original matrix $A$.

The default behavior is to build the coefficient matrix explicitly in a copy of $A$ (option "`copy`"). The user can change this as in the following example

```{code} console
$ ./program -st_type sinvert -eps_target 10 -st_ksp_type cg -st_pc_type jacobi -st_matmode shell
```

As always, the procedural equivalent is also available for specifying this option in the code of the program:

```{code} c
STSetMatMode(ST st,STMatMode mode);
```

The user must consider which approach is the most appropriate for the particular application. The different options have advantages and drawbacks. The "`shell`" approach is the simplest one but severely restricts the number of possibilities available for solving the system, in particular most of the PETSc preconditioners would not be available, including direct methods. The only preconditioners that can be used in this case are Jacobi (only if matrices $A$ and $B$ have the operation `MATOP_GET_DIAGONAL`) or a user-defined one.

The second approach ("`inplace`") can be much faster, specially in the generalized case. A more important advantage of this approach is that, in this case, the linear system solver can be combined with any of the preconditioners available in PETSc, including those which need to access internal matrix data-structures such as ILU. The main drawback is that, in the generalized problem, this approach probably makes sense only in the case that $A$ and $B$ have the same sparse pattern, because otherwise the function {external:doc}`MatAXPY`() might be inefficient. If the user knows that the pattern is the same (or a subset), then this can be specified with the function:

```{code} c
STSetMatStructure(ST st,MatStructure str);
```

Note that when the value of the shift $\sigma$ is very close to an eigenvalue, then the linear system will be ill-conditioned and using iterative methods may be problematic. On the other hand, in symmetric definite problems, the coefficient matrix will be indefinite whenever $\sigma$ is a point in the interior of the spectrum.

The third approach ("`copy`") uses more memory but avoids a potential problem that could appear in the "`inplace`" approach: the recovered matrix might be slightly different from the original one (due to roundoff).

{#sec:symm}
### Preserving the Symmetry in Generalized Eigenproblems

As mentioned in section [](#sec:defprob), some eigensolvers can exploit symmetry and compute a solution for Hermitian problems with less storage and/or computational cost than other methods. Also, symmetric solvers can be more accurate in some cases. However, in the case of generalized eigenvalue problems in which both $A$ and $B$ are symmetric, it happens that, due to the spectral transformation, symmetry is lost since none of the transformed operators $B^{-1}\!A-\sigma I$, $(A-\sigma B)^{-1}B$, etc. is symmetric (the same applies in the Hermitian case for complex matrices).

The solution proposed in SLEPc is based on selecting different kinds of inner products. Currently, we have the following choice of inner products:

-   Standard Hermitian inner product: $\langle x,y\rangle=x^*y$.

-   $B$-inner product: $\langle x,y\rangle_B=x^*\!B\,y$.

The second one can be used for preserving the symmetry in symmetric definite generalized problems, as described below. Note that $\langle x,y\rangle_B$ is a genuine inner product only if $B$ is symmetric positive definite (for the case of symmetric positive semi-definite $B$ see section [](#sec:purif)).

It can be shown that $\mathbb{R}^n$ with the $\langle x,y\rangle_B$ inner product is isomorphic to the Euclidean $n$-space $\mathbb{R}^n$ with the standard Hermitian inner product. This means that if we use $\langle x,y\rangle_B$ instead of the standard inner product, we are just changing the way lengths and angles are measured, but otherwise all the algebraic properties are maintained and therefore algorithms remain correct. What is interesting to observe is that the transformed operators such as $B^{-1}\!A$ or $(A-\sigma B)^{-1}B$ are self-adjoint with respect to $\langle x,y\rangle_B$.

```{figure} ../../_static/images/manual/svg/fig-abstr.svg
:alt: Abstraction used by SLEPc solvers
:name: fig:abstr

Abstraction used by SLEPc solvers
```

Internally, SLEPc operates with the scheme illustrated in figure [](#fig:abstr). The operations indicated by dashed arrows are implemented as virtual functions. From the user point of view, all the above explanation is transparent. The only thing he/she has to care about is to set the problem type appropriately with `EPSSetProblemType()` (see section [](#sec:defprob)). In the case of the Cayley transform, SLEPc is using $\langle x,y\rangle_{A+\nu B}$ as the inner product for preserving symmetry.

Using the $B$-inner product may be attractive also in the non-symmetric case ($A$ non-symmetric) as described in the next subsection.

Note that the above discussion is not directly applicable to `STPRECOND` and the preconditioned eigensolvers, in the sense that the goal is not to recover the symmetry of the operator. Still, the $B$-inner product is also used in generalized symmetric-definite problems.

{#sec:purif}
### Purification of Eigenvectors

In generalized eigenproblems, the case of singular $B$ deserves especial consideration. In this case the default spectral transformation (`STSHIFT`) cannot be used since $B^{-1}$ does not exist.

In shift-and-invert with operator matrix $T=(A-\sigma B)^{-1}B$, when $B$ is singular all the eigenvectors that belong to finite eigenvalues are also eigenvectors of $T$ and belong to the range of $T$, $\mathcal{R}(T)$. In this case, the bilinear form $\langle x,y\rangle_B$ introduced in section [](#sec:symm) is a semi-inner product and $\|x\|_B=\sqrt{\langle x,x\rangle_B}$ is a semi-norm. As before, $T$ is self-adjoint with respect to this inner product since $B\,T=T^*B$. Also, $\langle x,y\rangle_B$ is a true inner product on $\mathcal{R}(T)$.

The implication of all this is that, for singular $B$, if the $B$-inner product is used throughout the eigensolver then, assuming that the initial vector has been forced to lie in $\mathcal{R}(T)$, the computed eigenvectors should be correct, i.e., they should belong to $\mathcal{R}(T)$ as well. Nevertheless, finite precision arithmetic spoils this nice picture, and computed eigenvectors are easily corrupted by components of vectors in the null-space of $B$. Additional computation is required for achieving the desired property. This is usually referred to as *eigenvector purification*.

Although more elaborate purification strategies have been proposed (usually trying to reduce the computational effort, see {cite:p}`Nou87,Mee97`), the approach in SLEPc is simply to explicitly force the initial vector in the range of $T$, with $v_0\leftarrow Tv_0$, as well as the computed eigenvectors at the end, $x_i\leftarrow Tx_i$. Since this computation can be costly, it can be deactivated if the user knows that $B$ is non-singular, with:

```{code} c
EPSSetPurify(EPS eps,PetscBool purify);
```

A final comment is that eigenvector corruption happens also in the non-symmetric case. If $A$ is non-symmetric but $B$ is symmetric positive semi-definite, then the scheme presented above ($B$-inner product together with purification) can still be applied and is generally more successful than the straightforward approach with the standard inner product. For using this scheme in SLEPc, the user has to specify the special problem type `EPS_PGNHEP`, see table [](#tab:ptype).

{#sec:slice}
### Spectrum Slicing

In the context of symmetric-definite generalized eigenvalue problems (`EPS_GHEP`) it is often required to compute all eigenvalues contained in a given interval $[a,b]$. This poses some difficulties, such as:

-   The number of eigenvalues in the interval is not known a priori.

-   There might be many eigenvalues, in some applications a significant percentage of the spectrum (20%, say).

-   We must make certain that no eigenvalues are missed, and in particular all eigenvalues must be computed with their correct multiplicity.

-   In some applications, the interval is open in one end, i.e., either $a$ or $b$ can be infinite.

One possible strategy to solve this problem is to sweep the interval from one end to the other, computing chunks of eigenvalues with a spectral transformation that updates the shift dynamically. This is generally referred to as *spectrum slicing*. The method implemented in SLEPc is similar to that proposed by {cite:t}`Gri94`, where inertia information is used to validate sub-intervals. Given a symmetric-indefinite triangular factorization

```{math}
:label: eq-symindtri

A-\sigma B=LDL^T,
```

 by Sylvester's law of inertia we know that the number of eigenvalues on the left of $\sigma$ is equal to the number of negative eigenvalues of $D$,

```{math}
:label: eq-inertia

\nu(A-\sigma B)=\nu(D).
```

 A detailed description of the method available in SLEPc can be found in {cite:p}`Cam12`. The SLEPc interface hides all the complications of the algorithm. However, the user must be aware of all the restrictions for this technique to be employed:

-   This is currently implemented only in Krylov-Schur.

-   The method is based on shift-and-invert, so `STSINVERT` must be used. Furthermore, direct linear solvers are required.

-   The direct linear solver must provide the matrix inertia (see PETSc's {external:doc}`MatGetInertia`()).

An example command-line that sets up all the required options is:

```{code} console
$ ./ex2 -n 50 -eps_interval 0.4,0.8 -st_type sinvert -st_ksp_type preonly -st_pc_type cholesky
```

Note that PETSc's Cholesky factorization is not parallel, so for doing spectrum slicing in parallel it is required to use an external solver that supports inertia. For example, with MUMPS (see section [](#sec:lin) on how to use external linear solvers) we would do:

```{code} console
$ ./ex2 -n 50 -eps_interval 0.4,0.8 -st_type sinvert -st_ksp_type preonly -st_pc_type cholesky -st_pc_factor_mat_solver_type mumps -st_mat_mumps_icntl_13 1
```

The last option is required by MUMPS to compute the inertia. An alternative is to use SuperLU_DIST, in which case the required options would be:

```{code} console
$ ./ex2 -n 50 -eps_interval 0.4,0.8 -st_type sinvert -st_ksp_type preonly -st_pc_type cholesky -st_pc_factor_mat_solver_type superlu_dist -st_mat_superlu_dist_rowperm NOROWPERM
```

In the latter example, {external:doc}`MatSetOption`() must be used in both matrices to explicitly state that they are symmetric (or Hermitian in the complex case).

Apart from the above recommendations, the following must be taken into account:

-   The parameters `nev` and `ncv` from `EPSSetDimensions()` are determined internally (user-provided values are ignored, and set to the number of eigenvalues in the interval). It is possible for the user to specify the "local" `nev` and `ncv` parameters that will be used when computing eigenvalues around each shift, with `EPSKrylovSchurSetDimensions()`.

-   The user can also tune the computation by setting the value of `max_it`.

:::{note}
**Usage with Complex Scalars**:
Some external packages that provide inertia information (MUMPS, Pardiso) do so only in real scalars, but not in the case of complex scalars. Hence, with complex scalars spectrum slicing is available only sequentially (with PETSc's Cholesky factorization) or via SuperLU_DIST (as in the last example above). An alternative to spectrum slicing is to use the CISS solver with a region enclosing an interval on the real axis, see {cite:p}`Mae16` for details.
:::

#### Use of Multiple Communicators

Since spectrum slicing requires direct linear solves, parallel runs may suffer from bad scalability in the sense that increasing the number of MPI processes does not imply a performance gain. For this reason, SLEPc provides the option of using multiple communicators, that is, splitting the initial MPI communicator in several groups, each of them in charge of processing part of the interval.

The multi-communicator setting is activated with a value of `npart`\>1 in the following function:

```{code} c
EPSKrylovSchurSetPartitions(EPS eps,PetscInt npart);
```

The interval $[a,b]$ is then divided in `npart` subintervals of equal size, and the problem of computing all eigenvalues in $[a,b]$ is divided in `npart` independent subproblems. Each subproblem is solved using only a subset of the initial $p$ processes, with $\lceil p/\texttt{npart}\rceil$ processes at most. A final step will gather all computed solutions so that they are available in the whole `EPS` communicator.

The division of the interval in subintervals is done blindly, and this may result in load imbalance if some subintervals contain much more eigenvalues than others. This can be prevented by passing a list of subinterval boundaries, provided that the user has a priori information to roughly determine the eigenvalue distribution:

```{code} c
EPSKrylovSchurSetSubintervals(EPS eps,PetscReal *subint);
```

An additional benefit of multi-communicator support is that it enables parallel spectrum slicing runs without the need to install a parallel direct solver (MUMPS), by setting the number of partitions equal to the number of MPI processes. The following command-line example uses sequential linear solves in 4 partitions, one process each:

```{code} console
$ mpiexec -n 4 ./ex25 -eps_interval 0.4,0.8 -eps_krylovschur_partitions 4 -st_type sinvert -st_ksp_type preonly -st_pc_type cholesky
```

The analog example using MUMPS with 5 processes in each partition:

```{code} console
$ mpiexec -n 20 ./ex25 -eps_interval 0.4,0.8 -eps_krylovschur_partitions 4 -st_type sinvert -st_ksp_type preonly -st_pc_type cholesky -st_pc_factor_mat_solver_type mumps -st_mat_mumps_icntl_13 1
```

### Spectrum Folding

In SLEPc versions prior to 3.5, `ST` had another type intended to perform the spectrum folding technique described below. It is no longer available with `ST`, but it can be implemented directly in application code as illustrated in example `ex24.c`.

Spectrum folding involves squaring in addition to shifting. This makes sense for standard Hermitian eigenvalue problems, where the transformed problem to be addressed is

```{math}
:label: eq:spectrum-folding

(A-\sigma I)^2x=\theta x.
```

The following relation holds

```{math}
:label: eq:spectrum-folding-relation

\theta=(\lambda-\sigma)^2.
```

Note that the mapping between $\lambda$ and $\theta$ is not injective, and hence this cannot be considered a true spectral transformation.

```{figure} ../../_static/images/manual/svg/fig-fold.svg
:alt: Illustration of the effect of spectrum folding
:name: fig:fold

Illustration of the effect of spectrum folding
```

The effect is that the spectrum is folded around the value of $\sigma$. Thus, eigenvalues that are closest to the shift become the smallest eigenvalues in the folded spectrum, see figure [](#fig:fold). For this reason, spectrum folding is commonly used in combination with eigensolvers that compute the smallest eigenvalues, for instance in the context of electronic structure calculations {cite:p}`Can00`. This transformation can be an effective, low-cost alternative to shift-and-invert.

```{only} html
<p class="rubric">References</p>
```
```{bibliography}
:filter: docname in docnames
```

```{rubric} Footnotes
```

[^slepc-v3.5]: Note that the sign changed in SLEPc 3.5 with respect to previous versions.
