(ch:pep)=
# PEP: Polynomial Eigenvalue Problems

The Polynomial Eigenvalue Problem (`PEP`) solver object is intended for addressing polynomial eigenproblems of arbitrary degree, $P(\lambda)x=0$. A particular instance is the quadratic eigenvalue problem (degree 2), which is the case more often encountered in practice. For this reason, part of the description of this chapter focuses specifically on quadratic eigenproblems.

Currently, most `PEP` solvers are based on linearization, either implicit or explicit. The case of explicit linearization allows the use of eigensolvers from `EPS` to solve the linearized problem.

{#sec:pep label="sec:pep"}
## Overview of Polynomial Eigenproblems

In this section, we review some basic properties of the polynomial eigenvalue problem. The main goal is to set up the notation as well as to describe the linearization approaches that will be employed for solving via an `EPS` object. To simplify, we initially restrict the description to the case of quadratic eigenproblems, and then extend it to the general case of arbitrary degree. For additional background material, the reader is referred to {cite:p}`Tis01`. More information can be found in a paper by {cite:t}`Cam16a`, that focuses specifically on SLEPc implementation of methods based on Krylov iterations on the linearized problem.

{#sec:qep label="sec:qep"}
### Quadratic Eigenvalue Problems

In many applications, e.g., problems arising from second-order differential equations such as the analysis of damped vibrating systems, the eigenproblem to be solved is quadratic,

```{math}
:label: eq:eigquad

(K+\lambda C+\lambda^2M)x=0,
```

where $K,C,M\in\mathbb{C}^{n\times n}$ are the coefficients of a polynomial matrix of degree 2, $\lambda\in\mathbb{C}$ is the eigenvalue and $x\in\mathbb{C}^n$ is the eigenvector. As in the case of linear eigenproblems, the eigenvalues and eigenvectors can be complex even in the case that all three matrices are real.

It is important to point out some outstanding differences with respect to the linear eigenproblem. In the quadratic eigenproblem, the number of eigenvalues is $2n$, and the corresponding eigenvectors do not form a linearly independent set. If $M$ is singular, some eigenvalues are infinite. Even when the three matrices are symmetric and positive definite, there is no guarantee that the eigenvalues are real, but still methods can exploit symmetry to some extent. Furthermore, numerical difficulties are more likely than in the linear case, so the computed solution can sometimes be untrustworthy.

If equation {math:numref}`eq:eigquad` is written as $P(\lambda)x=0$, where $P$ is the polynomial matrix, then multiplication by $\lambda^{-2}$ results in $\operatorname{rev} P(\lambda^{-1})x=0$, where $\operatorname{rev} P$ denotes the polynomial matrix with the coefficients of $P$ in the reverse order. In other words, if a method is available for computing the largest eigenvalues, then reversing the roles of $M$ and $K$ results in the computation of the smallest eigenvalues. In general, it is also possible to formulate a spectral transformation for computing eigenvalues closest to a given target, as discussed in section [](#sec:qst).

#### Problem Types

As in the case of linear eigenproblems, there are some particular properties of the coefficient matrices that confer a certain structure to the quadratic eigenproblem, e.g., symmetry of the spectrum with respect to the real or imaginary axes. These structures are important as long as the solvers are able to exploit them.

-   Hermitian (symmetric) problems, when $M$, $C$, $K$ are all Hermitian (symmetric). Eigenvalues are real or come in complex conjugate pairs. Furthermore, if $M>0$ and $C,K\geq 0$ then the system is stable, i.e., $\text{Re}(\lambda)\leq 0$.

-   Hyperbolic problems, a particular class of Hermitian problems where $M>0$ and $(x^*Cx)^2>4(x^*Mx)(x^*Kx)$ for all nonzero $x\in\mathbb{C}^n$. All eigenvalues are real, and form two separate groups of $n$ eigenvalues, each of them having linearly independent eigenvectors. A particular subset of hyperbolic problems is the class of overdamped problems, where $C>0$ and $K\geq 0$, in which case all eigenvalues are non-positive.

-   Gyroscopic problems, when $M$, $K$ are Hermitian, $M>0$, and $C$ is skew-Hermitian, $C=-C^*$. The spectrum is symmetric with respect to the imaginary axis, and in the real case, it has a Hamiltonian structure, i.e., eigenvalues come in quadruples $(\lambda,\bar{\lambda},-\lambda,-\bar{\lambda})$.

:::{note}
Currently, the problem type is not exploited by `PEP` solvers, except for a few exceptions. In the future, we may add more support for structure-preserving solvers.
:::

{#sec:linearization}
#### Linearization

It is possible to transform the quadratic eigenvalue problem to a linear generalized eigenproblem $L_0y=\lambda L_1y$ by doubling the order of the system, i.e., $L_0,L_1\in\mathbb{C}^{2n\times 2n}$. There are many ways of doing this. For instance, consider the following two pencils $L(\lambda)=L_0-\lambda L_1$,

```{math}
:label: eq:n1

\begin{bmatrix}0 & I\\-K & -C\end{bmatrix}-\lambda\begin{bmatrix}I & 0\\0 & M\end{bmatrix},
```

```{math}
:label: eq:n2

\begin{bmatrix}-K & 0\\0 & I\end{bmatrix}-\lambda\begin{bmatrix}C & M\\I & 0\end{bmatrix}.
```

 Both of them have the same eigenvalues as the quadratic eigenproblem, and the corresponding eigenvectors can be expressed as

```{math}
:label: eq:linevec

y=\begin{bmatrix}x\\x\lambda\end{bmatrix},
```

 where $x$ is the eigenvector of the quadratic eigenproblem.

Other **non-symmetric** linearizations can be obtained by a linear combination of equations {math:numref}`eq:n1` and {math:numref}`eq:n2`,

```{math}
:label: eq:lingen

\begin{bmatrix}-\beta K & \alpha I\\-\alpha K & -\alpha C+\beta I\end{bmatrix}-\lambda\begin{bmatrix}\alpha I+\beta C & \beta M\\\beta I & \alpha M\end{bmatrix}.
```

 for any $\alpha,\beta\in\mathbb{R}$. The linearizations {math:numref}`eq:n1` and {math:numref}`eq:n2` are particular cases of equation {math:numref}`eq:lingen` taking $(\alpha,\beta)=(1,0)$ and $(0,1)$, respectively.

**Symmetric** linearizations are useful for the case that $M$, $C$, and $K$ are all symmetric (Hermitian), because the resulting matrix pencil is symmetric (Hermitian), although indefinite:

```{math}
:label: eq:linsym

\begin{bmatrix}\beta K&\alpha K\\\alpha K&\alpha C-\beta M\end{bmatrix}-\lambda\begin{bmatrix}\alpha K-\beta C&-\beta M\\-\beta M&-\alpha M\end{bmatrix},
```

And for gyroscopic problems, we can consider **Hamiltonian** linearizations,

```{math}
:label: eq:linham

\begin{bmatrix}\alpha K & -\beta K\\\alpha C+\beta M & \alpha K\end{bmatrix}-\lambda\begin{bmatrix}\beta M & \alpha K+\beta C\\-\alpha M & \beta M\end{bmatrix},
```

 where one of the matrices is Hamiltonian and the other one is skew-Hamiltonian if $(\alpha,\beta)$ is $(1,0)$ or $(0,1)$.

In SLEPc, the `PEPLINEAR` solver is based on using one of the above linearizations for solving the quadratic eigenproblem. This solver makes use of linear eigensolvers from the `EPS` package.

We could also consider the *reversed* forms, e.g., the reversed form of equation {math:numref}`eq:n2` is

```{math}
:label: eq:n2r

\begin{bmatrix}-C & -M\\I & 0\end{bmatrix}-\frac{1}{\lambda}\begin{bmatrix}K & 0\\0 & I\end{bmatrix},
```

 which is equivalent to the form {math:numref}`eq:n1` for the problem $\operatorname{rev} P(\lambda^{-1})x=0$. These reversed forms are not implemented in SLEPc, but the user can use them simply by reversing the roles of $M$ and $K$, and considering the reciprocals of the computed eigenvalues. Alternatively, this can be viewed as a particular case of the spectral transformation (with $\sigma=0$), see section [](#sec:qst).

{#sec:pep1 label="sec:pep1"}
### Polynomials of Arbitrary Degree

In general, the polynomial eigenvalue problem can be formulated as

```{math}
:label: eq:pep

P(\lambda)x=0,
```

where $P$ is an $n\times n$ polynomial matrix of degree $d$. An $n$-vector $x\neq 0$ satisfying this equation is called an eigenvector associated with the corresponding eigenvalue $\lambda$.

We start by considering the case where $P$ is expressed in terms of the monomial basis,

```{math}
:label: eq:pepmon

P(\lambda)=A_0+A_1 \lambda+A_2\lambda^2 +  \dotsb + A_d \lambda^d,
```

 where $A_0,\ldots,A_d$ are the $n\times n$ coefficient matrices. As before, the problem can be solved via some kind of linearization. One of the most commonly used ones is the first companion form

```{math}
:label: eq:firstcomp

L(\lambda)=L_0 -\lambda L_1,
```

 where the related linear eigenproblem is $L(\lambda)y=0$, with

```{math}
:label: eq:firstcompfull

L_0 =
\begin{bmatrix}
  & I \\
  & & \ddots \\
  & & & I \\
  -A_0 & -A_1 & \cdots  & -A_{d-1}
\end{bmatrix},\quad
L_1 =
\begin{bmatrix}
  I \\
  & \ddots \\
  & & I \\
  & & & A_d
\end{bmatrix}, \quad
y=
\begin{bmatrix}
  x \\ x\lambda\\ \vdots \\ x\lambda^{d-1}
\end{bmatrix}.
```

 This is the generalization of equation {math:numref}`eq:n1`.

The definition of vector $y$ above contains the successive powers of $\lambda$. For large polynomial degree, these values may produce overflow in finite precision computations, or at least lead to numerical instability of the algorithms due to the wide difference in magnitude of the eigenvector entries. For this reason, it is generally recommended to work with non-monomial polynomial bases whenever the degree is not small, e.g., for $d>5$.

In the most general formulation of the polynomial eigenvalue problem, $P$ is expressed as

```{math}
:label: eq:pepnonmon

P(\lambda)=A_0\phi_0(\lambda)+A_1\phi_1(\lambda)+\dots+A_d\phi_d(\lambda),
```

 where $\phi_i$ are the members of a given polynomial basis, for instance, some kind of orthogonal polynomials such as Chebyshev polynomials of the first kind. In that case, the expression of $y$ in equation {math:numref}`eq:firstcompfull` contains $\phi_0(\lambda),\dots,\phi_d(\lambda)$ instead of the powers of $\lambda$. Correspondingly, the form of $L_0$ and $L_1$ is different for each type of polynomial basis.

#### Avoiding the Linearization

An alternative to linearization is to directly perform a projection of the polynomial eigenproblem. These methods enforce a Galerkin condition on the polynomial residual, $P(\theta)u\perp \mathcal{K}$. Here, the subspace $\mathcal{K}$ can be built in various ways, for instance with the Jacobi-Davidson method. This family of methods need not worry about operating with vectors of dimension $dn$. The downside is that computing more than one eigenvalue is more difficult, since usual deflation strategies cannot be applied. For a detailed description of the polynomial Jacobi-Davidson method in SLEPc, see {cite:p}`Cam20a`.

## Basic Usage

The user interface of the `PEP` package is very similar to `EPS`. For basic usage, the most noteworthy difference is that all coefficient matrices $A_i$ have to be supplied in the form of an array of {external:doc}`Mat`.

A basic example code for solving a polynomial eigenproblem with `PEP` is shown in listing [](#fig:ex-pep), where the code for building matrices `A[0]`, `A[1]`, ... is omitted. The required steps are the same as those described in chapter [](#ch:eps) for the linear eigenproblem. As always, the solver context is created with `PEPCreate()`. The coefficient matrices are provided with `PEPSetOperators()`, and the problem type is specified with `PEPSetProblemType()`. Calling `PEPSetFromOptions()` allows the user to set up various options through the command line. The call to `PEPSolve()` invokes the actual solver. Then, the solution is retrieved with `PEPGetConverged()` and `PEPGetEigenpair()`. Finally, `PEPDestroy()` destroys the object.

```{code-block} c
:name: fig:ex-pep
:caption: Example code for basic solution with `PEP`

#define NMAT 5
PEP         pep;       /*  eigensolver context  */
Mat         A[NMAT];   /*  coefficient matrices */
Vec         xr, xi;    /*  eigenvector, x       */
PetscScalar kr, ki;    /*  eigenvalue, k        */
PetscInt    j, nconv;
PetscReal   error;

PEPCreate(PETSC_COMM_WORLD, &pep);
PEPSetOperators(pep, NMAT, A);
PEPSetProblemType(pep, PEP_GENERAL);  /* optional */
PEPSetFromOptions(pep);
PEPSolve(pep);
PEPGetConverged(pep, &nconv);
for (j=0;j<nconv;j++) {
  PEPGetEigenpair(pep, j, &kr, &ki, xr, xi);
  PEPComputeError(pep, j, PEP_ERROR_BACKWARD, &error);
}
PEPDestroy(&pep);
```

## Defining the Problem

:::{table} Polynomial bases available to represent the polynomial matrix in `PEP`
:name: tab:pepbasis

 | Polynomial Basis      | `PEPBasis`              | Options Database |
 |-----------------------|-------------------------|------------------|
 | Monomial              | `PEP_BASIS_MONOMIAL`    | `monomial`       |
 | Chebyshev (1st kind)  | `PEP_BASIS_CHEBYSHEV1`  | `chebyshev1`     |
 | Chebyshev (2nd kind)  | `PEP_BASIS_CHEBYSHEV2`  | `chebyshev2`     |
 | Legendre              | `PEP_BASIS_LEGENDRE`    | `legendre`       |
 | Laguerre              | `PEP_BASIS_LAGUERRE`    | `laguerre`       |
 | Hermite               | `PEP_BASIS_HERMITE`     | `hermite`        |

:::

As explained in section [](#sec:pep1), the polynomial matrix $P(\lambda)$ can be expressed in terms of the monomials $1$, $\lambda$, $\lambda^2,\ldots$, or in a non-monomial basis as in equation {math:numref}`eq:pepnonmon`. Hence, when defining the problem we must indicate which is the polynomial basis to be used as well as the coefficient matrices $A_i$ in that basis representation. By default, a monomial basis is used. Other possible bases are listed in table [](#tab:pepbasis), and can be set with:

```{code} c
PEPSetBasis(PEP pep,PEPBasis basis);
```

or with the command-line key `-pep_basis <name>`. The matrices are passed with:

```{code} c
PEPSetOperators(PEP pep,PetscInt nmat,Mat A[]);
```

:::{table} Problem types considered in `PEP`
:name: tab:ptypeq

 |Problem Type  |`PEPProblemType`  |Command line key
 |--------------|------------------|-------------------
 |General       |`PEP_GENERAL`     |`-pep_general`
 |Hermitian     |`PEP_HERMITIAN`   |`-pep_hermitian`
 |Hyperbolic    |`PEP_HYPERBOLIC`  |`-pep_hyperbolic`
 |Gyroscopic    |`PEP_GYROSCOPIC`  |`-pep_gyroscopic`

:::

As mentioned in section [](#sec:qep), it is possible to distinguish among different problem types. The problem types currently supported for `PEP` are listed in table [](#tab:ptypeq). The goal when choosing an appropriate problem type is to let the solver exploit the underlying structure, in order to possibly compute the solution more accurately with less floating-point operations. When in doubt, use the default problem type (`PEP_GENERAL`).

The problem type can be specified at run time with the corresponding command line key or, more usually, within the program with the function:

```{code} c
PEPSetProblemType(PEP pep,PEPProblemType type);
```

:::{note}
Currently, the problem type is ignored in most solvers and it is taken into account only in some cases for the quadratic eigenproblem only.
:::

Apart from the polynomial basis and the problem type, the definition of the problem is completed with the number and location of the eigenvalues to compute. This is done very much like in `EPS`, but with minor differences.

The number of eigenvalues (and eigenvectors) to compute, `nev`, is specified with the function:

```{code} c
PEPSetDimensions(PEP pep,PetscInt nev,PetscInt ncv,PetscInt mpd);
```

The default is to compute only one. This function also allows control over the dimension of the subspaces used internally. The second argument, `ncv`, is the number of column vectors to be used by the solution algorithm, that is, the largest dimension of the working subspace. The third argument, `mpd`, is the maximum projected dimension. These parameters can also be set from the command line with `-pep_nev`, `-pep_ncv` and `-pep_mpd`.

For the selection of the portion of the spectrum of interest, there are several alternatives listed in table [](#tab:portionq), to be selected with:

```{code} c
PEPSetWhichEigenpairs(PEP pep,PEPWhich which);
```

The default is to compute the largest magnitude eigenvalues. For the sorting criteria relative to a target value, the scalar $\tau$ must be specified with:

```{code} c
PEPSetTarget(PEP pep,PetscScalar target);
```

or in the command-line with `-pep_target`. As in `EPS`, complex values of $\tau$ are allowed only in complex scalar SLEPc builds. The criteria relative to a target must be used in combination with a spectral transformation as explained in section [](#sec:qst).

There is also support for spectrum slicing, that is, computing all eigenvalues in a given interval, see section [](#sec:qslice). For this, the user has to specify the computational interval with:

```{code} c
PEPSetInterval(PEP pep,PetscReal a,PetscReal b);
```

or equivalently with `-pep_interval a,b`.

Finally, we mention that the use of regions for filtering is also available in `PEP`, see section [](#sec:region).

:::{table} Available possibilities for selection of the eigenvalues of interest in `PEP`
:name: tab:portionq

 |`PEPWhich`                |Command line key           |Sorting criterion
 |--------------------------|---------------------------|----------------------------------------
 |`PEP_LARGEST_MAGNITUDE`   |`-pep_largest_magnitude`   |Largest $\|\lambda\|$
 |`PEP_SMALLEST_MAGNITUDE`  |`-pep_smallest_magnitude`  |Smallest $\|\lambda\|$
 |`PEP_LARGEST_REAL`        |`-pep_largest_real`        |Largest $\mathrm{Re}(\lambda)$
 |`PEP_SMALLEST_REAL`       |`-pep_smallest_real`       |Smallest $\mathrm{Re}(\lambda)$
 |`PEP_LARGEST_IMAGINARY`   |`-pep_largest_imaginary`   |Largest $\mathrm{Im}(\lambda)$[^eps-real]
 |`PEP_SMALLEST_IMAGINARY`  |`-pep_smallest_imaginary`  |Smallest $\mathrm{Im}(\lambda)$[^eps-real]
 |`PEP_TARGET_MAGNITUDE`    |`-pep_target_magnitude`    |Smallest $\|\lambda-\tau\|$
 |`PEP_TARGET_REAL`         |`-pep_target_real`         |Smallest $\|\mathrm{Re}(\lambda-\tau)\|$
 |`PEP_TARGET_IMAGINARY`    |`-pep_target_imaginary`    |Smallest $\|\mathrm{Im}(\lambda-\tau)\|$
 |`PEP_ALL`                 |`-pep_all`                 |All $\lambda\in[a,b]$
 |`PEP_WHICH_USER`          |` `                        |*user-defined*

:::

[^eps-real]: If SLEPc is compiled for real scalars, then the absolute value of the imaginary part, $\|\mathrm{Im}(\lambda)\|$, is used for eigenvalue selection and sorting.

## Selecting the Solver

The solution method can be specified procedurally with:

```{code} c
PEPSetType(PEP pep,PEPType method);
```

or via the options database command `-pep_type` followed by the name of the method. The methods currently available in `PEP` are listed in table [](#tab:solversp). All solvers except `PEPJD` are based on the linearization explained above, whereas `PEPJD` performs a projection on the polynomial problem (without linearizing).

The default solver is `PEPTOAR`. The Two-level Orthogonal Arnoldi (TOAR) method is a stable algorithm for building an Arnoldi factorization of the linearization {math:numref}`eq:firstcomp` without explicitly creating matrices $L_0,L_1$, and represents the Krylov basis in a compact way. Symmetric TOAR (STOAR) is a variant of TOAR that exploits symmetry (requires `PEP_HERMITIAN` or `PEP_HYPERBOLIC` problem types). Quadratic Arnoldi (Q-Arnoldi) is related to TOAR and follows a similar approach.

:::{table} Polynomial eigenvalue solvers available in the `PEP` module
:name: tab:solversp

 | Method                  | `PEPType`     | Options Database | Polynomial Degree | Polynomial Basis
 |-------------------------|---------------|------------------|-------------------|-----------------
 | TOAR                    | `PEPTOAR`     |`toar`            | Arbitrary         | Any
 | Symmetric TOAR          | `PEPSTOAR`    |`stoar`           | Quadratic         | Monomial
 | Q-Arnoldi               | `PEPQARNOLDI` |`qarnoldi`        | Quadratic         | Monomial
 | Linearization via `EPS` | `PEPLINEAR`   |`linear`          | Arbitrary         | Any
 | Jacobi-Davidson         | `PEPJD`       |`jd`              | Arbitrary         | Monomial
 | Contour integral SS     | `PEPCISS`     |`ciss`            | Arbitrary         | Monomial

:::

The `PEPLINEAR` method carries out an explicit linearization of the polynomial eigenproblem, as described in section [](#sec:pep), resulting in a generalized eigenvalue problem that is handled by an `EPS` object created internally. If required, this `EPS` object can be extracted with the operation:

```{code} c
PEPLinearGetEPS(PEP pep,EPS *eps);
```

This allows the application programmer to set any of the `EPS` options directly within the code. Also, it is possible to change the `EPS` options through the command-line, simply by prefixing the `EPS` options with `-pep_linear_`.

In `PEPLINEAR`, if the eigenproblem is quadratic, the expression used in the linearization is dictated by the problem type set with `PEPProblemType`, which chooses from non-symmetric {math:numref}`eq:lingen`, symmetric {math:numref}`eq:linsym`, and Hamiltonian {math:numref}`eq:linham` linearizations. The parameters $(\alpha,\beta)$ of these linearizations can be set with:

```{code} c
PEPLinearSetLinearization(PEP pep,PetscReal alpha,PetscReal beta);
```

For polynomial eigenproblems with degree $d>2$ the linearization is the one described in section [](#sec:pep1).

Another option of the `PEPLINEAR` solver is whether the matrices of the linearized problem are created explicitly or not. This is set with:

```{code} c
PEPLinearSetExplicitMatrix(PEP pep,PetscBool explicit);
```

The explicit matrix option is available only for quadratic eigenproblems (higher degree polynomials are always handled implicitly). In the case of explicit creation, matrices $L_0$ and $L_1$ are created as true {external:doc}`Mat`'s, with explicit storage, whereas the implicit option works with *shell* {external:doc}`Mat`'s that operate only with the constituent blocks $M$, $C$ and $K$ (or $A_i$ in the general case). The explicit case requires more memory but gives more flexibility, e.g., for choosing a preconditioner. Some examples of usage via the command line are shown at the end of next section.

{#sec:qst label="sec:qst"}
## Spectral Transformation for PEP

For computing eigenvalues in the interior of the spectrum (closest to a target $\tau$), it is necessary to use a spectral transformation. In `PEP` solvers this is handled via an `ST` object as in the case of linear eigensolvers. It is possible to proceed with no spectral transformation (shift) or with shift-and-invert. Every `PEP` object has an `ST` object internally.

The spectral transformation can be applied either to the polynomial problem or its linearization. We illustrate it first for the quadratic case.

Given the quadratic eigenproblem in equation {math:numref}`eq:eigquad`, it is possible to define the transformed problem

```{math}
:label: eq:sinvquad

(K_\sigma+\theta C_\sigma+\theta^2M_\sigma)x=0,
```

 where the coefficient matrices are

```{math}
:label: eq:coef-matrices

K_\sigma&= M,\\
C_\sigma&= C+2\sigma M,\\
M_\sigma&= \sigma^2 M+\sigma C+K,
```

and the relation between the eigenvalue of the original eigenproblem, $\lambda$, and the transformed one, $\theta$, is $\theta=(\lambda-\sigma)^{-1}$ as in the case of the linear eigenvalue problem. See chapter [](#ch:st) for additional details.

The polynomial eigenvalue problem of equation {math:numref}`eq:sinvquad` corresponds to the reversed form of the shifted polynomial, $\operatorname{rev} P(\theta)$. The extension to polynomial matrices of arbitrary degree is also possible, where the coefficients of $\operatorname{rev} P(\theta)$ have the general form

```{math}
:label: eq:sinvpep

T_k=\sum_{j=0}^{d-k}\binom{j+k}{k}\sigma^{j}A_{j+k},\qquad k=0,\ldots,d.
```

 The way this is implemented in SLEPc is that the `ST` object is in charge of computing the $T_k$ matrices, so that the `PEP` solver operates with these matrices as it would with the original $A_i$ matrices, without changing its behavior. We say that `ST` performs the transformation.

An alternative would be to apply the shift-and-invert spectral transformation to the linearization {math:numref}`eq:firstcomp` in a smart way, making the polynomial eigensolver aware of this fact so that it can exploit the block structure of the linearization. Let $S_\sigma:=(L_0-\sigma L_1)^{-1}L_1$, then when the solver needs to extend the Arnoldi basis with an operation such as $z=S_\sigma w$, a linear solve is required with the form

```{math}
:label: eq:sinvpeplin

\begin{bmatrix}
  -\sigma I  & I \\
  & -\sigma I & \ddots \\
  & & \ddots & I \\
  & & & -\sigma I & I \\
  -A_0 & -A_1 & \cdots  & -\tilde{A}_{d-2} & -\tilde{A}_{d-1}
\end{bmatrix}
\begin{bmatrix}
  z^0\\z^1\\\vdots\\z^{d-2}\\z^{d-1}
\end{bmatrix}
  =
\begin{bmatrix}
  w^0\\w^1\\\vdots\\w^{d-2}\\A_dw^{d-1}
\end{bmatrix},
```

 with $\tilde{A}_{d-2}=A_{d-2}+\sigma I$ and $\tilde{A}_{d-1}=A_{d-1}+\sigma A_d$. From the block LU factorization, it is possible to derive a simple recurrence to compute $z^i$, with one of the steps involving a linear solve with $P(\sigma)$.

Implementing the latter approach is more difficult (especially if different polynomial bases must be supported), and requires an intimate relation with the `PEP` solver. That is why it is only available currently in the default solver (TOAR) and in `PEPLINEAR` without explicit matrix. In order to choose between the two approaches, the user can set a flag with:

```{code} c
STSetTransform(ST st,PetscBool flg);
```

(or in the command line `-st_transform`) to activate the first one (`ST` performs the transformation). Note that this flag belongs to `ST`, not `PEP` (use `PEPGetST()` to extract it).

In terms of overall computational cost, both approaches are roughly equivalent, but the advantage of the second one is not having to store the $T_k$ matrices explicitly. It may also be slightly more accurate. Hence, the `STSetTransform()` flag is turned off by default.

A command line example would be:

```{code} console
$ ./ex16 -pep_nev 12 -pep_type toar -pep_target 0 -st_type sinvert
```

The example computes 12 eigenpairs closest to the origin with TOAR and shift-and-invert. The `-st_transform` could be added optionally to switch to `ST` being in charge of the transformation. The same example with Q-Arnoldi would be

```{code} console
$ ./ex16 -pep_nev 12 -pep_type qarnoldi -pep_target 0 -st_type sinvert -st_transform
```

where in this case `-st_transform` would be set as default if not specified.

As a complete example of how to solve a quadratic eigenproblem via explicit linearization with explicit construction of the $L_0$ and $L_1$ matrices, consider the following command line:

```{code} console
$ ./sleeper -pep_type linear -pep_target -10 -pep_linear_st_type sinvert -pep_linear_st_ksp_type preonly -pep_linear_st_pc_type lu -pep_linear_st_pc_factor_mat_solver_type mumps -pep_linear_st_mat_mumps_icntl_14 100 -pep_linear_explicitmatrix
```

This example uses MUMPS for solving the associated linear systems, see section [](#sec:lin) for details. The following command line example illustrates how to solve the same problem without explicitly forming the matrices. Note that in this case the `ST` options are not prefixed with `-pep_linear_` since now they do not refer to the `ST` within the `PEPLINEAR` solver but the general `ST` associated to `PEP`.

```{code} console
$ ./sleeper -pep_type linear -pep_target -10 -st_type sinvert -st_ksp_type preonly -st_pc_type lu -st_pc_factor_mat_solver_type mumps -st_mat_mumps_icntl_14 100
```

{#sec:qslice label="sec:qslice"}
### Spectrum Slicing for PEP

Similarly to the spectrum slicing technique available in linear symmetric-definite eigenvalue problems (cf. [](#sec:slice)), it is possible to compute all eigenvalues in a given interval $[a,b]$ for the case of hyperbolic quadratic eigenvalue problems (`PEP_HYPERBOLIC`). In more general symmetric (or Hermitian) quadratic eigenproblems (`PEP_HERMITIAN`), it may also be possible to do spectrum slicing provided that computing inertia is feasible, which essentially means that all eigenvalues in the interval must be real and of the same definite type.

This computation is available only in the `PEPSTOAR` solver. The spectrum slicing mechanism implemented in `PEP` is very similar to the one described in section [](#sec:slice) for linear problems, except for the multi-communicator option which is not implemented yet.

A command line example is the following:

```{code} console
$ ./spring -n 300 -pep_hermitian -pep_interval -10.1,-9.5 -pep_type stoar -st_type sinvert -st_ksp_type preonly -st_pc_type cholesky
```

In hyperbolic problems, where eigenvalues form two separate groups of $n$ eigenvalues, it will be necessary to explicitly set the problem type to `-pep_hyperbolic` if the interval $[a,b]$ includes eigenvalues from both groups.

Additional details can be found in {cite:p}`Cam20b`.

## Retrieving the Solution

After the call to `PEPSolve()` has finished, the computed results are stored internally. The procedure for retrieving the computed solution is exactly the same as in the case of `EPS`. The user has to call `PEPGetConverged()` first, to obtain the number of converged solutions, then call `PEPGetEigenpair()` repeatedly within a loop, once per each eigenvalue-eigenvector pair. The same considerations relative to complex eigenvalues apply, see section [](#sec:retrsol) for additional details.

**Controlling and Monitoring Convergence**:
As in the case of `EPS`, in `PEP` the number of iterations carried out by the solver can be determined with `PEPGetIterationNumber()`, and the tolerance and maximum number of iterations can be set with `PEPSetTolerances()`. Also, convergence can be monitored with command-line keys `-pep_monitor`, `-pep_monitor_all`, `-pep_monitor_conv`, `-pep_monitor draw::draw_lg`, or `-pep_monitor_all draw::draw_lg`. See section [](#sec:monitor) for additional details.

**Viewing the Solution**:
Likewise to linear eigensolvers, there is support for various kinds of viewers for the solution. One can for instance use `-pep_view_values`, `-pep_view_vectors`, `-pep_error_relative`, or `-pep_converged_reason`. See description in section [](#sec:epsviewers).

### Reliability of the Computed Solution

:::{table} Possible expressions for computing error bounds
:name: tab:peperrors

 |Error type      |`PEPErrorType`        |Command line key       |Error bound
 |----------------|----------------------|-----------------------|--------------------------------------
 |Absolute error  |`PEP_ERROR_ABSOLUTE`  |`-pep_error_absolute`  |$\\|r\\|$
 |Relative error  |`PEP_ERROR_RELATIVE`  |`-pep_error_relative`  |$\\|r\\|/\|\lambda\|$
 |Backward error  |`PEP_ERROR_BACKWARD`  |`-pep_error_backward`  |$\\|r\\|/(\sum_j\\|A_j\\|\|\lambda_i\|^j)$

:::

As in the case of linear problems, the following function:

```{code} c
PEPComputeError(PEP pep,PetscInt j,PEPErrorType type,PetscReal *error);
```

is available to assess the accuracy of the computed solutions. This error is based on the computation of the 2-norm of the residual vector, defined as

```{math}
:label: eq:respol

r=P(\tilde{\lambda})\tilde{x},
```

 where $\tilde{\lambda}$ and $\tilde{x}$ represent any of the `nconv` computed eigenpairs delivered by `PEPGetEigenpair()`. From the residual norm, the error bound can be computed in different ways, see table [](#tab:peperrors). It is usually recommended to assess the accuracy of the solution using the backward error, defined as

```{math}
:label: eq:backward

\eta(\tilde{\lambda},\tilde{x})=\frac{\|r\|}{\sum_{j=0}^d\|A_j\||\tilde\lambda|^j\|\tilde{x}\|},
```

 where $d$ is the degree of the polynomial. Note that the eigenvector is always assumed to have unit norm.

Similar expressions can be used in the convergence criterion used to accept converged eigenpairs internally by the solver. The convergence test can be set via the corresponding command-line switch (see table [](#tab:pepconv)) or with `PEPSetConvergenceTest()`.

:::{table} Available possibilities for the convergence criterion
:name: tab:pepconv

 |Convergence criterion     |`PEPConv`        |Command line key  |Error bound
 |--------------------------|-----------------|------------------|--------------------------------------
 |Absolute                  |`PEP_CONV_ABS`   |`-pep_conv_abs`   |$\\|r\\|$
 |Relative to eigenvalue    |`PEP_CONV_REL`   |`-pep_conv_rel`   |$\\|r\\|/\|\lambda\|$
 |Relative to matrix norms  |`PEP_CONV_NORM`  |`-pep_conv_norm`  |$\\|r\\|/(\sum_j\\|A_j\\|\|\lambda\|^j)$
 |User-defined              |`PEP_CONV_USER`  |`-pep_conv_user`  |*user function*

:::

{#sec:scaling label="sec:scaling"}
### Scaling

When solving a quadratic eigenproblem via linearization, an accurate solution of the generalized eigenproblem does not necessarily imply a similar level of accuracy for the quadratic problem. {cite:t}`Tis00` shows that in the case of the linearization {math:numref}`eq:n1`, a small backward error in the generalized eigenproblem guarantees a small backward error in the quadratic eigenproblem. However, this holds only if $M$, $C$ and $K$ have a similar norm.

When the norm of $M$, $C$ and $K$ vary widely, {cite:t}`Tis00` recommends to solve the scaled problem, defined as

```{math}
:label: eq:scaled

(\mu^2M_\alpha+\mu C_\alpha+K)x=0,
```

 with $\mu=\lambda/\alpha$, $M_\alpha=\alpha^2M$ and $C_\alpha=\alpha C$, where $\alpha$ is a scaling factor. Ideally, $\alpha$ should be chosen in such a way that the norms of $M_\alpha$, $C_\alpha$ and $K$ have similar magnitude. A tentative value would be $\alpha=\sqrt{\frac{\|K\|_\infty}{\|M\|_\infty}}$.

In the general case of polynomials of arbitrary degree, a similar scheme is also possible, but it is not clear how to choose $\alpha$ to achieve the same goal. {cite:t}`Bet08` proposes such a scaling scheme as well as more general diagonal scalings $D_\ell P(\lambda)D_r$. In SLEPc, we provide these types of scalings, whose settings can be tuned with:

```{code} c
PEPSetScale(PEP pep,PEPScale scale,PetscReal alpha,Vec Dl,Vec Dr,PetscInt its,PetscReal lambda);
```

See the manual page for details and the description in {cite:p}`Cam16a`.

{#sec:pepextr label="sec:pepextr"}
### Extraction

Some of the eigensolvers provided in the `PEP` package are based on solving the linearized eigenproblem of equation {math:numref}`eq:firstcompfull`. From the eigenvector $y$ of the linearization, it is possible to extract the eigenvector $x$ of the polynomial eigenproblem. The most straightforward way is to take the first block of $y$, but there are other, more elaborate extraction strategies. For instance, one may compute the norm of the residual {math:numref}`eq:respol` for every block of $y$, and take the one that gives the smallest residual. The different extraction techniques may be selected with:

```{code} c
PEPSetExtract(PEP pep,PEPExtract extract);
```

For additional information, see {cite:p}`Cam16a`.

{#sec:refine label="sec:refine"}
### Iterative Refinement

As mentioned above, scaling can sometimes improve the accuracy of the computed solution considerably, in the case that the coefficient matrices $A_i$ are very different in norm. Still, even when the matrix norms are well balanced the accuracy can sometimes be unacceptably low. The reason is that methods based on linearization are not always backward stable, that is, even if the computation of the eigenpairs of the linearization is done in a stable way, there is no guarantee that the extracted polynomial eigenpairs satisfy the given tolerance.

If good accuracy is required, one possibility is to perform a few steps of iterative refinement on the solution computed by the polynomial eigensolver algorithm. Iterative refinement can be seen as the Newton method applied to a set of nonlinear equations related to the polynomial eigenvalue problem {cite:p}`Bet11`. It is well known that global convergence of Newton's iteration is guaranteed only if the initial guess is close enough to the exact solution, so we still need an eigensolver such as TOAR to compute this initial guess.

Iterative refinement can be very costly (sometimes a single refinement step is more expensive than the whole iteration to compute the initial guess with TOAR), that is why in SLEPc it is disabled by default. When the user activates it, the computation of Newton iterations will take place within `PEPSolve()` as a final stage (identified as `PEPRefine()` in the `-log_view` report).

```{code} c
PEPSetRefine(PEP pep,PEPRefine refine,PetscInt npart,PetscReal tol,PetscInt its,PEPRefineScheme scheme);
```

There are two types of refinement, identified as *simple* and *multiple*. The first one performs refinement on each eigenpair individually, while the second one considers the computed invariant pair as a whole. This latter approach is more costly but it is expected to be more robust in the presence of multiple eigenvalues.

In `PEPSetRefine()`, the argument `npart` indicates the number of partitions in which the communicator must be split. This can sometimes improve the scalability when refining many eigenpairs.

Additional details can be found in {cite:p}`Cam16b`.

```{only} html
<p class="rubric">References</p>
```
```{bibliography}
:filter: docname in docnames
```

```{rubric} Footnotes
```
