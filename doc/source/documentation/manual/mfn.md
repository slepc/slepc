(ch:mfn)=
# MFN: Matrix Function

The Matrix Function (`MFN`) solver object provides algorithms that compute the action of a matrix function on a given vector, without evaluating the matrix function itself. This is not an eigenvalue problem, but some methods rely on approximating eigenvalues (for instance with Krylov subspaces) and that is why we have this in SLEPc.

{#sec:mfn label="sec:mfn"}
## The Problem $f(A)v$

The need to evaluate a function $f(A)\in\mathbb{C}^{n\times n}$ of a matrix $A\in\mathbb{C}^{n\times n}$ arises in many applications. There are many methods to compute matrix functions, see for instance the survey by {cite:t}`Higham:2010:CMF`. Here, we focus on the case that $A$ is large and sparse, or is available only as a matrix-vector product subroutine. In such cases, it is the action of $f(A)$ on a vector, $f(A)v$, that is required and not $f(A)$. For this, it is possible to adapt some of the methods used to approximate eigenvalues, such as those based on Krylov subspaces or on the concept of contour integral. The description below will be restricted to the case of Krylov methods.

In the sequel, we concentrate on the exponential function, which is one of the most demanded in applications, although the concepts are easily generalizable to other functions as well. Using the Taylor series expansion of $e^A$, we have

```{math}
:label: eq:taylor-series

y=e^Av=v+\frac{A}{1!}v+\frac{A^2}{2!}v+\cdots,
```

so, in principle, the vector $y$ can be approximated by an element of the Krylov subspace $\mathcal{K}_m(A,v)$ defined in equation {math:numref}`eq:krylov`. This is the basis of the method implemented in EXPOKIt {cite:p}`Sidje:1998:ESP`. Let $AV_m=V_{m+1}\underline{H}_m$ be an Arnoldi decomposition, where the columns of $V_m$ form an orthogonal basis of the Krylov subspace, then the approximation can be computed as

```{math}
:label: eq:aprox-solution

\tilde y=\beta V_m\exp(H_m)e_1,
```

 where $\beta=\|v\|_2$ and $e_1$ is the first coordinate vector. Hence, the problem of computing the exponential of a large matrix $A$ of order $n$ is reduced to computing the exponential of a small matrix $H_m$ of order $m$. For the latter task, we employ algorithms implemented in the `FN` auxiliary class, see section [](#sec:sys).

## Basic Usage

The user interface of the `MFN` package is simpler than the interface of eigensolvers. In some ways, it is more similar to {external:doc}`KSP`, in the sense that the solver maps a vector $v$ to a vector $y$.

```{code-block} c
:name: fig:ex-mfn
:caption: Example code for basic solution with `MFN`

MFN         mfn;       /*  MFN solver context                  */
Mat         A;         /*  problem matrix                      */
FN          f;         /*  the function, exp() in this example */
PetscScalar alpha;     /*  to compute exp(alpha*A)             */
Vec         v, y;      /*  right vector and solution           */

MFNCreate( PETSC_COMM_WORLD, &mfn );
MFNSetOperator( mfn, A );
MFNGetFN( mfn, &f );
FNSetType( f, FNEXP );
FNSetScale( f, alpha, 1.0 );
MFNSetFromOptions( mfn );
MFNSolve( mfn, v, y );
MFNDestroy( &mfn );
```

Figure [](#fig:ex-mfn) shows a simple example with the basic steps for computing $y=\exp(\alpha A)v$. After creating the solver context with `MFNCreate`, the problem matrix has to be passed with `MFNSetOperator` and the function to compute $f(\cdot)$ must be specified with the aid of the auxiliary class `FN`, see details in section [](#sec:sys). Then, a call to `MFNSolve` runs the solver on a given vector $v$, returning the computed result $y$. Finally, `MFNDestroy` is used to reclaim memory. We give a few more details below.

### Defining the Problem

Defining the problem consists in specifying the matrix, $A$, and the function to compute, $f(\cdot)$. The problem matrix is provided with `MFNSetOperator`

```{code} c
MFNSetOperator(MFN mfn,Mat A);
```

where `A` should be a square matrix, stored in any allowed PETSc format including the matrix-free mechanism (see section [](#sec:supported)). The function $f(\cdot)$ is defined with an `FN` object. One possibility is to extract the `FN` object handled internally by `MFN`: `MFNGetFN`

```{code} c
MFNGetFN(MFN mfn,FN *f);
```

An alternative would be to create a standalone `FN` object and pass it with `MFNSetFN`. In any case, the function is defined via its type and the relevant parameters, see section [](#sec:sys) for details. The scaling parameters can be used for instance for the exponential when used in the context of ODE integration, $y=e^{tA}v$, where $t$ represents the elapsed time. Note that some `MFN` solvers may be restricted to only some types of `FN` functions.

In `MFN` it makes no sense to specify the number of eigenvalues. However, there is a related operation that allows the user to specify the size of the subspace that will be used internally by the solver (`ncv`, the number of column vectors of the basis): `MFNSetDimensions`

```{code} c
MFNSetDimensions(EPS eps,PetscInt ncv);
```

This parameter can also be set at run time with the option `-mfn_ncv`.

### Selecting the Solver

:::{table} List of solvers available in the `MFN` module.
:name: tab:mfnsolvers

 |Method                   |`MFNType`     |Options Database Name  |Supported Functions
 |-------------------------|--------------|-----------------------|-----------------------
 |Restarted Krylov solver  |`MFNKRYLOV`   |        `krylov`       |          Any
 |Expokit algorithm        |`MFNEXPOKIT`  |        `expokit`      |          Exponential

:::

The methods available in `MFN` are shown in table [](#tab:mfnsolvers). The solution method can be specified procedurally with `MFNSetType`

```{code} c
MFNSetType(MFN mfn,MFNType method);
```

or via the options database command `-mfn_type` followed by the method name (see table [](#tab:mfnsolvers)).

Currently implemented methods are:

-   A Krylov method with restarts as proposed by {cite:t}`Eiermann:2006:RKS`.

-   The method implemented in EXPOKIt {cite:p}`Sidje:1998:ESP` for the matrix exponential.

### Accuracy and Monitors

In the $f(A)v$ problem, there is no clear definition of residual, as opposed to the case of linear systems or eigenproblems. Still, the solvers have different ways of assessing the accuracy of the computed solution. The user can provide a tolerance and maximum number of iterations with `MFNSetTolerances`, but there is no guarantee that an analog of the residual is below the tolerance.

After the solver has finished, the number of performed (outer) iterations can be obtained with `MFNGetIterationNumber`. There are also monitors that display the error estimate, which can be activated with command-line keys `-mfn_monitor`, or `-mfn_monitor draw::draw_lg`. See section [](#sec:monitor) for additional details.

```{rubric} Footnotes
```

```{eval-rst}
.. bibliography::
   :filter: docname in docnames
```
