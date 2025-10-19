(ch:aux)=
# Auxiliary Classes

Apart from the main solver classes listed in table [](#tab:modules), SLEPc contains several auxiliary classes:

-   `ST`: Spectral Transformation, fully described in chapter [](#ch:st).

-   `FN`: Mathematical Function, required in application code to represent the constituent functions of the nonlinear operator in split form (chapter [](#ch:nep)), as well as the function to be used when computing the action of a matrix function on a vector (chapter [](#ch:mfn)).

-   `RG`: Region, a way to define a region of the complex plane.

-   `DS`: Direct Solver (or Dense System), can be seen as a wrapper to LAPACK functions used within SLEPc. It is mostly an internal object that need not be called by end users.

-   `BV`: Basis Vectors, provides the concept of a block of vectors that represent the basis of a subspace.

The first three classes, `ST`, `FN` and `RG`, are relevant for end users, while `DS` and `BV` are intended mainly for SLEPc developers. Below we provide a brief description of `FN`, `RG`, `DS` and `BV`.

{#sec:fn}
## FN: Mathematical Functions

The `FN` class provides a few predefined mathematical functions, including rational functions (of which polynomials are a particular case) and exponentials. Objects of this class are instantiated by providing the values of the relevant parameters. `FN` objects are created with `FNCreate` and it is necessary to select the type of function (rational, exponential, etc.) with `FNSetType`. Table [](#tab:fn) lists available functions.

:::{table} Mathematical functions available as `FN` objects
:name: tab:fn

 |Function                 |`FNType`      |Expression
 |-------------------------|--------------|-------------------------------------
 |Polynomial and rational  |`FNRATIONAL`  |$p(x)/q(x)$
 |Exponential              |`FNEXP`       |$e^x$
 |Logarithm                |`FNLOG`       |$\log x$
 |$\varphi$-functions      |`FNPHI`       |$\varphi_0(x)$, $\varphi_1(x)$, ...
 |Square root              |`FNSQRT`      |$\sqrt{x}$
 |Inverse square root      |`FNINVSQRT`   |$x^{-\frac{1}{2}}$
 |Combine two functions    |`FNCOMBINE`   |See text

:::

Parameters common to all `FN` types are the scaling factors, which are set with:

```{code} c
FNSetScale(FN fn,PetscScalar alpha,PetscScalar beta);
```

where `alpha` multiplies the argument and `beta` multiplies the result. With this, the actual function is $\beta\cdot f(\alpha\cdot x)$ for a given function $f(\cdot)$. For instance, an exponential function $f(x)=e^x$ will turn into

```{math}
:label: eq:exp-function

g(x)=\beta e^{\alpha x}.
```

In a rational function there are specific parameters, namely the coefficients of the numerator and denominator,

```{math}
:label: eq:num-denum-coefficients

r(x)=\frac{p(x)}{q(x)}
=\frac{\nu_{n-1}x^{n-1}+\cdots+\nu_1x+\nu_0}{\delta_{m-1}x^{m-1}+\cdots+\delta_1x+\delta_0}.
```

These parameters are specified with:

```{code} c
FNRationalSetNumerator(FN fn,PetscInt np,PetscScalar *pcoeff);
FNRationalSetDenominator(FN fn,PetscInt nq,PetscScalar *qcoeff);
```

Here, polynomials are passed as an array with high order coefficients appearing in low indices.

The $\varphi$-functions are given by

```{math}
:label: eq:phi-functions

\varphi_0(x)=e^x,\qquad \varphi_1(x)=\frac{e^x-1}{x},\qquad \varphi_k(x)=\frac{\varphi_{k-1}(x)-1/(k-1)!}{x},
```

where the index $k$ must be specified with `FNPhiSetIndex`.

Whenever the solvers need to compute $f(x)$ or $f'(x)$ on a given scalar $x$, the following functions are invoked:

```{code} c
FNEvaluateFunction(FN fn,PetscScalar x,PetscScalar *y)
FNEvaluateDerivative(FN fn,PetscScalar x,PetscScalar *y)
```

The function can also be evaluated as a matrix function, $B=f(A)$, where $A,B$ are small, dense, square matrices. This is done with `FNEvaluateFunctionMat`. Note that for a rational function, the corresponding expression would be $q(A)^{-1}p(A)$. For computing functions such as the exponential of a small matrix $A$, several methods are available. When the matrix $A$ is symmetric, the default is to compute $f(A)$ using the eigendecomposition $A=Q\Lambda Q^*$, for instance the exponential would be computed as $\exp(A)=Q\,\mathrm{diag}(e^{\lambda_i})Q^*$. In the general case, it is necessary to have recourse to one of the methods discussed in, e.g., {cite:p}`Hig10`.

Finally, there is a mechanism to combine simple functions in order to create more complicated functions. For instance, the function

```{math}
:label: eq:combined-funct

f(x) = (1-x^2) \exp\left( \frac{-x}{1+x^2} \right)
```

can be represented with an expression tree with three leaves (one exponential function and two rational functions) and two interior nodes (one of them is the root, $f(x)$). Interior nodes are simply `FN` objects of type `FNCOMBINE` that specify how the two children must be combined (with either addition, multiplication, division or function composition):

```{code} c
FNCombineSetChildren(FN fn,FNCombineType comb,FN f1,FN f2)
```

The combination of $f_1$ and $f_2$ with division will result in $f_1(x)/f_2(x)$ and $f_2(A)^{-1}f_1(A)$ in the case of matrices.

{#sec:rg}
## RG: Region

The `RG` object defines a region of the complex plane, that can be used to specify where eigenvalues must be sought. Currently, the following types of regions are available:

-   A (generalized) interval, defined as $[a,b]\times[c,d]$, where the four parameters can be set with `RGIntervalSetEndpoints`. This covers the particular cases of an interval on the real axis (setting $c=d=0$), the left halfplane $[-\infty,0]\times[-\infty,+\infty]$, a quadrant, etc. See figure [](#fig:rg-interval).
```{figure} ../../_static/images/manual/svg/fig-rg-interval.svg
:alt: Interval region defined via de RG class
:name: fig:rg-interval

Interval region defined via de RG class
```

-   A polygon defined by its vertices, given via `RGPolygonSetVertices`. See figure [](#fig:rg-polygon).
```{figure} ../../_static/images/manual/svg/fig-rg-polygon.svg
:alt: Polygon region defined via de RG class
:name: fig:rg-polygon

Polygon region defined via de RG class
```

-   An ellipse defined by its center, radius and vertical scale (1 by default), specified with `RGEllipseSetParameters`. See figure [](#fig:rg-ellipse).
```{figure} ../../_static/images/manual/svg/fig-rg-ellipse.svg
:alt: Ellipse region defined via de RG class
:name: fig:rg-ellipse

Ellipse region defined via de RG class
```

-   A ring region similar to an ellipse but consisting of a thin stripe along the ellipse with optional start and end angles. See figure [](#fig:rg-ring). The parameters are set with `RGRingSetParameters`.
```{figure} ../../_static/images/manual/svg/fig-rg-ring.svg
:alt: Ring region defined via de RG class
:name: fig:rg-ring

Ring region defined via de RG class
```

Check table [](tab:rg) for the names that should be used in each case.

:::{table} Regions available as `RG` objects
:name: tab:rg

 |Region Type             |`RGType`      |Options Database
 |------------------------|--------------|------------------
 |(Generalized) Interval  |`RGINTERVAL`  |   `interval`
 |Polygon                 |`RGPOLYGON`   |   `polygon`
 |Ellipse                 |`RGELLIPSE`   |   `ellipse`
 |Ring                    |`RGRING`      |   `ring`

:::

Sometimes it is useful to specify the complement of a certain region, e.g., the part of the complex plane outside an ellipse. This can be achieved with:

```{code} c
RGSetComplement(RG rg,PetscBool flg)
```

or in the command line with `-rg_complement`.

By default, a newly created `RG` object that is not set a type nor parameters must represent the whole complex plane (the same as `RGINTERVAL` with values $[-\infty,+\infty]\times[-\infty,+\infty]$). We call this the *trivial* region, and provide a function to test this situation:

```{code} c
RGIsTrivial(RG rg,PetscBool *trivial)
```

Another useful operation is to check whether a given point of the complex plane is inside the region or not:

```{code} c
RGCheckInside(RG rg,PetscInt n,PetscScalar *ar,PetscScalar *ai,PetscInt *inside)
```

Note that the point is represented as two {external:doc}`PetscScalar`'s, similarly to eigenvalues in SLEPc.

{#sec:bv}
## BV: Basis Vectors

The `BV` class may be useful for advanced users, so we briefly describe it here for completeness. `BV` is a convenient way of handling a collection of vectors that often operate together, rather than working with an array of {external:doc}`Vec`. It can be seen as a generalization of {external:doc}`Vec` to a tall-skinny matrix with several columns.

:::{table} Operations available for `BV` objects
:name: tab:bv

 |Operation              |Block version      |Column version           |Vector version
 |-----------------------|-------------------|-------------------------|----------------------
 |$Y=X$                  |`BVCopy`           |`BVCopyColumn`           |`BVCopyVec`
 |$Y=\beta Y+\alpha XQ$  |`BVMult`           |`BVMultColumn`           |`BVMultVec`
 |$M=Y^*\!AX$            |`BVMatProject`     |--                       |--
 |$M=Y^*X$               |`BVDot`            |`BVDotColumn`            |`BVDotVec`
 |$Y=\alpha Y$           |`BVScale`          |`BVScaleColumn`          |--
 |$r=\\|X\\|_{type}$     |`BVNorm`           |`BVNormColumn`           |`BVNormVec`
 |Set to random values   |`BVSetRandom`      |`BVSetRandomColumn`      |--
 |Orthogonalize          |`BVOrthogonalize`  |`BVOrthogonalizeColumn`  |`BVOrthogonalizeVec`

:::

Table [](#tab:bv) shows a summary of the operations offered by the `BV` class, with variants that operate on the whole `BV`, on a single column, or on an external {external:doc}`Vec` object. Missing variants can be achieved simply with {external:doc}`Vec` and {external:doc}`Mat` operations. Other available variants not shown in the table are `BVMultInPlace`, `BVMultInPlaceHermitianTranspose` and `BVOrthogonalizeSomeColumn`.

Most SLEPc solvers use a `BV` object to represent the working subspace basis. In particular, orthogonalization operations are mostly confined within `BV`. Hence, `BV` provides options for specifying the method of orthogonalization of vectors (Gram-Schmidt) as well as the method of block orthogonalization, see `BVSetOrthogonalization`.

```{only} html
<p class="rubric">References</p>
```
```{bibliography}
:filter: docname in docnames
```
