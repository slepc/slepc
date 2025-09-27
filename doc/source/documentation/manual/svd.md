(ch:svd)=
# SVD: Singular Value Decomposition

The Singular Value Decomposition (`SVD`) solver object can be used for computing a partial SVD of a rectangular matrix, and other related problems. It provides uniform and efficient access to several specific SVD solvers included in SLEPc, and also gives the possibility to compute the decomposition via the eigensolvers provided in the `EPS` package.

In many aspects, the user interface of `SVD` resembles that of `EPS`. For this reason, this chapter and chapter [](#ch:eps) have a very similar structure.

{#sec:svdback label="sec:svdback"}
## Mathematical Background

This section provides some basic concepts about the singular value decomposition and other related problems. The objective is to set up the notation and also to justify some of the solution approaches, particularly those based on the `EPS` object. As in the case of eigensolvers, some of the implemented methods are described in detail in the SLEPc [technical reports](#str).

For background material about the SVD, see for instance {cite:p}`Bai:2000:TSA{ch. 6}`. Many other books such as {cite:p}`Bjorck:1996:NML` or {cite:p}`Hansen:1998:RDI` present the SVD from the perspective of its application to the solution of least squares problems and other related linear algebra problems.

{#sec:svd label="sec:svd"}
### The (Standard) Singular Value Decomposition (SVD)

The singular value decomposition (SVD) of an $m\times n$ matrix $A$ can be written as

```{math}
:label: eq:svd

A=U\Sigma V^*,
```

 where $U=[u_1,\ldots,u_m]$ is an $m\times m$ unitary matrix ($U^*U=I$), $V=[v_1,\ldots,v_n]$ is an $n\times n$ unitary matrix ($V^*V=I$), and $\Sigma$ is an $m\times n$ diagonal matrix with real diagonal entries $\Sigma_{ii}=\sigma_i$ for $i=1,\ldots,\min\{m,n\}$. If $A$ is real, $U$ and $V$ are real and orthogonal. The vectors $u_i$ are called the left singular vectors, the $v_i$ are the right singular vectors, and the $\sigma_i$ are the singular values.

```{figure} ../../_static/images/manual/svg/fig-svd.svg
:alt: Scheme of the thin SVD of a rectangular matrix $A$.
:name: fig:svd

Scheme of the thin SVD of a rectangular matrix $A$.
```

In the following, we will assume that $m\geq n$. If $m<n$ then $A$ should be replaced by $A^*$ (note that in SLEPc this is done transparently as described later in this chapter and the user need not worry about this). In the case that $m\geq n$, the top $n$ rows of $\Sigma$ contain $\mathrm{diag}(\sigma_1,\ldots,\sigma_n)$ and its bottom $m-n$ rows are zero. The relation equation {math:numref}`eq:svd` may also be written as $AV=U\Sigma$, or

```{math}
:label: eq:svdleft

Av_i=u_i\sigma_i\;,\quad i=1,\ldots,n,
```

 and also as $A^*U=V\Sigma^*$, or

```{math}
:label: eq:svdright

\begin{aligned}
A^*u_i&=v_i\sigma_i\;,\quad i=1,\ldots,n,\\
\end{aligned}
```

```{math}
:label: eq:svdright2

\begin{aligned}
A^*u_i&=0\;,\quad i=n+1,\ldots,m.
\end{aligned}
```

 The last left singular vectors corresponding to equation {math:numref}`eq:svdright2` are often not computed, especially if $m\gg n$. In that case, the resulting factorization is sometimes called the *thin* SVD, $A=U_n\Sigma_n V_n^*$, and is depicted in figure [](#fig:svd). This factorization can also be written as

```{math}
:label: eq:svdouter

A=\sum_{i=1}^{n}\sigma_iu_iv_i^*.
```

 Each $(\sigma_i,u_i,v_i)$ is called a singular triplet.

The singular values are real and nonnegative, $\sigma_1\geq\sigma_2\geq\ldots\geq\sigma_r>\sigma_{r+1}=\ldots=\sigma_n=0$, where $r=\mathrm{rank}(A)$. It can be shown that $\{u_1,\ldots,u_r\}$ span the range of $A$, $\mathcal{R}(A)$, whereas $\{v_{r+1},\ldots,v_n\}$ span the null space of $A$, $\mathcal{N}(A)$.

If the zero singular values are dropped from the sum in equation {math:numref}`eq:svdouter`, the resulting factorization, $A=\sum_{i=1}^{r}\sigma_iu_iv_i^*$, is called the *compact* SVD, $A=U_r\Sigma_r V_r^*$.

In the case of a very large and sparse $A$, it is usual to compute only a subset of $k\leq r$ singular triplets. We will refer to this decomposition as the *truncated* SVD of $A$. It can be shown that the matrix $A_k=U_k\Sigma_k V_k^*$ is the best rank-$k$ approximation to matrix $A$, in the least squares sense.

In general, one can take an arbitrary subset of the summands in equation {math:numref}`eq:svdouter`, and the resulting factorization is called the *partial* SVD of $A$. As described later in this chapter, SLEPc allows the computation of a partial SVD corresponding to either the $k$ largest or smallest singular triplets.

#### Equivalent Eigenvalue Problems

It is possible to formulate the problem of computing the singular triplets of a matrix $A$ as an eigenvalue problem involving a Hermitian matrix related to $A$. There are two possible ways of achieving this:

1.  With the *cross product* matrix, either $A^*A$ or $AA^*$.

2.  With the *cyclic* matrix, $H(A)=\left[\begin{smallmatrix}0&A\\A^*&0\end{smallmatrix}\right]$.

In SLEPc, the computation of the SVD is usually based on one of these two alternatives, either by passing one of these matrices to an `EPS` object or by performing the computation implicitly.

By pre-multiplying equation {math:numref}`eq:svdleft` by $A^*$ and then using equation {math:numref}`eq:svdright`, the following relation results

```{math}
:label: eq:eigleft

A^*Av_i=\sigma_i^2v_i,
```

 that is, the $v_i$ are the eigenvectors of matrix $A^*A$ with corresponding eigenvalues equal to $\sigma_i^2$. Note that after computing $v_i$ the corresponding left singular vector, $u_i$, is readily available through equation {math:numref}`eq:svdleft` with just a matrix-vector product, $u_i=\frac{1}{\sigma_i}Av_i$.

Alternatively, one could first compute the left vectors and then the right ones. For this, pre-multiply equation {math:numref}`eq:svdright` by $A$ and then use equation {math:numref}`eq:svdleft` to get

```{math}
:label: eq:eigright

AA^*u_i=\sigma_i^2u_i.
```

 In this case, the right singular vectors are obtained as $v_i=\frac{1}{\sigma_i}A^*u_i$.

The two approaches represented in equations {math:numref}`eq:eigleft` and {math:numref}`eq:eigright` are very similar. Note however that $A^*A$ is a square matrix of order $n$ whereas $AA^*$ is of order $m$. In cases where $m\gg n$, the computational effort will favor the $A^*A$ approach. On the other hand, the eigenproblem equation {math:numref}`eq:eigleft` has $n-r$ zero eigenvalues and the eigenproblem equation {math:numref}`eq:eigright` has $m-r$ zero eigenvalues. Therefore, continuing with the assumption that $m\geq n$, even in the full rank case the $AA^*$ approach may have a large null space resulting in difficulties if the smallest singular values are sought. In SLEPc, this will be referred to as the cross product approach and will use whichever matrix is smaller, either $A^*A$ or $AA^*$.

Computing the SVD via the cross product approach may be adequate for determining the largest singular triplets of $A$, but the loss of accuracy can be severe for the smallest singular triplets. The cyclic matrix approach is an alternative that avoids this problem, but at the expense of significantly increasing the cost of the computation. Consider the eigendecomposition of

```{math}
:label: eq:cyclic

H(A)=\begin{bmatrix}0&A\\A^*&0\end{bmatrix},
```

 which is a Hermitian matrix of order $(m+n)$. It can be shown that $\pm\sigma_i$ is a pair of eigenvalues of $H(A)$ for $i=1,\ldots,r$ and the other $m+n-2r$ eigenvalues are zero. The unit eigenvectors associated with $\pm\sigma_i$ are $\frac{1}{\sqrt{2}}\left[\begin{smallmatrix}\pm u_i\\v_i\end{smallmatrix}\right]$. Thus it is possible to extract the singular values and the left and right singular vectors of $A$ directly from the eigenvalues and eigenvectors of $H(A)$. Note that in this case the singular values are not squared, and therefore the computed values will be more accurate (especially the small ones). The drawback in this case is that small eigenvalues are located in the interior of the spectrum.

{#sec:gsvd label="sec:gsvd"}
### The Generalized Singular Value Decomposition (GSVD)

An extension of the SVD to the case of two matrices is the generalized singular value decomposition (GSVD), which can be applied in constrained least squares problems, among others. An overview of the problem can be found in {cite:p}`Golub:1996:MC{8.7.3}`.

Consider two matrices, $A\in\mathbb{C}^{m\times n}$ with $m\geq n$ and $B\in\mathbb{C}^{p\times n}$. Note that both matrices must have the same column dimension. Then there exist two unitary matrices $U\in\mathbb{C}^{m\times m}$ and $V\in\mathbb{C}^{p\times p}$ and an invertible matrix $X\in\mathbb{C}^{n\times n}$ such that

```{math}
:label: eq:gsvd

U^*AX=C,\qquad V^*BX=S,
```

 where $C=\mathrm{diag}(c_1,\dots,c_n)$ and $S=\mathrm{diag}(s_{n-q+1},\dots,s_n)$ with $q=\min(p,n)$. The values $c_i$ and $s_i$ are real and nonnegative, and the ratios define the generalized singular values,

```{math}
:label: eq:gsvd-values

\sigma(A,B)\equiv\{c_1/s_1,\dots,c_q/s_q\},
```

and if $p<n$ we can consider that the first $n-p$ generalized singular values are infinite, as if $s_1=\dots=s_{n-p}=0$. Note that if $B$ is the identity matrix, $X$ can be taken to be unitary and then we recover the standard SVD, $\sigma(A,I)=\sigma(A)$, that is why equation {math:numref}`eq:gsvd` is considered a generalization of the SVD.

The diagonal matrices $C$ and $S$ satisfy $C^*C+S^*S=I$, and are related to the CS decomposition {cite:p}`Golub:1996:MC{2.6.4}` associated with the orthogonal factor of the QR factorization of matrices $A$ and $B$ stacked, that is, if

```{math}
:label: eq:qr

Z:=\begin{bmatrix}A\\B\end{bmatrix}=\begin{bmatrix}Q_1\\Q_2\end{bmatrix}R,
```

 then $C$ and $S$ can be obtained from the singular values of $Q_1$ and $Q_2$, respectively. The matrix $Z$ is relevant for algorithms and is often built explicitly.

```{figure} ../../_static/images/manual/svg/fig-gsvd.svg
:alt: Scheme of the thin GSVD of two matrices $A$ and $B$, for the case $p < n$.
:name: fig:gsvd

Scheme of the thin GSVD of two matrices $A$ and $B$, for the case $p < n$.
```

The above description assumes that matrix $Z$ has full column rank, and the rank of $B$ is also $q$. In the general case where these assumptions do not hold, the structure of matrices $C$ and $S$ is a bit more complicated. This includes also the case where $m<n$. A detailed description of those cases can be found in {cite:p}`Anderson:1999:LUG{2.3.5.3}`.

As in the case of the SVD, one can consider thin, compact, truncated, and partial versions of the GSVD. A pictorial representation of the thin GSVD is shown in figure [](#fig:gsvd).

The columns of $X$, $x_i$, are called the (right) generalized singular vectors. The left vectors in this case would correspond to the columns of $U$ and $V$. In SLEPc, the user interface will provide these vectors stacked on top of each other, as a single $(m+p)$-vector $\begin{bmatrix}u_i\\v_i\end{bmatrix}$.

#### Equivalent Eigenvalue Problems

In the GSVD it is also possible to formulate the problem as an eigenvalue problem, which opens the door to approach its solution via `EPS`. The columns of $X$ satisfy

```{math}
:label: eq:gsvdeigcross

s_i^2A^*Ax_i=c_i^2B^*Bx_i,
```

and so if $s_i\neq 0$ then $A^*Ax_i=\sigma_i^2B^*Bx_i$, a generalized eigenvalue problem for the matrix pair $(A^*A,B^*B)$. This is the analog of the cross product matrix eigenproblem of equation {math:numref}`eq:eigleft`.

The formulation that is analog to the eigenproblem associated with the cyclic matrix equation {math:numref}`eq:cyclic` is to solve the generalized eigenvalue problem defined by any of the matrix pairs

```{math}
:label: eq:gsvdeigcyclic

\left(
\begin{bmatrix}0&A\\A^*&0\end{bmatrix},
\begin{bmatrix}I&0\\0&B^*B\end{bmatrix}
\right),
\qquad\text{or}\qquad
\left(
\begin{bmatrix}0&B\\B^*&0\end{bmatrix},
\begin{bmatrix}I&0\\0&A^*A\end{bmatrix}
\right).
```

{#sec:hsvd label="sec:hsvd"}
### The Hyperbolic Singular Value Decomposition (HSVD)

The hyperbolic singular value decomposition (HSVD) was introduced in {cite:p}`Onn:1991:HSV`, motivated by some signal processing applications such as the so-called covariance differencing problem. The formulation of the HSVD is similar to that of the SVD, except that $U$ is orthogonal with respect to a signature matrix,

```{math}
:label: eq:hsvd

A=U\Sigma V^*,\qquad U^*\Omega U=\tilde\Omega,
```

 where $\Omega=\mathrm{diag}(\pm 1)$ is an $m\times m$ signature matrix provided by the user, while $\tilde\Omega$ is another signature matrix obtained as part of the solution. Sometimes $U$ is said to be a hyperexchange matrix, or also an $(\Omega,\tilde\Omega)$-orthogonal matrix. Note that in the problem definition normally found in the literature it is $V$ that is $(\Omega,\tilde\Omega)$-orthogonal and not $U$. We choose this definition for consistency with respect to the generalized HSVD of two matrices. If the user wants to compute the HSVD according to the alternative definition, then it suffices to (conjugate) transpose the input matrix $A$, using for instance {external:doc}`MatHermitianTranspose`.

As in the case of the SVD, the solution of the problem consists in singular triplets $(\sigma_i,u_i,v_i)$, with $\sigma_i$ real and nonnegative and sorted in nonincreasing order. Note that these quantities are different from those of section [](#sec:svd), even though we use the same notation here. With each singular triplet, there is an associated sign $\tilde\omega_i$ (either 1 or $-1$), the corresponding diagonal element of $\tilde\Omega$. In SLEPc, this value is not returned by the user interface, but if required it can be easily computed as $u_i^*\Omega u_i$.

The relations between left and right singular vectors are slightly different from those of the standard SVD. We have $AV=U\Sigma$ and $A^*\Omega U=V\Sigma^*\tilde\Omega$, so for $m\geq n$ the following relations hold, together with equation {math:numref}`eq:svdleft`:

```{math}
:label: eq:hsvdright

\begin{aligned}
A^*\Omega u_i&=v_i\sigma_i\tilde\omega_i\;,\quad i=1,\ldots,n,\\
\end{aligned}
```

```{math}
:label: eq:hsvdright2

\begin{aligned}
A^*\Omega u_i&=0\;,\quad i=n+1,\ldots,m.
\end{aligned}
```

 In SLEPc we will compute a partial HSVD consisting of either the largest or smallest hyperbolic singular triplets. Note that the sign $\tilde\omega_i$ is not used when sorting for largest or smallest $\sigma_i$.

#### Equivalent Eigenvalue Problems

Once again, we can derive cross and cyclic schemes to compute the decomposition by solving an eigenvalue problem. The cross product matrix approach has two forms, to be selected depending on whether $m\geq n$ or not, as discussed in section [](#sec:svd). The first form,

```{math}
:label: eq:heigleft

A^*\Omega Av_i=\sigma_i^2\tilde\omega_iv_i,
```

is derived by pre-multiplying equation {math:numref}`eq:svdleft` by $A^*\Omega$ and then using equation {math:numref}`eq:hsvdright`. This eigenproblem can be solved as a HEP (cf. section [](#sec:defprob)) and may have both positive and negative eigenvalues, corresponding to $\tilde\omega_i=1$ and $\tilde\omega_i=-1$, respectively. Once the right vector $v_i$ has been computed, the corresponding left vector can be obtained using equation {math:numref}`eq:svdleft` with just a matrix-vector product, $u_i=\sigma_i^{-1}Av_i$.

The second form of cross computes the left vectors first, by pre-multiplying equation {math:numref}`eq:hsvdright` by $A$ and then using equation {math:numref}`eq:svdleft`,

```{math}
:label: eq:heigright

AA^*\Omega u_i=\sigma_i^2\tilde\omega_iu_i.
```

 In this case, the right singular vectors are obtained as $v_i=(\sigma_i\tilde\omega_i)^{-1}A^*\Omega u_i$. The coefficient matrix of {math:numref}`eq:heigright` is non-Hermitian, so the eigenproblem has to be solved as non-Hermitian, or alternatively it can be formulated as a generalized eigenvalue problem of GHIEP type (cf. section [](#sec:defprob)) for the indefinite pencil $(AA^*,\Omega)$,

```{math}
:label: eq:heigright2

AA^*\hat{u}_i=\sigma_i^2\tilde\omega_i\Omega \hat{u}_i,
```

 with $\hat{u}_i=\Omega u_i$. The eigenvectors obtained from equation {math:numref}`eq:heigright` or equation {math:numref}`eq:heigright2` must be normalized so that $U^*\Omega U=\tilde\Omega$ holds.

Finally, in the cyclic matrix approach for the HSVD we must solve a generalized eigenvalue problem defined by the matrices of order $(m+n)$

```{math}
:label: eq:hcyclic

H(A)=\begin{bmatrix}0&A\\A^*&0\end{bmatrix},\qquad
\hat\Omega=\begin{bmatrix}\Omega&0\\0&I\end{bmatrix}.
```

As in the case of equation {math:numref}`eq:heigright2`, this pencil is Hermitian-indefinite and hence it may have complex eigenvalues. However, it can be shown that nonzero eigenvalues of the pencil $(H(A),\hat\Omega)$ are either real (equal to $\pm\sigma_i$ for a certain hyperbolic singular value $\sigma_i$) or purely imaginary (equal to $\pm\sigma_ij$ for a certain $\sigma_i$ with $j=\sqrt{-1}$). The associated eigenvectors are $\left[\begin{smallmatrix}\varsigma_i u_i\\v_i\end{smallmatrix}\right]$, where $\varsigma_i$ is either $\pm 1$ or $\pm j$.

## Basic Usage

From the perspective of the user interface, the `SVD` package is very similar to `EPS`, with some differences that will be highlighted shortly.

```{code-block} c
:name: fig:ex-svd
:caption: Example code for basic solution with `SVD`.

SVD       svd;       /*  SVD solver context  */
Mat       A;         /*  problem matrix      */
Vec       u, v;      /*  singular vectors    */
PetscReal sigma;     /*  singular value      */
PetscInt  j, nconv;
PetscReal error;

SVDCreate( PETSC_COMM_WORLD, &svd );
SVDSetOperators( svd, A, NULL );
SVDSetProblemType( svd, SVD_STANDARD );
SVDSetFromOptions( svd );
SVDSolve( svd );
SVDGetConverged( svd, &nconv );
for (j=0; j<nconv; j++) {
  SVDGetSingularTriplet( svd, j, &sigma, u, v );
  SVDComputeError( svd, j, SVD_ERROR_RELATIVE, &error );
}
SVDDestroy( &svd );
```

The basic steps for computing a partial SVD with SLEPc are illustrated in figure [](#fig:ex-svd). The steps are more or less the same as those described in chapter [](#ch:eps) for the eigenvalue problem. First, the solver context is created with `SVDCreate`. Then the problem matrices have to be specified with `SVDSetOperators` and the type of problem can be selected via `SVDSetProblemType`. Then, a call to `SVDSolve` invokes the actual solver. After that, `SVDGetConverged` is used to determine how many solutions have been computed, which are retrieved with `SVDGetSingularTriplet`. Finally, `SVDDestroy` gets rid of the object.

If one compares this example code with the `EPS` example in figure [](#fig:ex-eps), the most outstanding differences are the following:

-   The singular value is a {external:doc}`PetscReal`, not a {external:doc}`PetscScalar`.

-   Each singular vector is defined with a single {external:doc}`Vec` object, not two as was the case for eigenvectors.

## Defining the Problem

Defining the problem consists in specifying the problem matrix $A$ (and also a second matrix $B$ in case of the GSVD), the problem type, and the portion of the spectrum to be computed.

The problem matrices are provided with the following function `SVDSetOperators`

```{code} c
SVDSetOperators(SVD svd,Mat A,Mat B);
```

where `A` can be any matrix, not necessarily square, stored in any allowed PETSc format including the matrix-free mechanism (see section [](#sec:supported) for a detailed discussion), and `B` is generally set to `NULL` except in the case of computing the GSVD. If the problem is hyperbolic (section [](#sec:hsvd)) then in addition a signature matrix $\Omega$ must be provided with `SVDSetSignature`

```{code} c
SVDSetSignature(SVD svd,Vec omega);
```

Note that $\Omega$ is represented as a vector containing values equal to 1 or $-1$.

The problem type can be specified with the function `SVDSetProblemType`

```{code} c
SVDSetProblemType(SVD svd,SVDProblemType type);
```

Note that in `SVD` calling it is currently not strictly required, since the problem type can be deduced from the information passed by the user: if only one matrix is passed to `SVDSetOperators` the problem type will default to a standard SVD (or hyperbolic SVD if a signature has been provided), while if two matrices are passed then it defaults to GSVD. Still, it is recommended to always call `SVDSetProblemType`. Table [](#tab:ptypesvd) lists the supported problem types.

:::{table} Problem types considered in `SVD`.
:name: tab:ptypesvd

 |Problem Type            |`SVDProblemType`   |Command line key
 |------------------------|-------------------|--------------------
 |Standard SVD            |`SVD_STANDARD`     |`-svd_standard`
 |Generalized SVD (GSVD)  |`SVD_GENERALIZED`  |`-svd_generalized`
 |Hyperbolic SVD (HSVD)   |`SVD_HYPERBOLIC`   |`-svd_hyperbolic`

:::

It is important to note that all SVD solvers in SLEPc make use of both $A$ and $A^*$, as suggested by the description in section [](#sec:svd). $A^*$ is not explicitly passed as an argument to `SVDSetOperators`, therefore it will have to stem from $A$. There are two possibilities for this: either $A$ is transposed explicitly and $A^*$ is created as a distinct matrix, or $A^*$ is handled implicitly via {external:doc}`MatMultTranspose` (or {external:doc}`MatMultHermitianTranspose` in the complex case) operations whenever a matrix-vector product is required in the algorithm. The default is to build $A^*$ explicitly, but this behavior can be changed with `SVDSetImplicitTranspose`

```{code} c
SVDSetImplicitTranspose(SVD svd,PetscBool impl);
```

In section [](#sec:svd), it was mentioned that in SLEPc the cross product approach chooses the smallest of the two possible cases $A^*A$ or $AA^*$, that is, $A^*A$ is used if $A$ is a tall, thin matrix ($m\geq n$), and $AA^*$ is used if $A$ is a fat, short matrix ($m<n$). In fact, what SLEPc does internally is that if $m<n$ the roles of $A$ and $A^*$ are reversed. This is equivalent to transposing all the SVD factorization, so left singular vectors become right singular vectors and vice versa. This is actually done in all singular value solvers, not only the cross product approach. The objective is to simplify the number of cases to be treated internally by SLEPc, as well as to reduce the computational cost in some situations. Note that this is done transparently and the user need not worry about transposing the matrix, only to indicate how the transpose has to be handled, as explained above.

The user can specify how many singular values and vectors to compute. The default is to compute only one singular triplet. The function `SVDSetDimensions`

```{code} c
SVDSetDimensions(SVD svd,PetscInt nsv,PetscInt ncv,PetscInt mpd);
```

allows the specification of the number of singular values to compute, `nsv`. The next argument can be set to prescribe the number of column vectors to be used by the solution algorithm, `ncv`, that is, the largest dimension of the working subspace. These two parameters can also be set at run time with the options `-svd_nsv` and `-svd_ncv`. For example, the command line

```{code} console
$ ./program -svd_nsv 10 -svd_ncv 24
```

requests 10 singular values and instructs to use 24 column vectors. Note that `ncv` must be at least equal to `nsv`, although in general it is recommended (depending on the method) to work with a larger subspace, for instance $\mathtt{ncv}\geq2\cdot\mathtt{nsv}$ or even more. As in the case of the `EPS` object, the last argument, `mpd`, can be used to limit the maximum dimension of the projected problem, as discussed in section [](#sec:large-nev). Using this parameter is especially important in the case that a large number of singular values are requested.

:::{table} Available possibilities for selection of the singular values of interest.
:name: tab:whichsvd

 |`SVDWhich`      |Command line key  |Sorting criterion
 |----------------|------------------|-------------------
 |`SVD_LARGEST`   |`-svd_largest`    |Largest $\sigma$
 |`SVD_SMALLEST`  |`-svd_smallest`   |Smallest $\sigma$

:::

For the selection of the portion of the spectrum of interest, there are only two possibilities in the case of SVD: largest and smallest singular values, see table [](#tab:whichsvd). The default is to compute the largest ones, but this can be changed with `SVDSetWhichSingularTriplets`

```{code} c
SVDSetWhichSingularTriplets(SVD svd,SVDWhich which);
```

which can also be specified at the command line. This criterion is used both for configuring how the algorithm seeks singular values and also for sorting the computed values. In contrast to the case of `EPS`, computing singular values located in the interior part of the spectrum is difficult, the only possibility is to use an `EPS` object combined with a spectral transformation (this possibility is explained in detail in the next section). Note that in this case, the value of `which` applies to the transformed spectrum.

{#sec:thres label="sec:thres"}
### Using a threshold to specify wanted singular values

In some applications, the number of wanted singular values is not known a priori. For instance, suppose that matrix $A$ is known to have low rank, and we want to approximate it by a truncated SVD containing the leading singular values. The goal is to have a good approximation, that is, discard only small singular values. The numerical rank $k$ might be difficult to estimate, due to discretization error or other reasons. But since we assume that singular values decay abruptly around some unknown value $k$, we can configure SLEPc to detect this decay.

The threshold $\delta$ can be set with `SVDSetThreshold`

```{code} c
SVDSetThreshold(SVD svd,PetscReal thres,PetscBool rel);
```

When `rel` is {external:doc}`PETSC_TRUE`, the solver will compute $k$ singular values, with $\sigma_j\geq\delta\cdot\sigma_1$, for $j=1,\dots,k-1$, where $k$ is computed internally. In this case, the threshold is relative to the largest singular value $\sigma_1$, i.e., the matrix norm. It can be interpreted as a percentage, for instance $\delta=0.8$ means compute singular values that are at least 80% of the matrix norm. Alternatively, an absolute threshold can be used by setting `rel` to {external:doc}`PETSC_FALSE`

```{figure} ../../_static/images/manual/svg/fig-thres.svg
:alt: Illustration of threshold usage with the singular values of the `rdb200` matrix. The red line represents a threshold of {math}`\delta=0.8`.
:name: fig:thres

Illustration of threshold usage with the singular values of the `rdb200` matrix. The red line represents a threshold of {math}`\delta=0.8`.
```

In the low-rank example suggested above, the singular values will decay fast for a relatively small $k$, so a relative threshold $\delta=0.5$ or less could do the job. But care must be taken in matrices whose singular values decay progressively, as in the example of figure [Illustration of threshold usage with the singular values of the `rdb200` matrix](#fig:thres), since a small $\delta$ would imply computing many singular triplets and hence a very high cost, both computationally and in memory. Since the number of computed singular values is not known a priori, the solver will need to reallocate the basis of vectors internally, to have enough room to accommodate all the singular vectors, so this option must be used with caution to avoid out-of-memory problems. The recommendation is to set the value of `ncv` to be larger than the estimated number of singular values, to minimize the number of reallocations. You can also use the `nsv` parameter in combination with the threshold to stop in case the number of computed singular triplets exceeds that value.

An absolute threshold is also available when computing smallest singular values, in which case the solver will compute the $k$ smallest singular values, where $\sigma_j<\delta$, $j=n-k+1,\dots,n$.

## Selecting the SVD Solver

The available methods for computing the partial SVD are shown in table [](#tab:svdsolvers). These methods can be classified in the following categories:

-   Solvers based on `EPS`. These solvers set up an `EPS` object internally, thus using the available eigensolvers for solving the SVD problem. The two possible approaches in this case are the cross product matrix and the cyclic matrix, as described in section [](#sec:svdback).

-   Specific SVD solvers. These are typically eigensolvers that have been adapted algorithmically to exploit the structure of the SVD or GSVD problems. There are two Lanczos-type solvers in this category: Lanczos and thick-restart Lanczos, see {cite:t}`str-8` for a detailed description of these methods. In this category, we could also add the randomized SVD (RSVD), a special solver that does not compute individual singular vectors accurately, but rather a low-rank approximation of $A$ by means of randomized techniques.

-   The LAPACK solver. This is an interface to some LAPACK routines, analog of those in the case of eigenproblems. These routines operate in dense mode with only one processor and therefore are suitable only for moderate size problems. This solver should be used only for debugging purposes.

-   External packages: SCALAPAck and ELEMENTAL are dense packages and compute the complete SVD, while PRIMME offers Davidson-type methods to compute only a few singular triplets.

The default solver is the one that uses the cross product matrix (`cross`), usually the fastest and most memory-efficient approach, at least for the standard SVD. See a more detailed explanation below.

:::{table} List of solvers available in the `SVD` module. In the column of supported problems, 'G' means GSVD and 'H' means HSVD (the standard SVD is supported by all solvers).
:name: tab:svdsolvers

 |Method                 |`SVDType`        |Options Database Name  |Supported  Problems   |Default
 |-----------------------|-----------------|-----------------------|----------------------|---------
 |Cross Product          |`SVDCROSS`       |        `cross`        |           G,H        |$\star$
 |Cyclic Matrix          |`SVDCYCLIC`      |        `cyclic`       |           G,H
 |Lanczos                |`SVDLANCZOS`     |        `lanczos`      |
 |Thick-restart Lanczos  |`SVDTRLANCZOS`   |        `trlanczos`    |           G,H
 |Randomized SVD (RSVD)  |`SVDRANDOMIZED`  |        `randomized`   |
 |LAPACK solver          |`SVDLAPACK`      |        `lapack`       |           G
 |Wrapper to SCALAPAck   |`SVDSCALAPACK`   |        `scalapack`    |
 |Wrapper to KSVD        |`SVDKSVD`        |        `ksvd`         |
 |Wrapper to ELEMENTAL   |`SVDELEMENTAL`   |        `elemental`    |
 |Wrapper to PRIMME      |`SVDPRIMME`      |        `primme`       |

:::

The solution method can be specified procedurally or via the command line. The application programmer can set it by means of the command `SVDSetType`

```{code} c
SVDSetType(SVD svd,SVDType method);
```

while the user writes the options database command `-svd_type` followed by the name of the method (see table [](#tab:svdsolvers)).

Note that not all solvers in the `SVD` module support non-standard problems (sections [](#sec:gsvd)-[](#sec:hsvd)). Table [](#tab:svdsolvers) indicates which solvers can be used for these.

The `EPS`-based solvers deserve some additional comments. These SVD solvers work by creating an `EPS` object internally and setting up an eigenproblem of type `EPS_HEP` (or `EPS_GHEP` in the case of the GSVD, or `EPS_GHIEP` in the case of the HSVD). These solvers implement the cross product matrix and the cyclic matrix approaches as described in section [](#sec:svdback). Therefore, the operator matrix associated with the `EPS` object will be $A^*A$ in the case of the `cross` solver and $H(A)$ in the case of the `cyclic` solver (with variations in the case of non-standard problems).

In the case of the `cross` solver, the matrix $A^*A$ is not built explicitly by default, since the resulting matrix may be much less sparse than the original matrix. By default, a shell matrix is created internally in the `SVD` object and passed to the `EPS` object. Still, the user may choose to force the computation of $A^*A$ explicitly, by means of PETSc's sparse matrix-matrix product subroutines. This is set with `SVDCrossSetExplicitMatrix`

```{code} c
SVDCrossSetExplicitMatrix(SVD svd,PetscBool explicit);
```

In the `cyclic` solver the user can also choose between handling $H(A)$ implicitly as a shell matrix (the default), or forming it explicitly, that is, storing its elements in a distinct matrix. The function for setting this option is `SVDCyclicSetExplicitMatrix`

```{code} c
SVDCyclicSetExplicitMatrix(SVD svd,PetscBool explicit);
```

The `EPS` object associated with the `cross` and `cyclic` SVD solvers is created with a set of reasonable default parameters. However, it may sometimes be necessary to change some of the `EPS` options such as the eigensolver. To allow application programmers to set any of the `EPS` options directly within the code, the following routines are provided to extract the `EPS` context from the `SVD` object, `SVDCrossGetEPS` `SVDCyclicGetEPS`

```{code} c
SVDCrossGetEPS(SVD svd,EPS *eps);
SVDCyclicGetEPS(SVD svd,EPS *eps);
```

A more convenient way of changing `EPS` options is through the command-line. This is achieved simply by prefixing the `EPS` options with `-svd_cross_` or `-svd_cyclic_` as in the following example:

```{code} console
$ ./program -svd_type cross -svd_cross_eps_type gd
```

At this point, one may consider changing also the options of the `ST` object associated with the `EPS` object in `cross` and `cyclic` SVD solvers, for example to compute singular values located at the interior of the spectrum via a shift-and-invert transformation. This is indeed possible, but some considerations have to be taken into account. When $A^*A$ or $H(A)$ are managed as shell matrices, then the potential of the spectral transformation is limited seriously, because some of the required operations will not be defined if working with implicit matrices (this is discussed briefly in sections [](#sec:supported) and [](#sec:explicit)). The following example illustrates the computation of interior singular values with the `cyclic` solver with explicit $H(A)$ matrix:

```{code} console
$ ./program -svd_type cyclic -svd_cyclic_explicitmatrix
            -svd_cyclic_st_type sinvert -svd_cyclic_eps_target 12.0
            -svd_cyclic_st_ksp_type preonly -svd_cyclic_st_pc_type lu
```

Similarly, in the case of GSVD the thick-restart Lanczos solver uses a {external:doc}`KSP` solver internally, that can be configured by accessing it with `SVDTRLanczosGetKSP`

```{code} c
SVDTRLanczosGetKSP(SVD svd,KSP *ksp);
```

or with the corresponding command-line options prefixed with `-svd_trlanczos_`. This {external:doc}`KSP` object is used to solve a linear least squares problem at each Lanczos step with coefficient matrix $Z$ equation {math:numref}`eq:qr`, which by default is a shell matrix but the user can choose to create it explicitly with the function `SVDTRLanczosSetExplicitMatrix`.

## Retrieving the Solution

Once the call to `SVDSolve` is complete, all the data associated with the computed partial SVD is kept internally in the `SVD` object. This information can be obtained by the calling program by means of a set of functions described below.

As in the case of eigenproblems, the number of computed singular triplets depends on the convergence and, therefore, it may be different from the number of solutions requested by the user. So the first task is to find out how many solutions are available, with `SVDGetConverged`

```{code} c
SVDGetConverged(SVD svd,PetscInt *nconv);
```

Usually, the number of converged solutions, `nconv`, will be equal to `nsv`, but in general it can be a number ranging from 0 to `ncv` (here, `nsv` and `ncv` are the arguments of function `SVDSetDimensions`).

Normally, the user is interested in the singular values only, or the complete singular triplets. The function `SVDGetSingularTriplet`

```{code} c
SVDGetSingularTriplet(SVD svd,PetscInt j,PetscReal *sigma,Vec u,Vec v);
```

returns the $j$-th computed singular triplet, $(\sigma_j,u_j,v_j)$, where both $u_j$ and $v_j$ are normalized to have unit norm[^hsvd]. Typically, this function is called inside a loop for each value of `j` from 0 to `nconv`--1. Note that singular values are ordered according to the same criterion specified with function `SVDSetWhichSingularTriplets` for selecting the portion of the spectrum of interest.

In some applications, it may be enough to compute only the right singular vectors. This is especially important in cases in which memory requirements are critical (remember that both $U_k$ and $V_k$ are dense matrices, and $U_k$ may require much more storage than $V_k$, see figure [](#fig:svd)). In SLEPc, there is no general option for specifying this, but the default behavior of some solvers is to compute only right vectors and allocate/compute left vectors only in the case that the user requests them. This is done in the `cross` solver and in some special variants of other solvers such as one-sided Lanczos (consult the {cite:t}`str-8` technical report for specific solver options).

In the case of the GSVD, the `sigma` argument of `SVDGetSingularTriplet` contains $\sigma_i=c_i/s_i$ and the second {external:doc}`Vec` argument (`v`) contains the right singular vectors ($x_i$), while the first {external:doc}`Vec` argument (`u`) contains the other vectors of the decomposition stacked on top of each other, as a single $(m+p)$-vector: $\left[\begin{smallmatrix}u_i\\v_i\end{smallmatrix}\right]$.

### Reliability of the Computed Solution

In SVD computations, a-posteriori error bounds are much the same as in the case of Hermitian eigenproblems, due to the equivalence discussed in section [](#sec:svd). The residual vector is defined in terms of the cyclic matrix, $H(A)$, so its norm is

```{math}
:label: eq:svd-residual-norm

\|r_\mathrm{SVD}\|_2=\left(\|A\tilde{v}-\tilde{\sigma}\tilde{u}\|_2^2+\|A^*\tilde{u}-\tilde{\sigma}\tilde{v}\|_2^2\right)^{\frac{1}{2}},
```

where $\tilde{\sigma}$, $\tilde{u}$ and $\tilde{v}$ represent any of the `nconv` computed singular triplets delivered by `SVDGetSingularTriplet`.

Given the above definition, the following relation holds

```{math}
:label: eq:svd-residual-norm-relation

|\sigma-\tilde{\sigma}|\leq \|r_\mathrm{SVD}\|_2,
```

where $\sigma$ is an exact singular value. The associated error can be obtained in terms of $\|r_\mathrm{SVD}\|_2$ with the following function: `SVDComputeError`

```{code} c
SVDComputeError(SVD svd,PetscInt j,SVDErrorType type,PetscReal *error);
```

In the case of the GSVD, the function `SVDComputeError` will compute a residual norm based on the two relations equation {math:numref}`eq:gsvd`,

```{math}
:label: eq:gsvd-residual-norm

\|r_\mathrm{GSVD}\|_2=\left(\|\tilde{s}^2A^*\tilde{u}-\tilde{c}B^*B\tilde{x}\|_2^2+\|\tilde{c}^2B^*\tilde{v}-\tilde{s}A^*A\tilde{x}\|_2^2\right)^{\frac{1}{2}},
```

 where $\tilde{x}$, $\tilde{u}$, $\tilde{v}$ are the computed singular vectors corresponding to $\tilde{\sigma}$, and $\tilde{c}$, $\tilde{s}$ are obtained from $\tilde{\sigma}$ as $\tilde{s}=1/\sqrt{1+\tilde{\sigma}^2}$ and $\tilde{c}=\tilde{\sigma}\tilde{s}$. See {cite:p}`Alvarruiz:2024:TLB` for details.

Similarly, in the HSVD we employ a modified residual

```{math}
:label: eq:hsvd-residual-norm

\|r_\mathrm{HSVD}\|_2=\left(\|A\tilde{v}-\tilde{\sigma}\tilde{u}\|_2^2+\|A^*\Omega\tilde{u}-\tilde{\sigma}\tilde{\omega}\tilde{v}\|_2^2\right)^{\frac{1}{2}},
```

 where $\tilde\omega$ is the corresponding element of the signature $\tilde\Omega$ of the definition equation {math:numref}`eq:hsvd`.

### Controlling and Monitoring Convergence

Similarly to the case of eigensolvers, in `SVD` the number of iterations carried out by the solver can be determined with `SVDGetIterationNumber`, and the tolerance and maximum number of iterations can be set with `SVDSetTolerances`. Also, convergence can be monitored with command-line keys `-svd_monitor`, `-svd_monitor_all`, `-svd_monitor_conv`, or graphically with `-svd_monitor draw::draw_lg`, or alternatively with `-svd_monitor_all draw::draw_lg`. See section [](#sec:monitor) for additional details.

:::{table} Available possibilities for the convergence criterion, with $Z=A$ in the standard SVD or as defined in equation {math:numref}`eq:qr` for the GSVD.
:name: tab:svdconvergence

 |Convergence criterion       |`SVDConv`         |Command line key   |Error bound
 |----------------------------|------------------|-------------------|-------------------
 |Absolute                    |`SVD_CONV_ABS`    |`-svd_conv_abs`    |$\|r\|$
 |Relative to eigenvalue      |`SVD_CONV_REL`    |`-svd_conv_rel`    |$\|r\|/\|\lambda\|$
 |Relative to matrix norms    |`SVD_CONV_NORM`   |`-svd_conv_norm`   |$\|r\|/\|Z\|$
 |Fixed number of iterations  |`SVD_CONV_MAXIT`  |`-svd_conv_maxit`  |$\infty$
 |User-defined                |`SVD_CONV_USER`   |`-svd_conv_user`   |user function

:::

As in the case of eigensolvers, the user can choose different convergence tests, based on an error bound obtained from the computed residual norm, $\|r\|$. Table [](#tab:svdconvergence) lists the available options. It is worth mentioning the `SVD_CONV_MAXIT` convergence criterion, which is a bit special. With this criterion, the solver will not compute any residual norms and will stop with a successful status when the maximum number of iterations is reached. This is intended for the `SVDRANDOMIZED` solver in cases when a low-rank approximation of a matrix needs to be computed instead of accurate singular vectors.

### Viewing the Solution

There is support for different kinds of viewers for the solution, as in the case of eigensolvers. One can for instance use `-svd_view_values`, `-svd_view_vectors`, `-svd_error_relative`, or `-svd_converged_reason`. See the description in section [](#sec:epsviewers).

```{rubric} Footnotes
```

[^hsvd]: The exception is in the HSVD, where $u_j$ is normalized so that $U^*\Omega U=\tilde\Omega$.

```{eval-rst}
.. bibliography::
   :filter: docname in docnames
```
