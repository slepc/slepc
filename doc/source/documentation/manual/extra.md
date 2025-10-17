(ch:add)=
# Additional Information

This chapter contains miscellaneous information as a complement to the previous chapters, that can be regarded as less important.

## Supported PETSc Features

SLEPc relies on PETSc for most features that are not directly related to eigenvalue problems. All functionality associated with vectors and matrices as well as systems of linear equations is provided by PETSc. Also, low level details are inherited directly from PETSc. In particular, the parallelism within SLEPc methods is handled almost completely by PETSc's vector and matrix modules.

SLEPc mainly contains high level objects, as depicted in figure [](#fig:slepc). These object classes have been designed and implemented following the philosophy of other high level objects in PETSc. In this way, SLEPc benefits from a number of PETSc's good properties such as the following (see {{'[PETSc Users Guide](https://petsc.org/{}/manual/)'.format(branch)}} for details):

-   Portability and scalability in a wide range of platforms. Different architecture builds can coexist in the same installation. Where available, shared libraries are used to reduce disk space of executable files.

-   Support for profiling of programs:

    -   Display performance statistics with `-log_view`, including also SLEPc's objects. The collected data are *flops*, memory usage and execution times as well as information about parallel performance, for individual subroutines and the possibility of user-defined stages.

    -   Event logging, including user-defined events.

    -   Direct wall-clock timing with {external:doc}`PetscTime`.

    -   Display detailed profile information and trace of events.

-   Convergence monitoring, both textual and graphical.

-   Support for debugging of programs:

    -   Debugger startup and attachment of parallel processes.

    -   Automatic generation of back-traces of the call stack.

    -   Detection of memory leaks.

-   A number of viewers for visualization of data, including built-in graphics capabilities that allow for sparse pattern visualization, graphic convergence monitoring, operator's spectrum visualization and display of regions of the complex plane.

-   Easy handling of runtime options.

-   Support for Fortran programming using Fortran modules. See section [](#sec:fortran) for details.

{#sec:supported}
## Supported Matrix Types

Methods implemented in SLEPc merely require vector operations and matrix-vector products. In PETSc, mathematical objects such as vectors and matrices have an interface that is independent of the underlying data structures. SLEPc manipulates vectors and matrices via this interface and, therefore, it can be used with any of the matrix representations provided by PETSc, including dense, sparse, and symmetric formats, either sequential or parallel.

The above statement must be reconsidered when using `EPS` in combination with `ST`. As explained in chapter [](#ch:st), in many cases the operator associated with a spectral transformation not only consists in pure matrix-vector products but also other operations may be required as well, most notably a linear system solve (see Table [](#tab:op)). In this case, the limitation is that there must be support for the requested operation for the selected matrix representation.

### Shell Matrices

In many applications, the matrices that define the eigenvalue problem are not available explicitly. Instead, the user knows a way of applying these matrices to a vector.

An intermediate case is when the matrices have some block structure and the different blocks are stored separately. There are numerous situations in which this occurs, such as the discretization of equations with a mixed finite-element scheme. An example is the eigenproblem arising in the stability analysis associated with Stokes problems,

```{math}
:label: eq:stokes

\begin{bmatrix}A & C\\C^* & 0\end{bmatrix}\begin{bmatrix}x\\p\end{bmatrix}
=\lambda\begin{bmatrix}B & 0\\0 & 0\end{bmatrix}\begin{bmatrix}x\\p\end{bmatrix}\;\;,
```

where $x$ and $p$ denote the velocity and pressure fields. Similar formulations also appear in many other situations.

In some cases, these problems can be solved by reformulating them as a reduced-order standard or generalized eigensystem, in which the matrices are equal to certain operations of the blocks. These matrices are not computed explicitly to avoid losing sparsity.

All these cases can be easily handled in SLEPc by means of shell matrices. These are matrices that do not require explicit storage of the matrix entries. Instead, the user must provide subroutines for all the necessary matrix operations, typically only the application of the linear operator to a vector.

Shell matrices, also called matrix-free matrices, are created in PETSc with the function {external:doc}`MatCreateShell`. Then, the function {external:doc}`MatShellSetOperation` is used to provide any user-defined shell matrix operations (see the {{'[PETSc Users Guide](https://petsc.org/{}/manual/mat/#application-specific-custom-matrices)'.format(branch)}} for additional details). Several examples are available in SLEPc that illustrate how to solve a matrix-free eigenvalue problem.

In the simplest case, defining matrix-vector product operations (`MATOP_MULT`) is enough for using `EPS` with shell matrices. However, in the case of generalized problems, if matrix $B$ is also a shell matrix then it may be necessary to define other operations in order to be able to solve the linear system successfully, for example `MATOP_GET_DIAGONAL` to use an iterative linear solver with Jacobi preconditioning. On the other hand, if the shift-and-invert `ST` is to be used, then in addition it may also be necessary to define `MATOP_SHIFT` or `MATOP_AXPY` (see section [](#sec:explicit) for discussion).

In the case of `SVD`, both $A$ and $A^*$ are required to solve the problem. So when computing the SVD, the shell matrix needs to have the `MATOP_MULT_TRANSPOSE` operation (or `MATOP_MULT_HERMITIAN_TRANSPOSE` in the case of complex scalars) in addition to `MATOP_MULT`. Alternatively, if $A^*$ is to be built explicitly, `MATOP_TRANSPOSE` is then the required operation. For details, see the manual page for `SVDSetImplicitTranspose`.

{#sec:gpu}
## GPU Computing

Support for graphics processing unit (GPU) computing is included in SLEPc. This is related to section [](#sec:supported) because GPU support in PETSc is based on using special types of {external:doc}`Mat` and {external:doc}`Vec`. GPU support in SLEPc has been tested in all solver classes and most solvers should work, although the performance gain to be expected depends on the particular algorithm. Regarding PETSc, all iterative linear solvers are prepared to run on the GPU, but this is not the case for direct solvers and preconditioners. The user must not expect a spectacular performance boost, but in general moderate gains can be achieved by running the eigensolver on the GPU instead of the CPU (in some cases a 10-fold improvement).

SLEPc currently provides support for NVIDIA GPUs using CUDA[^cuda] as well as AMD GPUs using HIP and ROCm[^rocm].

CUDA provides a C/C++ compiler with CUDA extensions as well as the cuBLAS and cuSPARSE libraries that implement dense and sparse linear algebra operations. For instance, to configure PETSc with GPU support in single precision arithmetic use the following options:

```{code} console
$ ./configure --with-precision=single --with-cuda
```

{external:doc}`VECCUDA` and {external:doc}`MATAIJCUSPARSE` are currently the mechanism in PETSc to run a computation on the GPU. {external:doc}`VECCUDA` is a special type of {external:doc}`Vec` whose array is mirrored in the GPU (and similarly for {external:doc}`MATAIJCUSPARSE`). PETSc takes care of keeping memory coherence between the two copies of the array, and performs the computation on the GPU when possible, trying to avoid unnecessary copies between the host and the device. For maximum efficiency, the user has to make sure that all vectors and matrices are of these types. If they are created in the standard way ({external:doc}`VecCreate` plus {external:doc}`VecSetFromOptions`) then it is sufficient to run the SLEPc program with

```{code} console
$ ./program -vec_type cuda -mat_type aijcusparse
```

Note that the first option is unnecessary if no {external:doc}`Vec` is created in the main program, or if all vectors are created via {external:doc}`MatCreateVecs` from a {external:doc}`MATAIJCUSPARSE`.

For AMD GPUs the procedure is very similar, with HIP providing the compiler and ROCm providing the analogue libraries hipBLAS and hipSPARSE. To configure PETSc with HIP do:

```{code} console
$ ./configure --with-precision=single --with-hip
```

Then the equivalent vector and matrix types are {external:doc}`VECHIP` and {external:doc}`MATAIJHIPSPARSE`, which can be used in the command line with

```{code} console
$ ./program -vec_type hip -mat_type aijhipsparse
```

{#sec:shell}
## Extending SLEPc

Shell matrices, presented in section [](#sec:supported), are a simple mechanism of extensibility, in the sense that the package is extended with new user-defined matrix objects. Once the new matrix has been defined, it can be used by SLEPc in the same way as the rest of the matrices as long as the required operations are provided.

A similar mechanism is available in SLEPc also for extending the system incorporating new spectral transformations (`ST`). This is done by using the `STSHELL` spectral transformation type, in a similar way as shell matrices or shell preconditioners. In this case, the user defines how the operator is applied to a vector and optionally how the computed eigenvalues are transformed back to the solution of the original problem. This tool is intended for simple spectral transformations. For more sophisticated transformations, the user should register a new `ST` type (see below).

The following function:

```{code} c
STShellSetApply(ST,PetscErrorCode(*)(ST,Vec,Vec));
```

has to be invoked after the creation of the `ST` object in order to provide a routine that applies the operator to a vector. And this function:

```{code} c
STShellSetBackTransform(ST,PetscErrorCode(*)(ST,PetscInt,PetscScalar*,PetscScalar*));
```

can be used optionally to specify the routine for the back-transformation of eigenvalues. The two functions provided by the user can make use of any required user-defined information via a context that can be retrieved with `STShellGetContext`. The example program {{'[ex10.c](https://slepc.upv.es/{}/src/eps/tutorials/ex10.c.html)'.format(branch)}} illustrates the use of shell transformations.

SLEPc further supports extensibility by allowing application programmers to code their own subroutines for unimplemented features such as new eigensolvers or new spectral transformations. It is possible to register these new methods to the system and use them as the rest of standard subroutines. For example, to implement a variant of the Subspace Iteration method, one could copy the SLEPc code associated with the `subspace` solver, modify it and register a new `EPS` type with the following line of code:

```{code} c
EPSRegister("newsubspace",EPSCreate_NEWSUB);
```

After this call, the new solver could be used in the same way as the rest of SLEPc solvers, e.g. with `-eps_type newsubspace` in the command line. A similar mechanism is available for registering new types of the other classes.

{#sec:sys}
## Auxiliary Classes

Apart from the main solver classes listed in table [](#tab:modules), SLEPc contains several auxiliary classes:

-   `ST`: Spectral Transformation, fully described in chapter [](#ch:st).

-   `FN`: Mathematical Function, required in application code to represent the constituent functions of the nonlinear operator in split form (chapter [](#ch:nep)), as well as the function to be used when computing the action of a matrix function on a vector (chapter [](#ch:mfn)).

-   `DS`: Direct Solver (or Dense System), can be seen as a wrapper to LAPACK functions used within SLEPc. It is mostly an internal object that need not be called by end users.

-   `BV`: Basis Vectors, provides the concept of a block of vectors that represent the basis of a subspace.

-   `RG`: Region, a way to define a region of the complex plane.

{#sec:fn}
### FN: Mathematical Functions

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

{#sec:bv}
### BV: Basis Vectors

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

{#sec:rg}
### RG: Region

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

## Directory Structure

The directory structure of the SLEPc software is very similar to that of PETSc. The root directory of SLEPc contains the following directories:

- `lib/slepc/conf` - Directory containing the base SLEPc makefile, to be included in application makefiles.
- `config` - SLEPc configuration scripts.
- `doc` - Directory containing the source from which all documentation of SLEPc is generated, including the manual and the website.
- `include` - All include files for SLEPc. The following subdirectories exist:
  - `slepc/finclude` - include files for Fortran programmers.
  - `slepc/private` - include files containing implementation details, for developer use only.
- `share/slepc` - Common files, including:
  - `datafiles` - data files used by some examples.
- `src` - The source code for all SLEPc components, which currently includes:
  - `sys` - system-related routines and auxiliary classes `bv`, `ds`, `fn`, `rg`, `st`.
  - `eps` - eigenvalue problem solver.
  - `svd` - singular value decomposition solver.
  - `pep` - polynomial eigenvalue problem solver.
  - `nep` - nonlinear eigenvalue problem solver.
  - `mfn` - matrix function.
  - `lme` - linear matrix equations.
  - `binding` - source code of slepc4py
- `$PETSC_ARCH` - For each value of `PETSC_ARCH`, a directory exists containing files generated during installation of that particular configuration. Among others, it includes the following subdirectories:
  - `lib` - all the generated libraries.
  - `lib/slepc/conf` - configuration parameters and log files.
  - `include` - automatically generated include files, such as Fortran `*.mod` files.

Each SLEPc source code component directory has the following subdirectories:

- `interface`: The calling sequences for the abstract interface to the components. Code here does not know about particular implementations.
- `impls`: Source code for the different implementations.
- `tutorials`: Example programs intended for learning to use SLEPc.
- `tests`: Example programs used by testing scripts.

{#sec:wrap}
## Wrappers to External Libraries

SLEPc interfaces to several external libraries for the solution of eigenvalue problems. This section provides a short description of each of these packages as well as some hints for using them with SLEPc, including pointers to the respective websites from which the software can be downloaded. The description may also include method-specific parameters, that can be set in the same way as other SLEPc options, either procedurally or via the command-line.

In order to use SLEPc together with an external library such as ARPACK, one needs to do the following.

1.  Install the external software, with the same compilers and MPI that will be used for PETSc/SLEPc.

2.  Enable the utilization of the external software from SLEPc by specifying configure options as explained in section [](#sec:opt-inst).

3.  Build the SLEPc libraries.

4.  Use the runtime option `-eps_type <type>` to select the solver.

Exceptions to the above rule are LAPACK, which should be enabled during PETSc's configuration, and BLOPEX, that must be installed with `--download-blopex` in SLEPc's configure. Other packages also support the download option.

### List of External Libraries

:::{warning}
This list might be incomplete. Check the output of `./configure --help` for other libraries not listed here.
:::

#### LAPACK

<https://www.netlib.org/lapack>

LAPACK {cite:p}`And99` is a software package for the solution of many different dense linear algebra problems, including various types of eigenvalue problems and singular value decompositions. SLEPc explicitly creates the operator matrix in dense form and then the appropriate LAPACK driver routine is invoked. Therefore, this interface should be used only for testing and validation purposes and not in a production code. The operator matrix is created by applying the operator to the columns of the identity matrix.

**Installation**: PETSc already depends on LAPACK. The SLEPc interface to LAPACK can be used directly. If SLEPc's configure script complains about missing LAPACK functions, then reconfigure PETSc with option `--download-f2cblaslapack`.

#### ARPACK

<https://www.arpack.org>

<https://github.com/opencollab/arpack-ng>

ARPACK {cite:p}`Leh98,Mas96` is a software package for the computation of a few eigenvalues and corresponding eigenvectors of a general $n\times n$ matrix $A$. It is most appropriate for large sparse matrices. ARPACK is based upon an algorithmic variant of the Arnoldi process called the Implicitly Restarted Arnoldi Method (IRAM). When the matrix $A$ is symmetric it reduces to a variant of the Lanczos process called the Implicitly Restarted Lanczos Method (IRLM). These variants may be viewed as a synthesis of the Arnoldi/Lanczos process with the Implicitly Shifted QR technique that is suitable for large scale problems. It can be used for standard and generalized eigenvalue problems, both in real and complex arithmetic. It is implemented in Fortran 77 and it is based on the reverse communication interface. A parallel version, PARPACK, is available with support for both MPI and BLACS.

**Installation**: To install ARPACK we recommend using the _arpack-ng_ distribution, available in `github.com`, that supports `configure`+`make` for installation. Also, SLEPc's `configure` allows to download this version automatically via the `--download-arpack` option. It is also possible to configure SLEPc with the serial version of ARPACK. For this, you have to configure PETSc with the option `--with-mpi=0`.

#### PRIMME

<https://www.cs.wm.edu/~andreas/software>

PRIMME {cite:p}`Sta10` is a C library for finding a number of eigenvalues and their corresponding eigenvectors of a real symmetric (or complex Hermitian) matrix. This library provides a multimethod eigensolver, based on Davidson/Jacobi-Davidson. Particular methods include GD+1, JDQMR, and LOBPCG. It supports preconditioning as well as the computation of interior eigenvalues.

**Installation**: Type `make lib` after customizing the file `Make_flags` appropriately. Alternatively, the `--download-primme` option is also available in SLEPc's `configure`.

**Specific options**: Since PRIMME contains preconditioned solvers, the SLEPc interface uses `STPRECOND`, as described in section [](#sec:precond). The SLEPc interface to this package allows the user to specify the maximum allowed block size with the function `EPSPRIMMESetBlockSize` or at run time with the option `-eps_primme_blocksize <size>`. For changing the particular algorithm within PRIMME, use the function `EPSPRIMMESetMethod`. PRIMME also provides a solver for the singular value decomposition that is interfaced in SLEPc's `SVD`, see `SVDPRIMMESetMethod`.

#### EVSL

<https://www-users.cse.umn.edu/~saad/software/EVSL/>

EVSL {cite:p}`Li19` is a sequential library that implements methods for computing all eigenvalues located in a given interval for real symmetric (standard or generalized) eigenvalue problems. Currently SLEPc only supports standard problems.

**Installation**: The option `--download-evsl` is available in SLEPc's configure for easy installation. Alternatively, one can use an already installed version.

#### BLOPEX

<https://github.com/lobpcg/blopex>

BLOPEX {cite:p}`Kny07` is a package that implements the Locally Optimal Block Preconditioned Conjugate Gradient (LOBPCG) method for computing several extreme eigenpairs of symmetric positive generalized eigenproblems. Numerical comparisons suggest that this method is a genuine analog for eigenproblems of the standard preconditioned conjugate gradient method for symmetric linear systems.

**Installation**: In order to use BLOPEX from SLEPc, it necessary to install it during SLEPc's configuration: `./configure --download-blopex`.

**Specific options**: Since BLOPEX contains preconditioned solvers, the SLEPc interface uses `STPRECOND`, as described in section [](#sec:precond).

#### ScaLAPACK

<https://www.netlib.org/scalapack>

ScaLAPACK {cite:p}`Bla97` is a library of high-performance linear algebra routines for parallel distributed memory machines. It contains eigensolvers for dense Hermitian eigenvalue problems, as well as solvers for the (dense) SVD.

**Installation**: For using ScaLAPACK from SLEPc it is necessary to select it during configuration of PETSc.

#### ELPA

<https://elpa.mpcdf.mpg.de/>

ELPA {cite:p}`Auc11` is a high-performance library for the parallel solution of dense symmetric (or Hermitian) eigenvalue problems on distributed memory computers. It uses a ScaLAPACK-compatible matrix distribution.

**Installation**: The SLEPc wrapper to ELPA can be activated at configure time with the option `--download_elpa`, in which case ScaLAPACK support must have been enabled during the configuration of PETSc.

#### KSVD

<https://github.com/ecrc/ksvd/>

KSVD {cite:p}`Suk19` is a high performance software framework for computing a dense SVD on distributed-memory manycore systems. The KSVD solver relies on the polar decomposition (PD) based on the QR Dynamically-Weighted Halley (QDWH) and ZOLO-PD algorithms.

**Installation**: The option `--download-ksvd` is available in SLEPc's configure for easy installation, which in turn requires adding `--download-polar` and `--download-elpa`.

#### ELEMENTAL

<https://github.com/elemental/Elemental>

ELEMENTAL {cite:p}`Pou13` is a distributed-memory, dense and sparse-direct linear algebra package. It contains eigensolvers for dense Hermitian eigenvalue problems, as well as solvers for the SVD.

**Installation**: For using ELEMENTAL from SLEPc it is necessary to select it during configuration of PETSc.

#### FEAST

<https://feast-solver.org/>

FEAST {cite:p}`Pol09` is a numerical library for solving the standard or generalized symmetric eigenvalue problem, and obtaining all the eigenvalues and eigenvectors within a given search interval. It is based on an innovative fast and stable numerical algorithm which deviates fundamentally from the traditional Krylov subspace based iterations or Davidson-Jacobi techniques. The FEAST algorithm takes its inspiration from the density-matrix representation and contour integration technique in quantum mechanics. Latest versions also support non-symmetric problems.

**Installation**: We only support the FEAST implementation included in Intel MKL. For using it from SLEPc it is necessary to configure PETSc with MKL by adding the corresponding option, e.g., `--with-blas-lapack-dir=$MKLROOT`.

**Specific options**: The SLEPc interface to FEAST allows the user to specify the number of contour integration points with the function `EPSFEASTSetNumPoints` or at run time with the option `-eps_feast_num_points <n>`.

#### CHASE

<https://github.com/ChASE-library/ChASE>

CHASE {cite:p}`Win19` is a modern and scalable library based on subspace iteration with polynomial acceleration to solve dense Hermitian (symmetric) algebraic eigenvalue problems, especially solving dense Hermitian eigenproblems arranged in a sequence. Novel to ChASE is the computation of the spectral estimates that enter in the filter and an optimization of the polynomial degree that further reduces the necessary floating-point operations.

**Installation**: Currently, the CHASE interface in SLEPc is based on the MPI version with block-cyclic distribution, i.e., ScaLAPACK matrix storage, so it is necessary to enable ScaLAPACK during configuration of PETSc.

{#sec:fortran}
## Fortran Interface

SLEPc provides an interface for Fortran programmers, very much like PETSc. As in the case of PETSc, there are slight differences between the C and Fortran SLEPc interfaces, due to differences in Fortran syntax. For instance, the error checking variable is the final argument of all the routines in the Fortran interface, in contrast to the C convention of providing the error variable as the routine's return value.

A detailed discussion can be found in the {{'[PETSc Users Guide](https://petsc.org/{}/manual/fortran/)'.format(branch)}}.

The following is a Fortran example. It is the Fortran equivalent of the program given in section [](#sec:simpleex) and can be found in `${SLEPC_DIR}/src/eps/tutorials` (file `ex1f.F90`).

```{literalinclude} /../../src/eps/tutorials/ex1f.F90
:name: ex1f.F90
:language: fortran
:start-at: program main
:end-before: '!/*TEST'
```

```{only} html
<p class="rubric">References</p>
```
```{bibliography}
:filter: docname in docnames
```

```{rubric} Footnotes
```

[^cuda]: <https://developer.nvidia.com/cuda-zone>

[^rocm]: <https://rocm.docs.amd.com>
