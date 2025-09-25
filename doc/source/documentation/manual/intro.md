(cap:int)=
# Getting Started

SLEPc, the *Scalable Library for Eigenvalue Problem Computations*, is a software library for the solution of large sparse eigenvalue problems on parallel computers.

Together with linear systems of equations, eigenvalue problems are a very important class of linear algebra problems. The need for the numerical solution of these problems arises in many situations in science and engineering, in problems associated with stability and vibration analysis in practical applications. These are usually formulated as large sparse eigenproblems.

Computing eigenvalues is essentially more difficult than solving linear systems of equations. This has resulted in a very active research activity in the area of computational methods for eigenvalue problems in the last years, with many remarkable achievements. However, these state-of-the-art methods and algorithms are not easily transferred to the scientific community, and, apart from a few exceptions, most user still rely on simpler, well-established techniques.

The reasons for this situation are diverse. First, new methods are increasingly complex and difficult to implement and therefore robust implementations must be provided by computational specialists, for example as software libraries. The development of such libraries requires to invest a lot of effort but sometimes they do not reach normal users due to a lack of awareness.

In the case of eigenproblems, using libraries is not straightforward. It is usually recommended that the user understands how the underlying algorithm works and typically the problem is successfully solved only after several cycles of testing and parameter tuning. Methods are often specific for a certain class of eigenproblems and this leads to an explosion of available algorithms from which the user has to choose. Not all these algorithms are available in the form of software libraries, even less frequently with parallel capabilities.

Another difficulty resides in how to represent the operator matrix. Unlike in dense methods, there is no widely accepted standard for basic sparse operations in the spirit of BLAS. This is due to the fact that sparse storage is more complicated, admitting of more variation, and therefore less standardized. For this reason, sparse libraries have an added level of complexity. This holds even more so in the case of parallel distributed-memory programming, where the data of the problem have to be distributed across the available processors.

The first implementations of algorithms for sparse matrices required a prescribed storage format for the sparse matrix, which is an obvious limitation. An alternative way of matrix representation is by means of a user-provided subroutine for the matrix-vector product. Apart from being format-independent, this approach allows the solution of problems in which the matrix is not available explicitly. The drawback is the restriction to a fixed-prototype subroutine.

A better solution for the matrix representation problem is the well-known reverse communication interface, a technique that allows the development of iterative methods disregarding the implementation details of various operations. Whenever the iterative method subroutine needs the results of one of the operations, it returns control to the user's subroutine that called it. The user's subroutine then invokes the module that performs the operation. The iterative method subroutine is invoked again with the results of the operation.

Several libraries with any of the interface schemes mentioned above are publicly available. For a survey of such software see the SLEPc Technical Report {cite:t}`str-6`, "A Survey of Software for Sparse Eigenvalue Problems", and references therein. Some of the most recent libraries are even prepared for parallel execution (some of them can be used from within SLEPc, see section [](#sec:wrap)). However, they still lack some flexibility or require too much programming effort from the user, especially in the case that the eigensolution requires to employ advanced techniques such as spectral transformations or preconditioning.

A further obstacle appears when these libraries have to be used in the context of large software projects carried out by inter-disciplinary teams. In this scenery, libraries must be able to interoperate with already existing software and with other libraries. In order to cope with the complexity associated with such projects, libraries must be designed carefully in order to overcome hurdles such as different storage formats or programming languages. In the case of parallel software, care must be taken also to achieve portability to a wide range of platforms with good performance and still retain flexibility and usability.

## SLEPc and PETSc

The SLEPc library is an attempt to provide a solution to the situation described in the previous paragraphs. It is intended to be a general library for the solution of eigenvalue problems that arise in different contexts, covering standard and generalized problems, both Hermitian and non-Hermitian, with either real or complex arithmetic. Issues such as usability, portability, efficiency and interoperability are addressed, and special emphasis is put on flexibility, providing data-structure neutral implementations and multitude of run-time options. SLEPc offers a growing number of eigensolvers as well as interfaces to integrate well-established eigenvalue packages such as ARPACK. In addition to the linear eigenvalue problem, SLEPc also includes other solver classes for nonlinear eigenproblems, SVD and the computation of the action of a matrix function.

SLEPc is based on PETSc, the Portable, Extensible Toolkit for Scientific Computation {cite:p}`Balay:PUM`, and, therefore, a large percentage of the software complexity is avoided since many PETSc developments are leveraged, including matrix storage formats and linear solvers, to name a few. SLEPc focuses on high level features for eigenproblems, structured around a few object classes as described below.

PETSc uses modern programming paradigms to ease the development of large-scale scientific application codes in Fortran, C, and C++ and provides a powerful set of tools for the numerical solution of partial differential equations and related problems on high-performance computers. Its approach is to encapsulate mathematical algorithms using object-oriented programming techniques, which allow to manage the complexity of efficient numerical message-passing codes. All the PETSc software is free and used around the world in a variety of application areas.

The design philosophy is not to try to completely conceal parallelism from the application programmer. Rather, the user initiates a combination of sequential and parallel phases of computations, but the library handles the detailed message passing required during the coordination of computations. Some of the design principles are described in {cite:p}`Balay:1997:EMP`.

PETSc is built around a variety of data structures and algorithmic objects. The application programmer works directly with these objects rather than concentrating on the underlying data structures. Each component manipulates a particular family of objects (for instance, vectors) and the operations one would like to perform on the objects. The three basic abstract data objects are index sets, vectors and matrices. Built on top of this foundation are various classes of solver objects, which encapsulate virtually all information regarding the solution procedure for a particular class of problems, including the local state and various options such as convergence tolerances, etc.

SLEPc can be considered an extension of PETSc providing all the functionality necessary for the solution of eigenvalue problems. Figure [](#fig:slepc) shows a diagram of all the different objects included in PETSc (on the left) and those added by SLEPc (on the right). PETSc is a prerequisite for SLEPc and users should be familiar with basic concepts such as vectors and matrices in order to use SLEPc. Therefore, together with this manual we recommend to use the PETSc Users Manual {cite:p}`Balay:PUM`.

```{figure} ../../_static/images/manual/svg/fig-slepc.svg
:alt: Numerical components of PETSc and SLEPc
:name: fig:slepc

Numerical components of PETSc and SLEPc
```

Each of these components consists of an abstract interface (simply a set of calling sequences) and one or more implementations using particular data structures. Both PETSc and SLEPc are written in C, which lacks direct support for object-oriented programming. However, it is still possible to take advantage of the three basic principles of object-oriented programming to manage the complexity of such large packages. PETSc uses data *encapsulation* in both vector and matrix data objects. Application code accesses data through function calls. Also, all the operations are supported through *polymorphism*. The user calls a generic interface routine, which then selects the underlying routine that handles the particular data structure. Finally, PETSc also uses *inheritance* in its design. All the objects are derived from an abstract base object. From this fundamental object, an abstract base object is defined for each PETSc object ({external:doc}`Mat`, {external:doc}`Vec` and so on), which in turn has a variety of instantiations that, for example, implement different matrix storage formats.

PETSc/SLEPc provide clean and effective codes for the various phases of solving PDEs, with a uniform approach for each class of problems. This design enables easy comparison and use of different algorithms (for example, to experiment with different Krylov subspace methods, preconditioners, or eigensolvers). Hence, PETSc, together with SLEPc, provides a rich environment for modeling scientific applications as well as for rapid algorithm design and prototyping.

Options can be specified by means of calls to subroutines in the source code and also as command-line arguments. Runtime options allow the user to test different tolerances, for example, without having to recompile the program. Also, since PETSc provides a uniform interface to all of its linear solvers ---the Conjugate Gradient, GMRES, etc.--- and a large family of preconditioners ---block Jacobi, overlapping additive Schwarz, etc.---, one can compare several combinations of method and preconditioner by simply specifying them at execution time. SLEPc shares this good property.

The components enable easy customization and extension of both algorithms and implementations. This approach promotes code reuse and flexibility, and separates the issues of parallelism from the choice of algorithms. The PETSc infrastructure creates a foundation for building large-scale applications.

{#sec:inst}
## Installation

This section describes SLEPc's installation procedure. Previously to the installation of SLEPc, the system must have an appropriate version of PETSc installed. Compatible versions of PETSc and SLEPc are those with coincident major and minor version number, the third number (patch level) being irrelevant for this. For instance, SLEPc {{env.config.release}} may be built with PETSc {{env.config.release}}. Also note that, if using git repositories, both PETSc and SLEPc must be either release versions or development versions, so make sure you select the appropriate branch in both repositories (`git checkout release` or `git checkout main`).

The installation process for SLEPc is very similar to PETSc, with two stages: configuration and compilation. SLEPc's configuration is much simpler because most of the configuration information is taken from PETSc, including compiler options and scalar type (real or complex). See section [](#sec:opt-inst) for a discussion of options that are most relevant for SLEPc. Several configurations can coexist in the same directory tree, so that for instance one can have SLEPc libraries compiled with real scalars as well as with complex scalars. This is explained in section [](#sec:mult-inst). Also, system-based installation is also possible with the `--prefix` option, as discussed in section [](#sec:prefix-inst).

{#sec:std-inst}
### Standard Installation

The basic steps for the installation are described next. Note that prior to these steps, optional packages must have been installed. If any of these packages is installed afterwards, reconfiguration and recompilation is necessary. Refer to sections [](#sec:opt-inst) and [](#sec:wrap) for details about installation of some of these packages.

1.  Unbundle the distribution file with

    ```{parsed-literal}
    tar xzf slepc-{{env.config.release}}.tar.gz
    ```

    or an equivalent command. This will create a directory and unpack the software there.

2.  Set the environment variable `SLEPC_DIR` to the full path of the SLEPc home directory. For example, under the `bash` shell:

    ```{parsed-literal}
    export SLEPC_DIR=/home/username/slepc-{{env.config.release}}
    ```

    In addition, the variables `PETSC_DIR` and `PETSC_ARCH` must also be set appropriately, e.g.

    ```{parsed-literal}
    export PETSC_DIR=/home/username/petsc-{{env.config.release}}
    export PETSC_ARCH=arch-darwin-c-debug
    ```

    The rationale for `PETSC_ARCH` is explained in section [](#sec:mult-inst) (see section [](#sec:prefix-inst) for a case in which `PETSC_ARCH` is not required).


(step-config)=
3.  Change to the SLEPc directory and run the configuration script:

    ```{code} console
    $ cd $SLEPC_DIR
    $ ./configure
    ```

4.  If the configuration was successful, build the libraries:

    ```{code} console
    $ make
    ```

5.  After the compilation, try running some test examples with

    ```{code} console
    $ make check
    ```

    Examine the output for any obvious errors or problems.

{#sec:opt-inst}
### Configuration Options

Several options are available in SLEPc's configuration script. To see all available options, type `./configure --help`.

In SLEPc, configure options have the following purposes:

-   Specify a directory for prefix-based installation, as explained in section [](#sec:prefix-inst).

-   Enable external eigensolver packages. For example, to use ARPACK, specify the following options (with the appropriate paths):

    ```{code} console
    $ ./configure --with-arpack-dir=/usr/software/ARPACK
    ```

    Some of the external packages also support the `--download-xxxx` option. Section [](#sec:wrap) provides more details related to use of external libraries.

Additionally, PETSc's configuration script provides a very long list of options that are relevant to SLEPc. Here is a list of options that may be useful. Note that these are options of PETSc that apply to both PETSc and SLEPc, in such a way that it is not possible to, e.g., build PETSc without debugging and SLEPc with debugging.

-   Add `--with-scalar-type=complex` to build complex scalar versions of all libraries. See below a note related to complex scalars.

-   Build single precision versions with `--with-precision=single`. In most applications, this can achieve a significant reduction of memory requirements, and a moderate reduction of computing time. Also, quadruple precision (128-bit floating-point representation) is also available using `--with-precision=__float128` on systems with GNU compilers (`gcc-4.6` or later).

-   Enable use from Fortran. By default, PETSc's configure looks for an appropriate Fortran compiler. If not required, this can be disabled: `--with-fc=0`. If required but not correctly detected, the compiler to be used can be specified with a configure option. It is also possible to configure with a Fortran compiler but do not build Fortran interfaces of PETSc and SLEPc, with `--with-fortran-bindings=0`.

-   If not detected, use `--with-blas-lapack-lib` to specify the location of BLAS and LAPACK. If SLEPc's configure complains about some missing LAPACK subroutines, reconfigure PETSc with option `--download-f2cblaslapack`.

-   Enable external libraries that provide direct linear solvers or preconditioners, such as MUMPS, hypre, or SuperLU; for example, `--download-mumps`. These are especially relevant for SLEPc in the case that a spectral transformation is used, see chapter [](#cap:st).

-   Add `--with-64-bit-indices=1` to use 8 byte integers (`long long`) for indexing in vectors and matrices. This is only needed when working with over roughly 2 billion unknowns.

-   Build static libraries, `--with-shared-libraries=0`. This is generally not recommended, since shared libraries produce smaller executables and the run time overhead is small.

-   Error-checking code can be disabled with `--with-debugging=0`, but this is only recommended in production runs of well-tested applications.

-   Enable GPU computing setting `--with-cuda=1` or `--with-hip=1`; see section [](#sec:gpu) for details.

-   The option `--with-mpi=0` allows building PETSc and SLEPc without MPI support (only sequential).

**Note about complex scalar versions**: PETSc supports the use of complex scalars by defining the data type {external:doc}`PetscScalar` either as a real or complex number. This implies that two different versions of the PETSc libraries can be built separately, one for real numbers and one for complex numbers, but they cannot be used at the same time. SLEPc inherits this property. In SLEPc it is not possible to completely separate real numbers and complex numbers because the solution of non-symmetric real-valued eigenvalue problems may be complex. SLEPc has been designed trying to provide a uniform interface to manage all the possible cases. However, there are slight differences between the interface in each of the two versions. In this manual, differences are clearly identified.

{#sec:mult-inst}
### Installing Multiple Configurations in a Single Directory Tree

Often, it is necessary to build two (or more) versions of the libraries that differ in a few configuration options. For instance, versions for real and complex scalars, or versions for double and single precision, or versions with debugging and optimized. In a standard installation, this is handled by building all versions in the same directory tree, as explained below, so that source code is not replicated unnecessarily. In contrast, in prefix-based installation where source code is not present, the issue of multiple configurations is handled differently, as explained in section [](#sec:prefix-inst).

In a standard installation, the different configurations are identified by a unique name that is assigned to the environment variable `PETSC_ARCH`. Let us illustrate how to set up PETSc with two configurations. First, set a value of `PETSC_ARCH` and proceed with the installation of the first one:

```{code} console
$ cd $PETSC_DIR
$ export PETSC_ARCH=arch-linux-gnu-c-debug-real
$ ./configure --with-scalar-type=real
$ make all
```

Note that if `PETSC_ARCH` is not given a value, PETSc suggests one for us. After this, a subdirectory named `$PETSC_ARCH` is created within `$PETSC_DIR`, that stores all information associated with that configuration, including the built libraries, configuration files, automatically generated source files, and log files. For the second configuration, proceed similarly:

```{code} console
$ cd $PETSC_DIR
$ export PETSC_ARCH=arch-linux-gnu-c-debug-complex
$ ./configure --with-scalar-type=complex
$ make all
```

The value of `PETSC_ARCH` in this case must be different than the previous one. It is better to set the value of `PETSC_ARCH` explicitly, because the name suggested by `configure` may coincide with an existing value, thus overwriting a previous configuration. After successful installation of the second configuration, two `$PETSC_ARCH` directories exist within `$PETSC_DIR`, and the user can easily choose to build his/her application with either configuration by simply changing the value of `PETSC_ARCH`.

The configuration of two versions of SLEPc in the same directory tree is very similar. The only important restriction is that the value of `PETSC_ARCH` used in SLEPc must exactly match an existing PETSc configuration, that is, a directory `$PETSC_DIR/$PETSC_ARCH` must exist.

{#sec:prefix-inst}
### Prefix-based Installation

Both PETSc and SLEPc allow for prefix-based installation. This consists in specifying a directory to which the files generated during the building process are to be copied.

In PETSc, if an installation directory has been specified during configuration (with option `--prefix` in step [configuration](#step-config) of section [](#sec:std-inst)), then after building the libraries the relevant files are copied to that directory by typing

```{code} console
$ make install
```

This is useful for building as a regular user and then copying the libraries and include files to the system directories as root.

To be more precise, suppose that the configuration was done with `--prefix=/opt/petsc-x.x-linux-gnu-c-debug`. Then, `make install` will create directory `/opt/petsc-x.x-linux-gnu-c-debug` if it does not exist, and several subdirectories containing the libraries, the configuration files, and the header files. Note that the source code files are not copied, nor the documentation, so the size of the installed directory will be much smaller than the original one. For that reason, it is no longer necessary to allow for several configurations to share a directory tree. In other words, in a prefix-based installation, variable `PETSC_ARCH` loses significance and must be unset. To maintain several configurations, one should specify different prefix directories, typically with a name that informs about the configuration options used.

In order to prepare a prefix-based installation of SLEPc that uses a prefix-based installation of PETSc, start by setting the appropriate value of `PETSC_DIR`. Then, run SLEPc's configure with a prefix directory.

```{parsed-literal}
export PETSC_DIR=/opt/petsc-{{env.config.release}}-linux-gnu-c-debug
unset PETSC_ARCH
cd $SLEPC_DIR
./configure --prefix=/opt/slepc-{{env.config.release}}-linux-gnu-c-debug
make
make install
export SLEPC_DIR=/opt/slepc-{{env.config.release}}-linux-gnu-c-debug
```

Note that the variable `PETSC_ARCH` has been unset before SLEPc's configure. SLEPc will use a temporary arch name during the build (this temporary arch is named `installed-arch-xxx`, where the `arch-xxx` string represents the configuration of the installed PETSc version). Although it is not a common case, it is also possible to configure SLEPc without prefix, in which case the `PETSC_ARCH` variable must still be empty and the arch directory `installed-xxx` is picked automatically (it is hardwired in file `$SLEPC_DIR/lib/slepc/conf/slepcvariables`). The combination PETSc without prefix and SLEPc with prefix is also allowed, in which case `PETSC_ARCH` should not be unset.

## Running SLEPc Programs

Before using SLEPc, the user must first set the environment variable `SLEPC_DIR`, indicating the full path of the directory containing SLEPc. For example, under the `bash` shell, a command of the form

```{parsed-literal}
export SLEPC_DIR=/software/slepc-{{env.config.release}}
```

can be placed in the user's `.bashrc` file. The `SLEPC_DIR` directory can be either a standard installation SLEPc directory, or a prefix-based installation directory, see section [](#sec:prefix-inst). In addition, the user must set the environment variables required by PETSc, that is, `PETSC_DIR`, to indicate the full path of the PETSc directory, and `PETSC_ARCH` to specify a particular architecture and set of options. Note that `PETSC_ARCH` should not be set in the case of prefix-based installations.

All PETSc programs use the MPI (Message Passing Interface) standard for message-passing communication {cite:p}`MPI-Forum:1994:MMI`. Thus, to execute SLEPc programs, users must know the procedure for launching MPI jobs on their selected computer system(s). Usually, the `mpiexec` command can be used to initiate a program as in the following example that uses eight processes:

```{code} console
$ mpiexec -n 8 slepc_program [command-line options]
```

Note that MPI may be deactivated during configuration of PETSc, if one wants to run only serial programs in a laptop, for example.

All PETSc-compliant programs support the use of the `-h` or `-help` option as well as the `-v` or `-version` option. In the case of SLEPc programs, specific information for SLEPc is also displayed.

## Writing SLEPc Programs

Most SLEPc programs begin with a call to `SlepcInitialize`

```{code} c
SlepcInitialize(int *argc,char ***argv,char *file,char *help);
```

which initializes SLEPc, PETSc and MPI. This subroutine is very similar to {external:doc}`PetscInitialize`, and the arguments have the same meaning. In fact, internally `SlepcInitialize` calls {external:doc}`PetscInitialize`.

After this initialization, SLEPc programs can use communicators defined by PETSc. In most cases users can employ the communicator {external:doc}`PETSC_COMM_WORLD` to indicate all processes in a given run and {external:doc}`PETSC_COMM_SELF` to indicate a single process. MPI provides routines for generating new communicators consisting of subsets of processes, though most users rarely need to use these features. SLEPc users need not program much message passing directly with MPI, but they must be familiar with the basic concepts of message passing and distributed memory computing.

All SLEPc programs should call `SlepcFinalize` as their final (or nearly final) statement

```{code} c
SlepcFinalize();
```

This routine handles operations to be executed at the conclusion of the program, and calls {external:doc}`PetscFinalize` if `SlepcInitialize` began PETSc.

**Note to Fortran Programmers**: In this manual all the examples and calling sequences are given for the C/C++ programming languages. However, Fortran programmers can use most of the functionality of SLEPc and PETSc from Fortran, with only minor differences in the user interface. For instance, the two functions mentioned above have their corresponding Fortran equivalent:

```{code} Fortran
call SlepcInitialize(file,ierr)
call SlepcFinalize(ierr)
```

Section [](#sec:fortran) provides a summary of the differences between using SLEPc from Fortran and C/C++, as well as a complete Fortran example.

{#sec:simpleex}
### Simple SLEPc Example

A simple example is listed next that solves an eigenvalue problem associated with the one-dimensional Laplacian operator discretized with finite differences. This example can be found in `${SLEPC_DIR}/src/eps/tutorials/ex1.c`. Following the code we highlight a few of the most important parts of this example.

```{include} ex1.c
:name: ex1.c
:code: c
```

#### Include Files.

The C/C++ include files for SLEPc should be used via statements such as

```{code} c
#include <slepceps.h>
```

where `slepceps.h` is the include file for the `EPS` component. Each SLEPc program must specify an include file that corresponds to the highest level SLEPc objects needed within the program; all of the required lower level include files are automatically included within the higher level files. For example, `slepceps.h` includes `slepcst.h` (spectral transformations), and `slepcsys.h` (base SLEPc file). Some PETSc header files are included as well, such as `PETScksp.h`. The SLEPc include files are located in the directory `${SLEPC_DIR}/include`.

#### The Options Database.

All the PETSc functionality related to the options database is available in SLEPc. This allows the user to input control data at run time very easily. In this example, the call {external:doc}`PetscOptionsGetInt``(NULL,NULL,"-n",&n,NULL)` checks whether the user has provided a command line option to set the value of `n`, the problem dimension. If so, the variable `n` is set accordingly; otherwise, `n` remains unchanged.

#### Vectors and Matrices.

Usage of matrices and vectors in SLEPc is exactly the same as in PETSc. The user can create a new parallel or sequential matrix, `A`, which has `M` global rows and `N` global columns, with

```{code} c
MatCreate(MPI_Comm comm,Mat *A);
MatSetSizes(Mat A,PetscInt m,PetscInt n,PetscInt M,PetscInt N);
MatSetFromOptions(Mat A);
```

where the matrix format can be specified at runtime. The example creates a matrix, sets the nonzero values with {external:doc}`MatSetValues` and then assembles it.

#### Eigensolvers.

Usage of eigensolvers is very similar to other kinds of solvers provided by PETSc. After creating the matrix (or matrices) that define the problem, $Ax = kx$ (or $Ax=kBx$), the user can then use `EPS` to solve the system with the following sequence of commands: `EPSCreate` `EPSSetOperators` `EPSSetProblemType` `EPSSetFromOptions` `EPSSolve` `EPSDestroy` `EPSGetConverged` `EPSGetEigenpair`

```{code} c
EPSCreate(MPI_Comm comm,EPS *eps);
EPSSetOperators(EPS eps,Mat A,Mat B);
EPSSetProblemType(EPS eps,EPSProblemType type);
EPSSetFromOptions(EPS eps);
EPSSolve(EPS eps);
EPSGetConverged(EPS eps,PetscInt *nconv);
EPSGetEigenpair(EPS eps,PetscInt i,PetscScalar *kr,PetscScalar *ki,Vec xr,Vec xi);
EPSDestroy(EPS *eps);
```

The user first creates the `EPS` context and sets the operators associated with the eigensystem as well as the problem type. The user then sets various options for customized solution, solves the problem, retrieves the solution, and finally destroys the `EPS` context. Chapter [](#cap:eps) describes in detail the `EPS` package, including the options database that enables the user to customize the solution process at runtime by selecting the solution algorithm and also specifying the convergence tolerance, the number of eigenvalues, the dimension of the subspace, etc.

#### Spectral Transformation.

In the example program shown above there is no explicit reference to spectral transformations. However, an `ST` object is handled internally so that the user is able to request different transformations such as shift-and-invert. Chapter [](#cap:st) describes the `ST` package in detail.

#### Error Checking.

All SLEPc routines return an integer indicating whether an error has occurred during the call. The error code is set to be nonzero if an error has been detected; otherwise, it is zero. The PETSc macro {external:doc}`PetscCall``(...)` checks the return value and calls the PETSc error handler upon error detection. {external:doc}`PetscCall``(...)` should be used in all subroutine calls to enable a complete error traceback. See the PETSc documentation for full details.

### Writing Application Codes with SLEPc

Several example programs demonstrate the software usage and can serve as templates for developing custom applications. They are scattered throughout the SLEPc directory tree, in particular in the `tutorials` directories under each class subdirectory.

To write a new application program using SLEPc, we suggest the following procedure:

1.  Install and test SLEPc according to the instructions given in the documentation.

2.  Copy the SLEPc example that corresponds to the class of problem of interest (e.g., singular value decomposition).

3.  Create a makefile as explained below, compile and run the example program.

4.  Use the example program as a starting point for developing a custom code.

Application program makefiles can be set up very easily just by including one file from the SLEPc makefile system. All the necessary PETSc definitions are loaded automatically. The following sample makefile illustrates how to build C and Fortran programs:

```{code} c
default: ex1

include ${SLEPC_DIR}/lib/slepc/conf/slepc_common

ex1: ex1.o
        -${CLINKER} -o ex1 ex1.o ${SLEPC_EPS_LIB}
        ${RM} ex1.o

ex1f: ex1f.o
        -${FLINKER} -o ex1f ex1f.o ${SLEPC_EPS_LIB}
        ${RM} ex1f.o
```

```{rubric} Footnotes
```

```{eval-rst}
.. bibliography::
   :filter: docname in docnames
```
