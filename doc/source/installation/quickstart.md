# Basic Installation Instructions

The following is a quick-start guide for installing SLEPc. For further details, the user is referred to the SLEPc Users Manual.

Previously to the installation of SLEPc, the system must have an appropriate version of PETSc installed (see the {{'[PETSc installation documentation](https://petsc.org/{}/install/)'.format(branch)}} for details).

The installation process for SLEPc is very similar to PETSc, with two main stages: configuration and compilation. SLEPc configuration is much simpler because most of the configuration information is taken from PETSc, including compiler options and scalar type (real or complex). Several configurations can coexist in the same directory tree, being selected by different values of `PETSC_ARCH`, so that one can, for instance, have a SLEPc compiled with real scalars and another one with complex scalars.

The main steps for the installation are:

  1. Unbundle the distribution file **[slepc-{{env.config.release}}.tar.gz](https://slepc.upv.es/download/distrib/slepc-{{env.config.release}}.tar.gz)** with a usual command such as

```{parsed-literal}
$ tar xzf slepc-{{env.config.release}}.tar.gz
```

This will create a directory and unpack the software there.

  2. Set the environment variable `SLEPC_DIR` to the full path of the SLEPc home directory, for example,

```{parsed-literal}
$ export SLEPC_DIR=/home/username/slepc-{{env.config.release}}
```

In addition to this variable, `PETSC_DIR` and `PETSC_ARCH` must also be set appropriately.

  3. In the SLEPc directory, execute

```{code} console
$ ./configure
```

:::{note}
In order to enable external packages (see below), this command must be run with appropiate options. To see all available options use `./configure --help`
:::

  4. In the SLEPc home directory, type

```{code} console
$ make
```

  5. Optionally, if an installation directory has been specified during configuration (with option `--prefix` in step 3 above), then type

```{code} console
$ make install
```

This is useful for building as a regular user and then copying the libraries and include files to the system directories as root.

  6. If the installation went smoothly, then try running some test examples with

```{code} console
$ make check
```

Examine the output for any obvious errors or problems.

## Optional Software

SLEPc provides an interface to several software packages. These should be installed before installing SLEPc. These packages are not developed, maintained, or supported by the SLEPc team; we merely provide an interface to them. To integrate one of these libraries in SLEPc:

  * First install the external package following its instructions. Make sure you use the same compilers and MPI that you plan to use with PETSc/SLEPc.
  * Enable the utilization of the external software from SLEPc by adding specific command-line parameters when executing `configure`. For example, to use ARPACK, specify the following options (with the appropriate paths):

```{code} console
$ ./configure --with-arpack-dir=/usr/software/ARPACK
```

  * Build the SLEPc libraries.

SLEPc currently interfaces to the following libraries:

  * [ARPACK](https://github.com/opencollab/arpack-ng) (Implicitly Restarted Arnoldi/Lanczos solver).
  * [PRIMME](http://www.cs.wm.edu/~andreas/software) (PReconditioned Iterative MultiMethod Eigensolver).
  * [BLOPEX](https://github.com/lobpcg/blopex) (Block Locally Optimal Preconditioned Eigenvalue Xolvers).

Additional information about these packages can be found in the SLEPc Users Manual.
