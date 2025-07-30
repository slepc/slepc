# Download

## Release Version {{env.config.release}}

### Git

The source repository is hosted at [gitlab.com](https://gitlab.com/slepc/slepc) and uses the `git` version control system. The development process follows the methodology used by PETSc, with two integration branches (`main` and `release`), together with various feature branches. One can expect SLEPc's `release` branch to be synchronized with PETSc's `release` branch.

In order to get the release version of SLEPc, you need to clone the repository indicating `release` as the branch name:

```{code} console
$ git clone -b release https://gitlab.com/slepc/slepc
```

At a later time, the repository can be refreshed simply by:

```{code} console
$ git pull
```
The `release` branch tracks the latest release version, including all the patches. It is also possible to check out a specific version with the tag name, for instance:

%```{code} console
```{parsed-literal}
$ git checkout v{{env.config.release}}
```

### Tarball

Alternatively, the following downloads correspond to the latest available release.

Description                                               |  File                                                                          |  MD5 checksum
---                                                       |  ---                                                                           |  ---
SLEPc distribution file (source code only)                |  **[slepc-{{env.config.release}}.tar.gz](https://slepc.upv.es/download/distrib/slepc-{{env.config.release}}.tar.gz)**                     |  XXXmd5sumXXX fe40ac5cd967044f9d2317a12482e797
SLEPc distribution file (source code with documentation)  |  **[slepc-with-docs-{{env.config.release}}.tar.gz](https://slepc.upv.es/download/distrib/slepc-with-docs-{{env.config.release}}.tar.gz)** | XXXmd5sumXXX 3af9b2dcb142f2122ecfe803c28fcaee
slepc4py distribution file (enables separate install)     |  **[slepc4py-{{env.config.release}}.tar.gz](https://slepc.upv.es/download/distrib/slepc4py-{{env.config.release}}.tar.gz)**               | XXXmd5sumXXX b544b9b9011055ef7592fa84154744c3

There are no separate patch files, the current fixes are included in the tar file. Patches are documented at: [slepc-{{env.config.version}} changelog](https://gitlab.com/slepc/slepc/commits/release)

:::{note}
slepc4py source (without documentation) is already included in the SLEPc tarball.
:::

:::{note}
Major releases are announced via the [slepc-announce](../contact/mail_list) mailing list, while patch releases are announced via the RSS news feed.
:::

## Development Version

%The development repository is hosted at [gitlab.com](https://gitlab.com/slepc/slepc) and uses the `git` version control system. Note that using the development version requires also a development version of PETSc.

%The development process follows the methodology used by PETSc, with two integration branches (`main` and `release`), together with various feature branches. One can expect SLEPc's `main` branch to be synchronized with PETSc's `main` branch.

In order to get the development version of SLEPc, you need to clone the repository:

```{code} console
$ git clone https://gitlab.com/slepc/slepc
```

Or, if already cloned, switch to the `main` branch:

```{code} console
$ git switch main
```

At a later time, the repository can be refreshed simply by:

```{code} console
$ git pull
```

<div style="float: right; width: 45%;">
<iframe src='https://www.openhub.net/p/slepc/widgets/project_basic_stats' style='height: 225px; width: 350px; border: none'></iframe>
</div>

The source code in the repository can be browsed at [gitlab.com](https://gitlab.com/slepc/slepc). Some statistics related to source code development can be found at [openhub.net](https://www.openhub.net/p/slepc).

Additional information can be found at [PETSc Developer's Documentation](https://petsc.org/release/developers/).

Users of the development version may also want to check the [nightly tests](https://gitlab.com/slepc/slepc/-/pipeline_schedules) for various builds of SLEPc's `main` branch with PETSc's `main` branch, as well as the source [coverage](https://slepc.upv.es/coverage/) of these tests.

## Previous Versions

Distribution files:
[[3.0.0]](https://slepc.upv.es/download/distrib/slepc-3.0.0-p7.tgz)
[[3.1]](https://slepc.upv.es/download/distrib/slepc-3.1-p6.tgz)
[[3.2]](https://slepc.upv.es/download/distrib/slepc-3.2-p5.tar.gz)
[[3.3]](https://slepc.upv.es/download/distrib/slepc-3.3-p4.tar.gz)
[[3.4]](https://slepc.upv.es/download/distrib/slepc-3.4.4.tar.gz)
[[3.5]](https://slepc.upv.es/download/distrib/slepc-3.5.4.tar.gz)
[[3.6]](https://slepc.upv.es/download/distrib/slepc-3.6.3.tar.gz)
[[3.7]](https://slepc.upv.es/download/distrib/slepc-3.7.4.tar.gz)
[[3.8]](https://slepc.upv.es/download/distrib/slepc-3.8.3.tar.gz)
[[3.9]](https://slepc.upv.es/download/distrib/slepc-3.9.2.tar.gz)
[[3.10]](https://slepc.upv.es/download/distrib/slepc-3.10.2.tar.gz)
[[3.11]](https://slepc.upv.es/download/distrib/slepc-3.11.2.tar.gz)
[[3.12]](https://slepc.upv.es/download/distrib/slepc-3.12.2.tar.gz)
[[3.13]](https://slepc.upv.es/download/distrib/slepc-3.13.4.tar.gz)
[[3.14]](https://slepc.upv.es/download/distrib/slepc-3.14.2.tar.gz)
[[3.15]](https://slepc.upv.es/download/distrib/slepc-3.15.2.tar.gz)
[[3.16]](https://slepc.upv.es/download/distrib/slepc-3.16.3.tar.gz)
[[3.17]](https://slepc.upv.es/download/distrib/slepc-3.17.2.tar.gz)
[[3.18]](https://slepc.upv.es/download/distrib/slepc-3.18.3.tar.gz)
[[3.19]](https://slepc.upv.es/download/distrib/slepc-3.19.2.tar.gz)
[[3.20]](https://slepc.upv.es/download/distrib/slepc-3.20.2.tar.gz)
[[3.21]](https://slepc.upv.es/download/distrib/slepc-3.21.2.tar.gz)
[[3.22]](https://slepc.upv.es/download/distrib/slepc-3.22.2.tar.gz)
[[3.23]](https://slepc.upv.es/download/distrib/slepc-3.23.0.tar.gz)

:::{note}
Users of previous versions are strongly recommended to upgrade to the latest one.
:::

:::{note}
Even older releases can be found [here](older).
:::

