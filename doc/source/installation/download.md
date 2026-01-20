# Download

## Release Version {{env.config.release}}

### Git

The source repository is hosted at [gitlab.com](https://gitlab.com/slepc/slepc) and uses the `git` version control system. The development process follows the methodology used by PETSc, with two integration branches (`main` and `release`), together with various feature branches. One can expect SLEPc's `release` branch to be synchronized with PETSc's `release` branch, and similarly for `main`.

In order to get the release version of SLEPc, you need to clone the repository indicating `release` as the branch name:

```{code} console
$ git clone -b release https://gitlab.com/slepc/slepc
```

At a later time, the repository can be refreshed simply by:

```{code} console
$ git pull
```
The `release` branch tracks the latest release version, including all the patches. It is also possible to check out a specific version with the tag name, for instance:

```{parsed-literal}
$ git checkout v{{env.config.release}}
```

### Tarball

Alternatively, the following downloads correspond to the latest available release.

Description                                               |  File                                                                          |  MD5 checksum
---                                                       |  ---                                                                           |  ---
SLEPc distribution file (source code only)                |  {{'**[slepc-{}.tar.gz](https://slepc.upv.es/download/distrib/slepc-{}.tar.gz)**'.format(env.config.release,env.config.release)}} | 3f46f3a7e376ed928cc4ba6281d0f48d
SLEPc distribution file (source code with documentation)  |  {{'**[slepc-with-docs-{}.tar.gz](https://slepc.upv.es/download/distrib/slepc-with-docs-{}.tar.gz)**'.format(env.config.release,env.config.release)}} | 6b8932dca011b7bb88dafaa42e290080
slepc4py distribution file (enables separate install)     |  {{'**[slepc4py-{}.tar.gz](https://slepc.upv.es/download/distrib/slepc4py-{}.tar.gz)**'.format(env.config.release,env.config.release)}}               | a853b642d0d89e826e9459d1d5e33602

There are no separate patch files, the current fixes are included in the tar file. Patches are documented at: [slepc-{{env.config.version}} changelog](https://gitlab.com/slepc/slepc/-/commits/release)

:::{note}
slepc4py source (without documentation) is already included in the SLEPc tarball.
:::

:::{note}
Major releases are announced via the [slepc-announce](../contact/mail_list) mailing list, while patch releases are announced via the RSS news feed.
:::

## Development Version

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

The source code in the repository can be browsed at [gitlab.com](https://gitlab.com/slepc/slepc). Some statistics related to source code development can be found at [openhub.net](https://openhub.net/p/slepc).

Additional information can be found at {{'[PETSc Developer\'s Documentation](https://petsc.org/{}/developers/)'.format(branch)}}.

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
[[3.23]](https://slepc.upv.es/download/distrib/slepc-3.23.3.tar.gz)

:::{warning}
Users of previous versions are strongly recommended to upgrade to the latest one.
:::

:::{admonition} Even Older SLEPc Releases
:class: dropdown

These are very old releases of SLEPc. They are here for historical reasons. They were released with a different license. By downloading the software, you are implicitly agreeing on the following conditions:

> The software is free for academic and research use. This means that if you work in an academic or research institution such as a University or a government laboratory then you can use the software without formally requiring a license.
>
> For commercial use, you need to sign a software license agreement. This includes all users working for a private company, even if the software is going to be used only for in-house research activities. A reasonable testing period is allowed before asking for the license. Normally, we license the software without a fee, but this has to be concerted individually.
>
> In both cases, please note the following. This software is provided 'as is', with absolutely no warranty, expressed or implied. Any use is at your own risk. In no event will the authors be liable for any direct or indirect damages arising in any way out of the use of this software. The user will acknowledge the contribution of the software in any publication of material dependent on its use. The user can modify the code but at no time shall the right or title to all or any part of this software pass to the user. A modified version of the software cannot be redistributed. This software (or a modified version) may not be sold.

Distribution files:
[[2.1.1]](https://slepc.upv.es/download/distrib/slepc-2.1.1.tgz)
[[2.1.5]](https://slepc.upv.es/download/distrib/slepc-2.1.5.tgz)
[[2.2.0]](https://slepc.upv.es/download/distrib/slepc-2.2.0.tgz)
[[2.2.1]](https://slepc.upv.es/download/distrib/slepc-2.2.1.tgz)
[[2.3.0]](https://slepc.upv.es/download/distrib/slepc-2.3.0.tgz)
[[2.3.1]](https://slepc.upv.es/download/distrib/slepc-2.3.1.tgz)
[[2.3.2]](https://slepc.upv.es/download/distrib/slepc-2.3.2.tgz)
[[2.3.3]](https://slepc.upv.es/download/distrib/slepc-2.3.3.tgz)

Patches:
[[2.1.5]](https://slepc.upv.es/download/distrib/patches/slepc_patch_all-2.1.5)
[[2.2.0]](https://slepc.upv.es/download/distrib/patches/slepc_patch_all-2.2.0)
[[2.2.1]](https://slepc.upv.es/download/distrib/patches/slepc_patch_all-2.2.1)
[[2.3.0]](https://slepc.upv.es/download/distrib/patches/slepc_patch_all-2.3.0)
[[2.3.1]](https://slepc.upv.es/download/distrib/patches/slepc_patch_all-2.3.1)
[[2.3.2]](https://slepc.upv.es/download/distrib/patches/slepc_patch_all-2.3.2)
[[2.3.3]](https://slepc.upv.es/download/distrib/patches/slepc_patch_all-2.3.3)

For 2.3.3 and older releases, apply the patch with GNU patch command in the SLEPc root directory as `patch -Np1 < patchfile`.
:::

## Changes

A change log of the different versions can be found in [CHANGELOG.md](https://gitlab.com/slepc/slepc/-/blob/main/CHANGELOG.md). It includes the release date of each version. For easy comparison, the following table lists all versions of SLEPc, showing the correspondence between SLEPc and PETSc releases, as well as the release date.

:::{admonition} Table of SLEPc/PETSc Releases
:class: dropdown

SLEPc version | PETSc versions      | Release date
---           | ---                 | ---
2.1.0         | 2.1.0               | Not released
2.1.1         | 2.1.1, 2.1.2, 2.1.3 | Dec 2002
2.1.5         | 2.1.5, 2.1.6        | May 2003
2.2.0         | 2.2.0               | Apr 2004
2.2.1         | 2.2.1               | Aug 2004
2.3.0         | 2.3.0               | Jun 2005
2.3.1         | 2.3.1               | Mar 2006
2.3.2         | 2.3.1, 2.3.2        | Oct 2006
2.3.3         | 2.3.3               | Jun 2007
3.0.0         | 3.0.0               | Feb 2009
3.1           | 3.1                 | Aug 2010
3.2           | 3.2                 | Oct 2011
3.3           | 3.3                 | Aug 2012
3.4           | 3.4                 | Jul 2013
3.5           | 3.5                 | Jul 2014
3.6           | 3.6                 | Jun 2015
3.7           | 3.7                 | May 2016
3.8           | 3.8                 | Oct 2017
3.9           | 3.9                 | Apr 2018
3.10          | 3.10                | Sep 2018
3.11          | 3.11                | Mar 2019
3.12          | 3.12                | Sep 2019
3.13          | 3.13                | Mar 2020
3.14          | 3.14                | Sep 2020
3.15          | 3.15                | Mar 2021
3.16          | 3.16                | Sep 2021
3.17          | 3.17                | Mar 2022
3.18          | 3.18                | Oct 2022
3.19          | 3.19                | Mar 2023
3.20          | 3.20                | Sep 2023
3.21          | 3.21                | Mar 2024
3.22          | 3.22                | Sep 2024
3.23          | 3.23                | Mar 2025
3.24          | 3.24                | Sep 2025

:::
