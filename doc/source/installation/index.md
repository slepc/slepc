(install)=

# Installation

SLEPc is available on some package managers. Select the tabs below to show the corresponding information. Note that some of them provide separate packages for the real and complex versions.

::::{tab-set}

:::{tab-item} Archlinux
<https://aur.archlinux.org/packages/slepc>
:::

:::{tab-item} Conda
<https://anaconda.org/conda-forge/slepc>

    $ conda install conda-forge::slepc
:::

:::{tab-item} Debian
<https://packages.debian.org/search?keywords=slepc-dev>

    $ sudo apt install slepc-dev
:::

:::{tab-item} E4S
<https://e4s.io/download.html>

Extreme-scale Scientific Software Stack (E4S) provides container images for docker and singularity.

E4S packages are available as pre-built Spack binaries.
:::

:::{tab-item} Homebrew
<https://formulae.brew.sh/formula/slepc>

    $ brew install slepc
:::

:::{tab-item} MacPorts
<https://ports.macports.org/port/slepc/>

    $ sudo port install slepc
:::

:::{tab-item} MSYS2 (Windows)
<https://packages.msys2.org/packages/mingw-w64-x86_64-slepc>

    $ pacman -S mingw-w64-x86_64-slepc
:::

:::{tab-item} openSUSE
<https://software.opensuse.org/package/slepc>

    $ sudo zypper install slepc
:::

:::{tab-item} python
<https://pypi.org/project/slepc/>

**release**

    $ python -m pip install numpy mpi4py
    $ python -m pip install petsc petsc4py
    $ python -m pip install slepc slepc4py

**development**

    $ python -m pip install Cython numpy mpi4py
    $ python -m pip install --no-deps https://gitlab.com/petsc/petsc/-/archive/main/petsc-main.tar.gz
    $ python -m pip install --no-deps https://gitlab.com/slepc/slepc/-/archive/main/slepc-main.tar.gz
:::

:::{tab-item} Spack

**debug install**

    $ spack install slepc +debug

**list available variants (configurations)**

    $ spack info
:::

:::{tab-item} Ubuntu
 <https://packages.ubuntu.com/search?keywords=slepc-dev>

    $ sudo apt install slepc-dev
:::

::::

```{raw} html
<br>
```

Alternatively, you can download and install it from source code, according to the following instructions:

```{toctree}
:maxdepth: 1

download
quickstart
```
