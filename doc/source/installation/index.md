# Installation

SLEPc is available on some package managers.

::::{tab-set}

:::{tab-item} Archlinux
<https://aur.archlinux.org/packages/slepc>
:::

:::{tab-item} Conda
<https://anaconda.org/conda-forge/slepc>

    $ conda install -c conda-forge slepc
:::

:::{tab-item} Debian
<https://packages.debian.org/slepc-dev>

    $ sudo apt install slepc-dev
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
 <https://packages.ubuntu.com/slepc-dev>

    $ sudo apt install slepc-dev
:::

::::

And as a direct [download](download)

```{toctree}
:maxdepth: 2
:hidden:

download
older
quickstart
```
