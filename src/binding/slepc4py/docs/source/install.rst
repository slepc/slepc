Installation
============

Using **pip** or **easy_install**
---------------------------------

You can use :program:`pip` to install :mod:`slepc4py` and its
dependencies (:mod:`mpi4py` is optional but highly recommended)::

  $ pip install [--user] numpy mpi4py
  $ pip install [--user] petsc petsc4py
  $ pip install [--user] slepc slepc4py

Alternatively, you can use :program:`easy_install` (deprecated)::

  $ easy_install [--user] slepc4py

If you already have working PETSc and SLEPc installs, set environment
variables :envvar:`SLEPC_DIR` and :envvar:`PETSC_DIR` (and perhaps
:envvar:`PETSC_ARCH` for non-prefix installs) to appropriate values
and next use :program:`pip`::

  $ export SLEPC_DIR=/path/to/slepc
  $ export PETSC_DIR=/path/to/petsc
  $ export PETSC_ARCH=arch-linux2-c-opt
  $ pip install [--user] petsc4py slepc4py


Using **distutils**
-------------------

Requirements
^^^^^^^^^^^^

You need to have the following software properly installed in order to
build *SLEPc for Python*:

* Any MPI_ implementation [#]_ (e.g., MPICH_ or `Open MPI`_),
  built with shared libraries.

* A matching version of PETSc_ built with shared libraries.

* A matching version of SLEPc_ built with shared libraries.

* NumPy_ package.

* petsc4py_ package.

.. [#] Unless you have appropriately configured and built SLEPc and
       PETSc without MPI (configure option ``--with-mpi=0``).

.. [#] You may need to use a parallelized version of the Python
       interpreter with some MPI-1 implementations (e.g. MPICH1).

.. include:: links.txt

Downloading
^^^^^^^^^^^

The *SLEPc for Python* package is available for download at the
Python Package Index. You can use
:program:`curl` or :program:`wget` to get a release tarball.

* Using :program:`curl`::

    $ curl -LO https://pypi.io/packages/source/s/slepc4py/slepc4py-X.Y.Z.tar.gz

* Using :program:`wget`::

    $ wget https://pypi.io/packages/source/s/slepc4py/slepc4py-X.Y.Z.tar.gz

Building
^^^^^^^^

After unpacking the release tarball::

  $ tar -zxf slepc4py-X.Y.tar.gz
  $ cd slepc4py-X.Y

the distribution is ready for building.

.. note:: **Mac OS X** users employing a Python distribution built
   with **universal binaries** may need to set the environment
   variables :envvar:`MACOSX_DEPLOYMENT_TARGET`, :envvar:`SDKROOT`,
   and :envvar:`ARCHFLAGS` to appropriate values. As an example,
   assume your Mac is running **Snow Leopard** on a **64-bit Intel**
   processor and you want to override the hard-wired cross-development
   SDK in Python configuration, your environment should be modified
   like this::

     $ export MACOSX_DEPLOYMENT_TARGET=10.6
     $ export SDKROOT=/
     $ export ARCHFLAGS='-arch x86_64'

Some environment configuration is needed to inform the location of
PETSc and SLEPc. You can set (using :command:`setenv`,
:command:`export` or what applies to you shell or system) the
environment variables :envvar:`SLEPC_DIR`, :envvar:`PETSC_DIR`, and
:envvar:`PETSC_ARCH` indicating where you have built/installed SLEPc
and PETSc::

  $ export SLEPC_DIR=/usr/local/slepc
  $ export PETSC_DIR=/usr/local/petsc
  $ export PETSC_ARCH=arch-linux2-c-opt

Alternatively, you can edit the file :file:`setup.cfg` and provide the
required information below the ``[config]`` section::

  [config]
  slepc_dir  = /usr/local/slepc
  petsc_dir  = /usr/local/petsc
  petsc_arch = arch-linux2-c-opt
  ...

Finally, you can build the distribution by typing::

  $ python setup.py build

Installing
^^^^^^^^^^

After building, the distribution is ready for installation.

If you have root privileges (either by log-in as the root user of by
using :command:`sudo`) and you want to install *SLEPc for Python* in
your system for all users, just do::

  $ python setup.py install

The previous steps will install the :mod:`slepc4py` package at standard
location :file:`{prefix}/lib/python{X}.{X}/site-packages`.

If you do not have root privileges or you want to install *SLEPc for
Python* for your private use, just do::

  $ python setup.py install --user
