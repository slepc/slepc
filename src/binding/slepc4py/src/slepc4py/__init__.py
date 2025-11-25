# -----------------------------------------------------------------------------

"""
SLEPc for Python
================

This package is an interface to SLEPc_ libraries.

SLEPc_ (the Scalable Library for Eigenvalue Problem Computations) is a
software library for the solution of large scale sparse eigenvalue
problems on parallel computers. It is an extension of PETSc_ and can
be used for either standard or generalized eigenproblems, with real or
complex arithmetic. It can also be used for computing a partial SVD of
a large, sparse, rectangular matrix, and to solve nonlinear eigenvalue
problems (polynomial or general). Additionally, SLEPc provides solvers
for the computation of the action of a matrix function on a vector.

.. _SLEPc: https://slepc.upv.es
.. _PETSc: https://petsc.org
"""

__author__  = 'Lisandro Dalcin'
__version__ = '3.24.1'
__credits__ = 'SLEPc Team <slepc-maint@upv.es>'

# -----------------------------------------------------------------------------

def init(args=None, arch=None, comm=None):
    """
    Initialize SLEPc.

    Parameters
    ----------
    args
        Command-line arguments, usually the `sys.argv` list.
    arch
        Specific configuration to use.
    comm
        MPI commmunicator.

    Notes
    -----
    This function should be called only once, typically at the very
    beginning of the bootstrap script of an application.
    """
    import slepc4py.lib
    SLEPc = slepc4py.lib.ImportSLEPc(arch)
    PETSc = slepc4py.lib.ImportPETSc(arch)
    args  = slepc4py.lib.getInitArgs(args)
    PETSc._initialize(args, comm)
    SLEPc._initialize(args)

# -----------------------------------------------------------------------------

def get_include():
    """
    Return the directory in the package that contains header files.

    Extension modules that need to compile against slepc4py should use
    this function to locate the appropriate include directory.

    Example
    -------
    Using Python distutils or NumPy distutils::

      import petscc4py, slepc4py
      Extension('extension_name', ...
                include_dirs=[...,
                              petsc4py.get_include(),
                              slepc4py.get_include(),])
    """
    from os.path import dirname, join
    return join(dirname(__file__), 'include')

# -----------------------------------------------------------------------------

def get_config():
    """Return a dictionary with information about SLEPc."""
    import os

    from io import StringIO
    from configparser import ConfigParser
    pgkdir = os.path.dirname(__file__)
    filename = os.path.join(pgkdir, 'lib', 'slepc.cfg')
    with open(filename) as fp:
        stream = StringIO('[slepc]\n' + fp.read())
    parser = ConfigParser()
    parser.optionxform = str
    parser.read_file(stream, filename)
    return dict(parser.items('slepc'))
