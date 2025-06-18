# -----------------------------------------------------------------------------

"""
Extension modules for different SLEPc configurations.

SLEPc can be configured with different options (eg. debug/optimized,
single/double precisionm, C/C++ compilers, external packages). Each
configuration variant is associated to a name, frequently available as
an environmental variable named ``PETSC_ARCH``.

This package is a holds all the available variants of the SLEPc
extension module built against specific SLEPc configurations. It also
provides a convenience function using of the builtin ``imp`` module
for easily importing any of the available extension modules depending
on the value of a user-provided configuration name, the ``PETSC_ARCH``
environmental variable, or a configuration file.
"""

# -----------------------------------------------------------------------------

from petsc4py.lib import ImportPETSc
from petsc4py.lib import Import, getPathArch, getInitArgs  # noqa: F401


def ImportSLEPc(arch=None):
    """
    Import the SLEPc extension module for a given configuration name.
    """
    path, arch = getPathArchSLEPc(arch)
    PETSc = ImportPETSc(arch)  # noqa: F841
    return Import('slepc4py', 'SLEPc', path, arch)


def getPathArchSLEPc(arch=None):
    """
    Undocumented.
    """
    import os
    import sys
    PETSc = sys.modules.get('petsc4py.PETSc')
    arch = getattr(PETSc, '__arch__', arch)
    path = os.path.dirname(__file__)
    rcvar, rcfile  =  'PETSC_ARCH', 'slepc.cfg'
    path, arch = getPathArch(path, arch, rcvar, rcfile)
    return (path, arch)

# -----------------------------------------------------------------------------
