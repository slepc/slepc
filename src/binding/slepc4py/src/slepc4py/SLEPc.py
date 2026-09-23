"""The SLEPc module."""

ARCH = None
from slepc4py.lib import ImportSLEPc  # noqa: E402
from slepc4py.lib import ImportPETSc  # noqa: E402

SLEPc = ImportSLEPc(ARCH)
PETSc = ImportPETSc(ARCH)
PETSc._initialize()
SLEPc._initialize()
del SLEPc, PETSc
del ImportSLEPc, ImportPETSc
del ARCH
