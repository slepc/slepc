# --------------------------------------------------------------------

cdef extern from * nogil:
    """
    #include "pyapicompat.h"
    #include "lib-slepc/compat.h"
    #include "lib-slepc/custom.h"

    /* Silence Clang warnings in Cython-generated C code */
    #if defined(__clang__)
      #pragma clang diagnostic ignored "-Wextra-semi-stmt"
      #pragma clang diagnostic ignored "-Wparentheses-equality"
      #pragma clang diagnostic ignored "-Wunreachable-code-fallthrough"
      #pragma clang diagnostic ignored "-Woverlength-strings"
      #pragma clang diagnostic ignored "-Wunreachable-code"
      #pragma clang diagnostic ignored "-Wundef"
    #elif defined(__GNUC__) || defined(__GNUG__)
      #pragma GCC diagnostic ignored "-Wstrict-aliasing"
      #pragma GCC diagnostic ignored "-Wtype-limits"
    #endif
    """

# -----------------------------------------------------------------------------

from petsc4py.PETSc import COMM_NULL
from petsc4py.PETSc import COMM_SELF
from petsc4py.PETSc import COMM_WORLD

# -----------------------------------------------------------------------------

from petsc4py.PETSc cimport MPI_Comm
from petsc4py.PETSc cimport PetscErrorCode, PetscErrorType
from petsc4py.PETSc cimport PETSC_SUCCESS, PETSC_ERR_PYTHON
from petsc4py.PETSc cimport CHKERR
from petsc4py.PETSc cimport PetscObject, PetscViewer
from petsc4py.PETSc cimport PetscRandom
from petsc4py.PETSc cimport PetscVec, PetscMat
from petsc4py.PETSc cimport PetscKSP, PetscPC

from petsc4py.PETSc cimport Comm
from petsc4py.PETSc cimport Object, Viewer
from petsc4py.PETSc cimport Random
from petsc4py.PETSc cimport Vec, Mat
from petsc4py.PETSc cimport KSP, PC

# -----------------------------------------------------------------------------

cdef inline object bytes2str(const char p[]):
     if p == NULL:
         return None
     cdef bytes s = <char*>p
     if isinstance(s, str):
         return s
     else:
         return s.decode()

cdef inline object str2bytes(object s, const char *p[]):
    if s is None:
        p[0] = NULL
        return None
    if not isinstance(s, bytes):
        s = s.encode()
    p[0] = <const char*>(<char*>s)
    return s

cdef inline object S_(const char p[]):
     if p == NULL: return None
     cdef object s = <char*>p
     return s if isinstance(s, str) else s.decode()

include "allocate.pxi"

# -----------------------------------------------------------------------------

cdef extern from * nogil:
    ctypedef bint      PetscBool
    const    PetscBool PETSC_TRUE
    const    PetscBool PETSC_FALSE
    ctypedef long      PetscInt
    ctypedef double    PetscReal
    ctypedef double    PetscScalar

cdef inline object toBool(PetscBool value):
    return True if value == PETSC_TRUE else False
cdef inline PetscBool asBool(object value) except? <PetscBool>0:
    return PETSC_TRUE if value else PETSC_FALSE

cdef inline object toInt(PetscInt value):
    return value
cdef inline PetscInt asInt(object value) except? -1:
    return value

cdef inline object toReal(PetscReal value):
    return value
cdef inline PetscReal asReal(object value) except? -1:
    return value

cdef extern from "<petsc4py/pyscalar.h>":
    object      PyPetscScalar_FromPetscScalar(PetscScalar)
    PetscScalar PyPetscScalar_AsPetscScalar(object) except? <PetscScalar>-1.0

cdef inline object toScalar(PetscScalar value):
    return PyPetscScalar_FromPetscScalar(value)
cdef inline PetscScalar asScalar(object value) except? <PetscScalar>-1.0:
    return PyPetscScalar_AsPetscScalar(value)

cdef extern from "Python.h":
     double PyComplex_RealAsDouble(object)
     double PyComplex_ImagAsDouble(object)

cdef inline object toComplex(PetscScalar rvalue, PetscScalar ivalue):
    return complex(toScalar(rvalue), toScalar(ivalue))

cdef inline PetscReal asComplexReal(object value) except? <PetscReal>-1.0:
    return <PetscReal>PyComplex_RealAsDouble(value)

cdef inline PetscReal asComplexImag(object value) except? <PetscReal>-1.0:
    return <PetscReal>PyComplex_ImagAsDouble(value)

cdef extern from * nogil:
    PetscReal PetscRealPart(PetscScalar v)
    PetscReal PetscImaginaryPart(PetscScalar v)

# --------------------------------------------------------------------

# NumPy support
# -------------

include "arraynpy.pxi"

import_array()

IntType     = PyArray_TypeObjectFromType(NPY_PETSC_INT)
RealType    = PyArray_TypeObjectFromType(NPY_PETSC_REAL)
ScalarType  = PyArray_TypeObjectFromType(NPY_PETSC_SCALAR)
ComplexType = PyArray_TypeObjectFromType(NPY_PETSC_COMPLEX)

# -----------------------------------------------------------------------------

cdef extern from "<string.h>"  nogil:
    void* memset(void*,int,size_t)
    void* memcpy(void*,void*,size_t)
    char* strdup(char*)

# -----------------------------------------------------------------------------

include "slepcmpi.pxi"
include "slepcsys.pxi"
include "slepcutil.pxi"
include "slepcst.pxi"
include "slepcbv.pxi"
include "slepcds.pxi"
include "slepcfn.pxi"
include "slepcrg.pxi"
include "slepceps.pxi"
include "slepcsvd.pxi"
include "slepcpep.pxi"
include "slepcnep.pxi"
include "slepcmfn.pxi"
include "slepclme.pxi"

# -----------------------------------------------------------------------------

__doc__ = u"""
Scalable Library for Eigenvalue Problem Computations
"""

DECIDE    = PETSC_DECIDE
DEFAULT   = PETSC_DEFAULT
DETERMINE = PETSC_DETERMINE
CURRENT   = PETSC_CURRENT

include "Sys.pyx"
include "Util.pyx"
include "ST.pyx"
include "BV.pyx"
include "DS.pyx"
include "FN.pyx"
include "RG.pyx"
include "EPS.pyx"
include "SVD.pyx"
include "PEP.pyx"
include "NEP.pyx"
include "MFN.pyx"
include "LME.pyx"

# -----------------------------------------------------------------------------

include "CAPI.pyx"

# -----------------------------------------------------------------------------

cdef extern from "Python.h":
    int Py_AtExit(void (*)() noexcept nogil)
    void PySys_WriteStderr(char*,...)

cdef extern from "<stdio.h>" nogil:
    ctypedef struct FILE
    FILE *stderr
    int fprintf(FILE *, char *, ...)

cdef int initialize(object args) except PETSC_ERR_PYTHON:
    if (<int>SlepcInitializeCalled): return 1
    if (<int>SlepcFinalizeCalled):   return 0
    # initialize SLEPC
    CHKERR( SlepcInitialize(NULL, NULL, NULL, NULL) )
    # register finalization function
    if Py_AtExit(finalize) < 0:
        PySys_WriteStderr(b"warning: could not register %s with Py_AtExit()",
                          b"SlepcFinalize()")
    return 1 # and we are done, enjoy !!

from petsc4py.PETSc cimport PyPetscType_Register

cdef extern from * nogil:
    PetscErrorCode SlepcInitializePackageAll()
    ctypedef int PetscClassId
    PetscClassId SLEPC_ST_CLASSID  "ST_CLASSID"
    PetscClassId SLEPC_BV_CLASSID  "BV_CLASSID"
    PetscClassId SLEPC_DS_CLASSID  "DS_CLASSID"
    PetscClassId SLEPC_FN_CLASSID  "FN_CLASSID"
    PetscClassId SLEPC_RG_CLASSID  "RG_CLASSID"
    PetscClassId SLEPC_EPS_CLASSID "EPS_CLASSID"
    PetscClassId SLEPC_SVD_CLASSID "SVD_CLASSID"
    PetscClassId SLEPC_PEP_CLASSID "PEP_CLASSID"
    PetscClassId SLEPC_NEP_CLASSID "NEP_CLASSID"
    PetscClassId SLEPC_MFN_CLASSID "MFN_CLASSID"
    PetscClassId SLEPC_LME_CLASSID "LME_CLASSID"

cdef PetscErrorCode register() except PETSC_ERR_PYTHON:
    # make sure all SLEPc packages are initialized
    CHKERR( SlepcInitializePackageAll() )
    # register Python types
    PyPetscType_Register(SLEPC_ST_CLASSID,  ST)
    PyPetscType_Register(SLEPC_BV_CLASSID,  BV)
    PyPetscType_Register(SLEPC_DS_CLASSID,  DS)
    PyPetscType_Register(SLEPC_FN_CLASSID,  FN)
    PyPetscType_Register(SLEPC_RG_CLASSID,  RG)
    PyPetscType_Register(SLEPC_EPS_CLASSID, EPS)
    PyPetscType_Register(SLEPC_SVD_CLASSID, SVD)
    PyPetscType_Register(SLEPC_PEP_CLASSID, PEP)
    PyPetscType_Register(SLEPC_NEP_CLASSID, NEP)
    PyPetscType_Register(SLEPC_MFN_CLASSID, MFN)
    PyPetscType_Register(SLEPC_LME_CLASSID, LME)
    return PETSC_SUCCESS

cdef void finalize() noexcept nogil:
    cdef PetscErrorCode ierr = PETSC_SUCCESS
    # manage SLEPc finalization
    if not (<int>SlepcInitializeCalled): return
    if (<int>SlepcFinalizeCalled): return
    # finalize SLEPc
    ierr = SlepcFinalize()
    if ierr != PETSC_SUCCESS:
        fprintf(stderr, "SlepcFinalize() failed "
                "[error code: %d]\n", <int>ierr)
    # and we are done, see you later !!

# -----------------------------------------------------------------------------

def _initialize(args=None):
    cdef int ready = initialize(args)
    if ready: register()

def _finalize():
    finalize()

# -----------------------------------------------------------------------------
