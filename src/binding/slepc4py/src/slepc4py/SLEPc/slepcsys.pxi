cdef extern from * nogil:

    enum: PETSC_DECIDE
    enum: PETSC_DEFAULT
    enum: PETSC_DETERMINE
    enum: PETSC_CURRENT

    ctypedef const char* PetscVecType "VecType"

    ctypedef enum  PetscNormType "NormType":
        PETSC_NORM_1          "NORM_1"
        PETSC_NORM_2          "NORM_2"
        PETSC_NORM_1_AND_2    "NORM_1_AND_2"
        PETSC_NORM_FROBENIUS  "NORM_FROBENIUS"
        PETSC_NORM_INFINITY   "NORM_INFINITY"
        PETSC_NORM_MAX        "NORM_MAX"

    ctypedef enum  PetscMatStructure "MatStructure":
        MAT_SAME_NONZERO_PATTERN      "SAME_NONZERO_PATTERN"
        MAT_DIFFERENT_NONZERO_PATTERN "DIFFERENT_NONZERO_PATTERN"
        MAT_SUBSET_NONZERO_PATTERN    "SUBSET_NONZERO_PATTERN"

cdef extern from * nogil:
    PetscErrorCode PetscMalloc(size_t,void*)
    PetscErrorCode PetscFree(void*)
    PetscErrorCode PetscMemcpy(void*,void*,size_t)
    PetscErrorCode PetscMemzero(void*,size_t)

cdef extern from * nogil:
    MPI_Comm PetscObjectComm(PetscObject)
    PetscErrorCode PetscObjectReference(PetscObject)
    PetscErrorCode PetscObjectDereference(PetscObject)
    PetscErrorCode PetscObjectDestroy(PetscObject*)
    PetscErrorCode PetscObjectTypeCompare(PetscObject,char[],PetscBool*)

cdef extern from * nogil:
    PetscErrorCode VecCopy(PetscVec,PetscVec)
    PetscErrorCode VecSet(PetscVec,PetscScalar)
    PetscErrorCode VecDestroy(PetscVec*)

cdef extern from * nogil:
    PetscErrorCode MatGetSize(PetscMat,PetscInt*,PetscInt*)
    PetscErrorCode MatGetLocalSize(PetscMat,PetscInt*,PetscInt*)

cdef extern from * nogil:
    const char SLEPC_AUTHOR_INFO[]
    PetscErrorCode SlepcGetVersion(char[],size_t)
    PetscErrorCode SlepcGetVersionNumber(PetscInt*,PetscInt*,PetscInt*,PetscInt*)

    PetscErrorCode SlepcInitialize(int*,char***,char[],char[])
    PetscErrorCode SlepcFinalize()
    PetscBool SlepcInitializeCalled
    PetscBool SlepcFinalizeCalled

    PetscErrorCode SlepcHasExternalPackage(const char[],PetscBool*)

cdef inline PetscMatStructure matstructure(object structure) \
    except <PetscMatStructure>(-1):
    if   structure is None:  return MAT_DIFFERENT_NONZERO_PATTERN
    elif structure is False: return MAT_DIFFERENT_NONZERO_PATTERN
    elif structure is True:  return MAT_SAME_NONZERO_PATTERN
    else:                    return structure

cdef inline PetscErrorCode PetscINCREF(PetscObject *obj):
    if obj    == NULL: return PETSC_SUCCESS
    if obj[0] == NULL: return PETSC_SUCCESS
    return PetscObjectReference(obj[0])

cdef inline PetscErrorCode SlepcCLEAR(PetscObject* obj):
    if obj    == NULL: return PETSC_SUCCESS
    if obj[0] == NULL: return PETSC_SUCCESS
    cdef PetscObject tmp
    tmp = obj[0]; obj[0] = NULL
    return PetscObjectDestroy(&tmp)

cdef inline PetscViewer def_Viewer(Viewer viewer):
   return viewer.vwr if viewer is not None else <PetscViewer>NULL

cdef inline KSP ref_KSP(PetscKSP ksp):
    cdef KSP ob = <KSP> KSP()
    ob.ksp = ksp
    CHKERR( PetscINCREF(ob.obj) )
    return ob

cdef inline Mat ref_Mat(PetscMat mat):
    cdef Mat ob = <Mat> Mat()
    ob.mat = mat
    CHKERR( PetscINCREF(ob.obj) )
    return ob

cdef inline Vec ref_Vec(PetscVec vec):
    cdef Vec ob = <Vec> Vec()
    ob.vec = vec
    CHKERR( PetscINCREF(ob.obj) )
    return ob

