cdef extern from * nogil:

    ctypedef char* SlepcDSType "const char*"
    SlepcDSType DSHEP
    SlepcDSType DSNHEP
    SlepcDSType DSGHEP
    SlepcDSType DSGHIEP
    SlepcDSType DSGNHEP
    SlepcDSType DSNHEPTS
    SlepcDSType DSSVD
    SlepcDSType DSHSVD
    SlepcDSType DSGSVD
    SlepcDSType DSPEP
    SlepcDSType DSNEP

    ctypedef enum SlepcDSStateType "DSStateType":
        DS_STATE_RAW
        DS_STATE_INTERMEDIATE
        DS_STATE_CONDENSED
        DS_STATE_TRUNCATED

    ctypedef enum SlepcDSMatType "DSMatType":
        DS_MAT_A
        DS_MAT_B
        DS_MAT_C
        DS_MAT_T
        DS_MAT_D
        DS_MAT_Q
        DS_MAT_Z
        DS_MAT_X
        DS_MAT_Y
        DS_MAT_U
        DS_MAT_V
        DS_MAT_W
        DS_NUM_MAT

    ctypedef enum SlepcDSParallelType "DSParallelType":
        DS_PARALLEL_REDUNDANT
        DS_PARALLEL_SYNCHRONIZED
        DS_PARALLEL_DISTRIBUTED

    PetscErrorCode DSCreate(MPI_Comm,SlepcDS*)
    PetscErrorCode DSView(SlepcDS,PetscViewer)
    PetscErrorCode DSDestroy(SlepcDS*)
    PetscErrorCode DSReset(SlepcDS)
    PetscErrorCode DSSetType(SlepcDS,SlepcDSType)
    PetscErrorCode DSGetType(SlepcDS,SlepcDSType*)

    PetscErrorCode DSSetOptionsPrefix(SlepcDS,char[])
    PetscErrorCode DSGetOptionsPrefix(SlepcDS,char*[])
    PetscErrorCode DSAppendOptionsPrefix(SlepcDS,char[])
    PetscErrorCode DSSetFromOptions(SlepcDS)
    PetscErrorCode DSDuplicate(SlepcDS,SlepcDS*)

    PetscErrorCode DSAllocate(SlepcDS,PetscInt)
    PetscErrorCode DSGetLeadingDimension(SlepcDS,PetscInt*)
    PetscErrorCode DSSetState(SlepcDS,SlepcDSStateType)
    PetscErrorCode DSGetState(SlepcDS,SlepcDSStateType*)
    PetscErrorCode DSSetDimensions(SlepcDS,PetscInt,PetscInt,PetscInt)
    PetscErrorCode DSGetDimensions(SlepcDS,PetscInt*,PetscInt*,PetscInt*,PetscInt*)
    PetscErrorCode DSTruncate(SlepcDS,PetscInt,PetscBool)
    PetscErrorCode DSSetBlockSize(SlepcDS,PetscInt)
    PetscErrorCode DSGetBlockSize(SlepcDS,PetscInt*)
    PetscErrorCode DSSetMethod(SlepcDS,PetscInt)
    PetscErrorCode DSGetMethod(SlepcDS,PetscInt*)
    PetscErrorCode DSSetParallel(SlepcDS,SlepcDSParallelType)
    PetscErrorCode DSGetParallel(SlepcDS,SlepcDSParallelType*)
    PetscErrorCode DSSetCompact(SlepcDS,PetscBool)
    PetscErrorCode DSGetCompact(SlepcDS,PetscBool*)
    PetscErrorCode DSSetExtraRow(SlepcDS,PetscBool)
    PetscErrorCode DSGetExtraRow(SlepcDS,PetscBool*)
    PetscErrorCode DSSetRefined(SlepcDS,PetscBool)
    PetscErrorCode DSGetRefined(SlepcDS,PetscBool*)
    PetscErrorCode DSGetMat(SlepcDS,SlepcDSMatType,PetscMat*)
    PetscErrorCode DSRestoreMat(SlepcDS,SlepcDSMatType,PetscMat*)
    PetscErrorCode DSGetArray(SlepcDS,SlepcDSMatType,PetscScalar**)
    PetscErrorCode DSRestoreArray(SlepcDS,SlepcDSMatType,PetscScalar**)
    PetscErrorCode DSMatGetSize(SlepcDS,SlepcDSMatType,PetscInt*,PetscInt*)
    PetscErrorCode DSSetIdentity(SlepcDS,SlepcDSMatType)
    PetscErrorCode DSVectors(SlepcDS,SlepcDSMatType,PetscInt*,PetscReal*)
    PetscErrorCode DSSolve(SlepcDS,PetscScalar*,PetscScalar*)
    PetscErrorCode DSSort(SlepcDS,PetscScalar*,PetscScalar*,PetscScalar*,PetscScalar*,PetscInt*)
    PetscErrorCode DSUpdateExtraRow(SlepcDS)
    PetscErrorCode DSCond(SlepcDS,PetscReal*)
    PetscErrorCode DSTranslateHarmonic(SlepcDS,PetscScalar,PetscReal,PetscBool,PetscScalar*,PetscReal*)
    PetscErrorCode DSTranslateRKS(SlepcDS,PetscScalar)
    PetscErrorCode DSCond(SlepcDS,PetscReal*)
    PetscErrorCode DSNormalize(SlepcDS,SlepcDSMatType,PetscInt)

    PetscErrorCode DSSVDSetDimensions(SlepcDS,PetscInt)
    PetscErrorCode DSSVDGetDimensions(SlepcDS,PetscInt*)
    PetscErrorCode DSHSVDSetDimensions(SlepcDS,PetscInt)
    PetscErrorCode DSHSVDGetDimensions(SlepcDS,PetscInt*)
    PetscErrorCode DSGSVDSetDimensions(SlepcDS,PetscInt,PetscInt)
    PetscErrorCode DSGSVDGetDimensions(SlepcDS,PetscInt*,PetscInt*)

    PetscErrorCode DSPEPSetDegree(SlepcDS,PetscInt)
    PetscErrorCode DSPEPGetDegree(SlepcDS,PetscInt*)
    PetscErrorCode DSPEPSetCoefficients(SlepcDS,PetscReal*)
    PetscErrorCode DSPEPGetCoefficients(SlepcDS,PetscReal**)

