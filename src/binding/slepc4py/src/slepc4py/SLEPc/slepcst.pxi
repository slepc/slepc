cdef extern from * nogil:

    ctypedef char* SlepcSTType "const char*"
    SlepcSTType STSHELL
    SlepcSTType STSHIFT
    SlepcSTType STSINVERT
    SlepcSTType STCAYLEY
    SlepcSTType STPRECOND
    SlepcSTType STFILTER

    ctypedef enum SlepcSTMatMode "STMatMode":
        ST_MATMODE_COPY
        ST_MATMODE_INPLACE
        ST_MATMODE_SHELL

    PetscErrorCode STView(SlepcST,PetscViewer)
    PetscErrorCode STDestroy(SlepcST*)
    PetscErrorCode STReset(SlepcST)
    PetscErrorCode STCreate(MPI_Comm,SlepcST*)
    PetscErrorCode STGetType(SlepcST,SlepcSTType*)
    PetscErrorCode STSetType(SlepcST,SlepcSTType)
    PetscErrorCode STGetOptionsPrefix(SlepcST,char*[])
    PetscErrorCode STSetOptionsPrefix(SlepcST,char[])
    PetscErrorCode STAppendOptionsPrefix(SlepcST,char[])
    PetscErrorCode STSetFromOptions(SlepcST)

    PetscErrorCode STGetShift(SlepcST,PetscScalar*)
    PetscErrorCode STSetShift(SlepcST,PetscScalar)

    PetscErrorCode STGetKSP(SlepcST,PetscKSP*)
    PetscErrorCode STSetKSP(SlepcST,PetscKSP)

    PetscErrorCode STGetOperator(SlepcST,PetscMat*)
    PetscErrorCode STRestoreOperator(SlepcST,PetscMat*)

    PetscErrorCode STGetNumMatrices(SlepcST,PetscInt*)
    PetscErrorCode STGetMatrix(SlepcST,PetscInt,PetscMat*)
    PetscErrorCode STSetMatrices(SlepcST,PetscInt,PetscMat*)
    PetscErrorCode STSetMatStructure(SlepcST,PetscMatStructure)
    PetscErrorCode STGetMatStructure(SlepcST,PetscMatStructure*)
    PetscErrorCode STSetPreconditionerMat(SlepcST,PetscMat)
    PetscErrorCode STGetPreconditionerMat(SlepcST,PetscMat*)
    PetscErrorCode STSetSplitPreconditioner(SlepcST,PetscInt,PetscMat*,PetscMatStructure)
    PetscErrorCode STGetSplitPreconditionerTerm(SlepcST,PetscInt,PetscMat*)
    PetscErrorCode STGetSplitPreconditionerInfo(SlepcST,PetscInt*,PetscMatStructure*)

    PetscErrorCode STSetTransform(SlepcST,PetscBool)
    PetscErrorCode STGetTransform(SlepcST,PetscBool*)

    PetscErrorCode STGetMatMode(SlepcST,SlepcSTMatMode*)
    PetscErrorCode STSetMatMode(SlepcST,SlepcSTMatMode)

    PetscErrorCode STSetUp(SlepcST)
    PetscErrorCode STApply(SlepcST,PetscVec,PetscVec)
    PetscErrorCode STApplyMat(SlepcST,PetscMat,PetscMat)
    PetscErrorCode STApplyTranspose(SlepcST,PetscVec,PetscVec)
    PetscErrorCode STApplyHermitianTranspose(SlepcST,PetscVec,PetscVec)

    PetscErrorCode STCayleySetAntishift(SlepcST,PetscScalar)
    PetscErrorCode STCayleyGetAntishift(SlepcST,PetscScalar*)

    ctypedef enum SlepcSTFilterType "STFilterType":
        ST_FILTER_FILTLAN
        ST_FILTER_CHEBYSHEV

    ctypedef enum SlepcSTFilterDamping "STFilterDamping":
        ST_FILTER_DAMPING_NONE
        ST_FILTER_DAMPING_JACKSON
        ST_FILTER_DAMPING_LANCZOS
        ST_FILTER_DAMPING_FEJER

    PetscErrorCode STFilterSetType(SlepcST,SlepcSTFilterType);
    PetscErrorCode STFilterGetType(SlepcST,SlepcSTFilterType*);
    PetscErrorCode STFilterSetInterval(SlepcST,PetscReal,PetscReal)
    PetscErrorCode STFilterGetInterval(SlepcST,PetscReal*,PetscReal*)
    PetscErrorCode STFilterSetRange(SlepcST,PetscReal,PetscReal)
    PetscErrorCode STFilterGetRange(SlepcST,PetscReal*,PetscReal*)
    PetscErrorCode STFilterSetDegree(SlepcST,PetscInt)
    PetscErrorCode STFilterGetDegree(SlepcST,PetscInt*)
    PetscErrorCode STFilterSetDamping(SlepcST,SlepcSTFilterDamping);
    PetscErrorCode STFilterGetDamping(SlepcST,SlepcSTFilterDamping*);

