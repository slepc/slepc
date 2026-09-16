cdef extern from * nogil:

    ctypedef char* SlepcMFNType "const char*"
    SlepcMFNType MFNKRYLOV
    SlepcMFNType MFNEXPOKIT

    ctypedef enum SlepcMFNConvergedReason "MFNConvergedReason":
        MFN_CONVERGED_TOL
        MFN_CONVERGED_ITS
        MFN_DIVERGED_ITS
        MFN_DIVERGED_BREAKDOWN
        MFN_CONVERGED_ITERATING

    ctypedef PetscErrorCode (*SlepcMFNCtxDel)(void*)
    ctypedef PetscErrorCode (*SlepcMFNMonitorFunction)(SlepcMFN,
                                                       PetscInt,
                                                       PetscReal,
                                                       void*) except PETSC_ERR_PYTHON

    PetscErrorCode MFNCreate(MPI_Comm, SlepcMFN*)
    PetscErrorCode MFNDestroy(SlepcMFN*)
    PetscErrorCode MFNReset(SlepcMFN)
    PetscErrorCode MFNView(SlepcMFN, PetscViewer)

    PetscErrorCode MFNSetType(SlepcMFN, SlepcMFNType)
    PetscErrorCode MFNGetType(SlepcMFN, SlepcMFNType*)
    PetscErrorCode MFNSetOperator(SlepcMFN, PetscMat)
    PetscErrorCode MFNGetOperator(SlepcMFN, PetscMat*)
    PetscErrorCode MFNSetOptionsPrefix(SlepcMFN, char*)
    PetscErrorCode MFNGetOptionsPrefix(SlepcMFN, char*[])
    PetscErrorCode MFNSetFromOptions(SlepcMFN)
    PetscErrorCode MFNAppendOptionsPrefix(SlepcMFN, char*)
    PetscErrorCode MFNSetUp(SlepcMFN)
    PetscErrorCode MFNSolve(SlepcMFN, PetscVec, PetscVec)
    PetscErrorCode MFNSolveTranspose(SlepcMFN, PetscVec, PetscVec)

    PetscErrorCode MFNSetBV(SlepcMFN, SlepcBV)
    PetscErrorCode MFNGetBV(SlepcMFN, SlepcBV*)
    PetscErrorCode MFNSetFN(SlepcMFN, SlepcFN)
    PetscErrorCode MFNGetFN(SlepcMFN, SlepcFN*)
    PetscErrorCode MFNSetTolerances(SlepcMFN, PetscReal, PetscInt)
    PetscErrorCode MFNGetTolerances(SlepcMFN, PetscReal*, PetscInt*)
    PetscErrorCode MFNSetDimensions(SlepcMFN, PetscInt)
    PetscErrorCode MFNGetDimensions(SlepcMFN, PetscInt*)

    PetscErrorCode MFNSetErrorIfNotConverged(SlepcMFN, PetscBool)
    PetscErrorCode MFNGetErrorIfNotConverged(SlepcMFN, PetscBool*)

    PetscErrorCode MFNMonitorSet(SlepcMFN, SlepcMFNMonitorFunction, void*, SlepcMFNCtxDel)
    PetscErrorCode MFNMonitorCancel(SlepcMFN)
    PetscErrorCode MFNGetIterationNumber(SlepcMFN, PetscInt*)

    PetscErrorCode MFNGetConvergedReason(SlepcMFN, SlepcMFNConvergedReason*)

# -----------------------------------------------------------------------------

cdef inline MFN ref_MFN(SlepcMFN mfn):
    cdef MFN ob = <MFN> MFN()
    ob.mfn = mfn
    CHKERR(PetscINCREF(ob.obj))
    return ob

# -----------------------------------------------------------------------------

cdef PetscErrorCode MFN_Monitor(
    SlepcMFN    mfn,
    PetscInt    it,
    PetscReal   errest,
    void        *ctx,
    ) except PETSC_ERR_PYTHON with gil:
    cdef MFN Mfn = ref_MFN(mfn)
    cdef object monitorlist = Mfn.get_attr('__monitor__')
    if monitorlist is None: return PETSC_SUCCESS
    for (monitor, args, kargs) in monitorlist:
        monitor(Mfn, toInt(it), toReal(errest), *args, **kargs)
    return PETSC_SUCCESS

# -----------------------------------------------------------------------------
