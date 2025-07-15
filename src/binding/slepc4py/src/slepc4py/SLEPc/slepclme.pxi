cdef extern from * nogil:

    ctypedef char* SlepcLMEType "const char*"
    SlepcLMEType LMEKRYLOV

    ctypedef enum SlepcLMEConvergedReason "LMEConvergedReason":
        LME_CONVERGED_TOL
        LME_DIVERGED_ITS
        LME_DIVERGED_BREAKDOWN
        LME_CONVERGED_ITERATING

    ctypedef enum SlepcLMEProblemType "LMEProblemType":
        LME_LYAPUNOV
        LME_SYLVESTER
        LME_GEN_LYAPUNOV
        LME_GEN_SYLVESTER
        LME_DT_LYAPUNOV
        LME_STEIN

    ctypedef PetscErrorCode (*SlepcLMECtxDel)(void*)
    ctypedef PetscErrorCode (*SlepcLMEMonitorFunction)(SlepcLME,
                                            PetscInt,
                                            PetscReal,
                                            void*) except PETSC_ERR_PYTHON

    PetscErrorCode LMECreate(MPI_Comm,SlepcLME*)
    PetscErrorCode LMEDestroy(SlepcLME*)
    PetscErrorCode LMEReset(SlepcLME)
    PetscErrorCode LMEView(SlepcLME,PetscViewer)

    PetscErrorCode LMESetType(SlepcLME,SlepcLMEType)
    PetscErrorCode LMEGetType(SlepcLME,SlepcLMEType*)
    PetscErrorCode LMESetProblemType(SlepcLME,SlepcLMEProblemType)
    PetscErrorCode LMEGetProblemType(SlepcLME,SlepcLMEProblemType*)

    PetscErrorCode LMESetCoefficients(SlepcLME,PetscMat,PetscMat,PetscMat,PetscMat)
    PetscErrorCode LMEGetCoefficients(SlepcLME,PetscMat*,PetscMat*,PetscMat*,PetscMat*)
    PetscErrorCode LMESetRHS(SlepcLME,PetscMat)
    PetscErrorCode LMEGetRHS(SlepcLME,PetscMat*)
    PetscErrorCode LMESetSolution(SlepcLME,PetscMat)
    PetscErrorCode LMEGetSolution(SlepcLME,PetscMat*)
    PetscErrorCode LMEGetErrorEstimate(SlepcLME,PetscReal*)
    PetscErrorCode LMEComputeError(SlepcLME,PetscReal*)

    PetscErrorCode LMESetOptionsPrefix(SlepcLME,char*)
    PetscErrorCode LMEGetOptionsPrefix(SlepcLME,char*[])
    PetscErrorCode LMESetFromOptions(SlepcLME)
    PetscErrorCode LMEAppendOptionsPrefix(SlepcLME,char*)
    PetscErrorCode LMESetUp(SlepcLME)
    PetscErrorCode LMESolve(SlepcLME)

    PetscErrorCode LMESetBV(SlepcLME,SlepcBV)
    PetscErrorCode LMEGetBV(SlepcLME,SlepcBV*)
    PetscErrorCode LMESetFN(SlepcLME,SlepcFN)
    PetscErrorCode LMEGetFN(SlepcLME,SlepcFN*)
    PetscErrorCode LMESetTolerances(SlepcLME,PetscReal,PetscInt)
    PetscErrorCode LMEGetTolerances(SlepcLME,PetscReal*,PetscInt*)
    PetscErrorCode LMESetDimensions(SlepcLME,PetscInt)
    PetscErrorCode LMEGetDimensions(SlepcLME,PetscInt*)

    PetscErrorCode LMESetErrorIfNotConverged(SlepcLME,PetscBool)
    PetscErrorCode LMEGetErrorIfNotConverged(SlepcLME,PetscBool*)

    PetscErrorCode LMEMonitorSet(SlepcLME,SlepcLMEMonitorFunction,void*,SlepcLMECtxDel)
    PetscErrorCode LMEMonitorCancel(SlepcLME)
    PetscErrorCode LMEGetIterationNumber(SlepcLME,PetscInt*)

    PetscErrorCode LMEGetConvergedReason(SlepcLME,SlepcLMEConvergedReason*)

# -----------------------------------------------------------------------------

cdef inline LME ref_LME(SlepcLME lme):
    cdef LME ob = <LME> LME()
    ob.lme = lme
    CHKERR( PetscINCREF(ob.obj) )
    return ob

# -----------------------------------------------------------------------------

cdef PetscErrorCode LME_Monitor(
    SlepcLME    lme,
    PetscInt    it,
    PetscReal   errest,
    void        *ctx,
    ) except PETSC_ERR_PYTHON with gil:
    cdef LME Lme = ref_LME(lme)
    cdef object monitorlist = Lme.get_attr('__monitor__')
    if monitorlist is None: return PETSC_SUCCESS
    for (monitor, args, kargs) in monitorlist:
        monitor(Lme, toInt(it), toReal(errest), *args, **kargs)
    return PETSC_SUCCESS

# -----------------------------------------------------------------------------
