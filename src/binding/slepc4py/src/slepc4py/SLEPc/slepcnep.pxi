cdef extern from * nogil:

    ctypedef char* SlepcNEPType "const char*"
    SlepcNEPType NEPRII
    SlepcNEPType NEPSLP
    SlepcNEPType NEPNARNOLDI
    SlepcNEPType NEPCISS
    SlepcNEPType NEPINTERPOL
    SlepcNEPType NEPNLEIGS

    ctypedef enum SlepcNEPProblemType "NEPProblemType":
        NEP_GENERAL
        NEP_RATIONAL

    ctypedef enum SlepcNEPWhich "NEPWhich":
        NEP_LARGEST_MAGNITUDE
        NEP_SMALLEST_MAGNITUDE
        NEP_LARGEST_REAL
        NEP_SMALLEST_REAL
        NEP_LARGEST_IMAGINARY
        NEP_SMALLEST_IMAGINARY
        NEP_TARGET_MAGNITUDE
        NEP_TARGET_REAL
        NEP_TARGET_IMAGINARY
        NEP_ALL
        NEP_WHICH_USER

    ctypedef enum SlepcNEPErrorType "NEPErrorType":
        NEP_ERROR_ABSOLUTE
        NEP_ERROR_RELATIVE
        NEP_ERROR_BACKWARD

    ctypedef enum SlepcNEPRefine "NEPRefine":
        NEP_REFINE_NONE
        NEP_REFINE_SIMPLE
        NEP_REFINE_MULTIPLE

    ctypedef enum SlepcNEPRefineScheme "NEPRefineScheme":
        NEP_REFINE_SCHEME_SCHUR
        NEP_REFINE_SCHEME_MBE
        NEP_REFINE_SCHEME_EXPLICIT

    ctypedef enum SlepcNEPConv "NEPConv":
        NEP_CONV_ABS
        NEP_CONV_REL
        NEP_CONV_NORM
        NEP_CONV_USER

    ctypedef enum SlepcNEPStop "NEPStop":
        NEP_STOP_BASIC
        NEP_STOP_USER

    ctypedef enum SlepcNEPConvergedReason "NEPConvergedReason":
        NEP_CONVERGED_TOL
        NEP_CONVERGED_USER
        NEP_DIVERGED_ITS
        NEP_DIVERGED_BREAKDOWN
        NEP_DIVERGED_LINEAR_SOLVE
        NEP_DIVERGED_SUBSPACE_EXHAUSTED
        NEP_CONVERGED_ITERATING

    ctypedef PetscErrorCode (*SlepcNEPFunction)(SlepcNEP,
                                     PetscScalar,
                                     PetscMat,
                                     PetscMat,
                                     void*) except PETSC_ERR_PYTHON

    ctypedef PetscErrorCode (*SlepcNEPJacobian)(SlepcNEP,
                                     PetscScalar,
                                     PetscMat,
                                     void*) except PETSC_ERR_PYTHON

    ctypedef PetscErrorCode (*SlepcNEPCtxDel)(void*)
    ctypedef PetscErrorCode (*SlepcNEPStoppingFunction)(SlepcNEP,
                                             PetscInt,
                                             PetscInt,
                                             PetscInt,
                                             PetscInt,
                                             SlepcNEPConvergedReason*,
                                             void*) except PETSC_ERR_PYTHON
    ctypedef PetscErrorCode (*SlepcNEPMonitorFunction)(SlepcNEP,
                                            PetscInt,
                                            PetscInt,
                                            PetscScalar*,
                                            PetscScalar*,
                                            PetscReal*,
                                            PetscInt,
                                            void*) except PETSC_ERR_PYTHON

    PetscErrorCode NEPCreate(MPI_Comm,SlepcNEP*)
    PetscErrorCode NEPDestroy(SlepcNEP*)
    PetscErrorCode NEPReset(SlepcNEP)
    PetscErrorCode NEPView(SlepcNEP,PetscViewer)

    PetscErrorCode NEPSetType(SlepcNEP,SlepcNEPType)
    PetscErrorCode NEPGetType(SlepcNEP,SlepcNEPType*)
    PetscErrorCode NEPSetTarget(SlepcNEP,PetscScalar)
    PetscErrorCode NEPGetTarget(SlepcNEP,PetscScalar*)
    PetscErrorCode NEPSetOptionsPrefix(SlepcNEP,char*)
    PetscErrorCode NEPGetOptionsPrefix(SlepcNEP,char*[])
    PetscErrorCode NEPSetFromOptions(SlepcNEP)
    PetscErrorCode NEPAppendOptionsPrefix(SlepcNEP,char*)
    PetscErrorCode NEPSetUp(SlepcNEP)
    PetscErrorCode NEPSolve(SlepcNEP)

    PetscErrorCode NEPSetFunction(SlepcNEP,PetscMat,PetscMat,SlepcNEPFunction,void*)
    PetscErrorCode NEPGetFunction(SlepcNEP,PetscMat*,PetscMat*,SlepcNEPFunction*,void**)
    PetscErrorCode NEPSetJacobian(SlepcNEP,PetscMat,SlepcNEPJacobian,void*)
    PetscErrorCode NEPGetJacobian(SlepcNEP,PetscMat*,SlepcNEPJacobian*,void**)
    PetscErrorCode NEPSetSplitOperator(SlepcNEP,PetscInt,PetscMat[],SlepcFN[],PetscMatStructure)
    PetscErrorCode NEPGetSplitOperatorTerm(SlepcNEP,PetscInt,PetscMat*,SlepcFN*)
    PetscErrorCode NEPGetSplitOperatorInfo(SlepcNEP,PetscInt*,PetscMatStructure*)
    PetscErrorCode NEPSetSplitPreconditioner(SlepcNEP,PetscInt,PetscMat[],PetscMatStructure)
    PetscErrorCode NEPGetSplitPreconditionerTerm(SlepcNEP,PetscInt,PetscMat*)
    PetscErrorCode NEPGetSplitPreconditionerInfo(SlepcNEP,PetscInt*,PetscMatStructure*)

    PetscErrorCode NEPSetBV(SlepcNEP,SlepcBV)
    PetscErrorCode NEPGetBV(SlepcNEP,SlepcBV*)
    PetscErrorCode NEPSetRG(SlepcNEP,SlepcRG)
    PetscErrorCode NEPGetRG(SlepcNEP,SlepcRG*)
    PetscErrorCode NEPSetDS(SlepcNEP,SlepcDS)
    PetscErrorCode NEPGetDS(SlepcNEP,SlepcDS*)
    PetscErrorCode NEPSetTolerances(SlepcNEP,PetscReal,PetscInt)
    PetscErrorCode NEPGetTolerances(SlepcNEP,PetscReal*,PetscInt*)

    PetscErrorCode NEPSetTwoSided(SlepcNEP,PetscBool)
    PetscErrorCode NEPGetTwoSided(SlepcNEP,PetscBool*)
    PetscErrorCode NEPApplyResolvent(SlepcNEP,SlepcRG,PetscScalar,PetscVec,PetscVec)

    PetscErrorCode NEPSetTrackAll(SlepcNEP,PetscBool)
    PetscErrorCode NEPGetTrackAll(SlepcNEP,PetscBool*)

    PetscErrorCode NEPSetDimensions(SlepcNEP,PetscInt,PetscInt,PetscInt)
    PetscErrorCode NEPGetDimensions(SlepcNEP,PetscInt*,PetscInt*,PetscInt*)

    PetscErrorCode NEPGetConverged(SlepcNEP,PetscInt*)
    PetscErrorCode NEPGetEigenpair(SlepcNEP,PetscInt,PetscScalar*,PetscScalar*,PetscVec,PetscVec)
    PetscErrorCode NEPGetLeftEigenvector(SlepcNEP,PetscInt,PetscVec,PetscVec)
    PetscErrorCode NEPComputeError(SlepcNEP,PetscInt,SlepcNEPErrorType,PetscReal*)
    PetscErrorCode NEPErrorView(SlepcNEP,SlepcNEPErrorType,PetscViewer)
    PetscErrorCode NEPValuesView(SlepcNEP,PetscViewer)
    PetscErrorCode NEPVectorsView(SlepcNEP,PetscViewer)
    PetscErrorCode NEPGetErrorEstimate(SlepcNEP,PetscInt,PetscReal*)

    PetscErrorCode NEPMonitorSet(SlepcNEP,SlepcNEPMonitorFunction,void*,SlepcNEPCtxDel)
    PetscErrorCode NEPMonitorCancel(SlepcNEP)
    PetscErrorCode NEPGetIterationNumber(SlepcNEP,PetscInt*)

    PetscErrorCode NEPSetInitialSpace(SlepcNEP,PetscInt,PetscVec*)
    PetscErrorCode NEPSetProblemType(SlepcNEP,SlepcNEPProblemType)
    PetscErrorCode NEPGetProblemType(SlepcNEP,SlepcNEPProblemType*)
    PetscErrorCode NEPSetWhichEigenpairs(SlepcNEP,SlepcNEPWhich)
    PetscErrorCode NEPGetWhichEigenpairs(SlepcNEP,SlepcNEPWhich*)

    PetscErrorCode NEPSetRefine(SlepcNEP,SlepcNEPRefine,PetscInt,PetscReal,PetscInt,SlepcNEPRefineScheme)
    PetscErrorCode NEPGetRefine(SlepcNEP,SlepcNEPRefine*,PetscInt*,PetscReal*,PetscInt*,SlepcNEPRefineScheme*)
    PetscErrorCode NEPRefineGetKSP(SlepcNEP,PetscKSP*)

    PetscErrorCode NEPGetConvergedReason(SlepcNEP,SlepcNEPConvergedReason*)
    PetscErrorCode NEPSetConvergenceTest(SlepcNEP,SlepcNEPConv)
    PetscErrorCode NEPGetConvergenceTest(SlepcNEP,SlepcNEPConv*)

    PetscErrorCode NEPSetStoppingTestFunction(SlepcNEP,SlepcNEPStoppingFunction,void*,SlepcNEPCtxDel)
    PetscErrorCode NEPStoppingBasic(SlepcNEP,PetscInt,PetscInt,PetscInt,PetscInt,SlepcNEPConvergedReason*,void*) except PETSC_ERR_PYTHON

    PetscErrorCode NEPRIISetLagPreconditioner(SlepcNEP,PetscInt)
    PetscErrorCode NEPRIIGetLagPreconditioner(SlepcNEP,PetscInt*)
    PetscErrorCode NEPRIISetConstCorrectionTol(SlepcNEP,PetscBool)
    PetscErrorCode NEPRIIGetConstCorrectionTol(SlepcNEP,PetscBool*)
    PetscErrorCode NEPRIISetMaximumIterations(SlepcNEP,PetscInt)
    PetscErrorCode NEPRIIGetMaximumIterations(SlepcNEP,PetscInt*)
    PetscErrorCode NEPRIISetHermitian(SlepcNEP,PetscBool)
    PetscErrorCode NEPRIIGetHermitian(SlepcNEP,PetscBool*)
    PetscErrorCode NEPRIISetDeflationThreshold(SlepcNEP,PetscReal)
    PetscErrorCode NEPRIIGetDeflationThreshold(SlepcNEP,PetscReal*)
    PetscErrorCode NEPRIISetKSP(SlepcNEP,PetscKSP)
    PetscErrorCode NEPRIIGetKSP(SlepcNEP,PetscKSP*)

    PetscErrorCode NEPSLPSetDeflationThreshold(SlepcNEP,PetscReal)
    PetscErrorCode NEPSLPGetDeflationThreshold(SlepcNEP,PetscReal*)
    PetscErrorCode NEPSLPSetEPS(SlepcNEP,SlepcEPS)
    PetscErrorCode NEPSLPGetEPS(SlepcNEP,SlepcEPS*)
    PetscErrorCode NEPSLPSetEPSLeft(SlepcNEP,SlepcEPS)
    PetscErrorCode NEPSLPGetEPSLeft(SlepcNEP,SlepcEPS*)
    PetscErrorCode NEPSLPSetKSP(SlepcNEP,PetscKSP)
    PetscErrorCode NEPSLPGetKSP(SlepcNEP,PetscKSP*)

    PetscErrorCode NEPNArnoldiSetKSP(SlepcNEP,PetscKSP)
    PetscErrorCode NEPNArnoldiGetKSP(SlepcNEP,PetscKSP*)
    PetscErrorCode NEPNArnoldiSetLagPreconditioner(SlepcNEP,PetscInt)
    PetscErrorCode NEPNArnoldiGetLagPreconditioner(SlepcNEP,PetscInt*)

    PetscErrorCode NEPInterpolSetPEP(SlepcNEP,SlepcPEP)
    PetscErrorCode NEPInterpolGetPEP(SlepcNEP,SlepcPEP*)
    PetscErrorCode NEPInterpolSetInterpolation(SlepcNEP,PetscReal,PetscInt)
    PetscErrorCode NEPInterpolGetInterpolation(SlepcNEP,PetscReal*,PetscInt*)

    PetscErrorCode NEPNLEIGSSetRestart(SlepcNEP,PetscReal)
    PetscErrorCode NEPNLEIGSGetRestart(SlepcNEP,PetscReal*)
    PetscErrorCode NEPNLEIGSSetLocking(SlepcNEP,PetscBool)
    PetscErrorCode NEPNLEIGSGetLocking(SlepcNEP,PetscBool*)
    PetscErrorCode NEPNLEIGSSetInterpolation(SlepcNEP,PetscReal,PetscInt)
    PetscErrorCode NEPNLEIGSGetInterpolation(SlepcNEP,PetscReal*,PetscInt*)
    PetscErrorCode NEPNLEIGSSetRKShifts(SlepcNEP,PetscInt,PetscScalar[])
    PetscErrorCode NEPNLEIGSGetRKShifts(SlepcNEP,PetscInt*,PetscScalar*[])
    PetscErrorCode NEPNLEIGSGetKSPs(SlepcNEP,PetscInt*,PetscKSP**)
    PetscErrorCode NEPNLEIGSSetFullBasis(SlepcNEP,PetscBool)
    PetscErrorCode NEPNLEIGSGetFullBasis(SlepcNEP,PetscBool*)
    PetscErrorCode NEPNLEIGSSetEPS(SlepcNEP,SlepcEPS)
    PetscErrorCode NEPNLEIGSGetEPS(SlepcNEP,SlepcEPS*)

    ctypedef enum SlepcNEPCISSExtraction "NEPCISSExtraction":
        NEP_CISS_EXTRACTION_RITZ
        NEP_CISS_EXTRACTION_HANKEL
        NEP_CISS_EXTRACTION_CAA

    PetscErrorCode NEPCISSSetExtraction(SlepcNEP,SlepcNEPCISSExtraction)
    PetscErrorCode NEPCISSGetExtraction(SlepcNEP,SlepcNEPCISSExtraction*)
    PetscErrorCode NEPCISSSetSizes(SlepcNEP,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscBool)
    PetscErrorCode NEPCISSGetSizes(SlepcNEP,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscBool*)
    PetscErrorCode NEPCISSSetThreshold(SlepcNEP,PetscReal,PetscReal)
    PetscErrorCode NEPCISSGetThreshold(SlepcNEP,PetscReal*,PetscReal*)
    PetscErrorCode NEPCISSSetRefinement(SlepcNEP,PetscInt,PetscInt)
    PetscErrorCode NEPCISSGetRefinement(SlepcNEP,PetscInt*,PetscInt*)
    PetscErrorCode NEPCISSGetKSPs(SlepcNEP,PetscInt*,PetscKSP**)

# -----------------------------------------------------------------------------

cdef inline NEP ref_NEP(SlepcNEP nep):
    cdef NEP ob = <NEP> NEP()
    ob.nep = nep
    CHKERR( PetscINCREF(ob.obj) )
    return ob

# -----------------------------------------------------------------------------

cdef PetscErrorCode NEP_Function(
    SlepcNEP    nep,
    PetscScalar mu,
    PetscMat    A,
    PetscMat    B,
    void*       ctx,
    ) except PETSC_ERR_PYTHON with gil:
    cdef NEP Nep  = ref_NEP(nep)
    cdef Mat Amat = ref_Mat(A)
    cdef Mat Bmat = ref_Mat(B)
    (function, args, kargs) = Nep.get_attr('__function__')
    retv = function(Nep, toScalar(mu), Amat, Bmat, *args, **kargs)
    cdef PetscMat Atmp = NULL, Btmp = NULL
    Atmp = A; A = Amat.mat; Amat.mat = Atmp
    Btmp = B; B = Bmat.mat; Bmat.mat = Btmp
    return PETSC_SUCCESS

# -----------------------------------------------------------------------------

cdef PetscErrorCode NEP_Jacobian(
    SlepcNEP    nep,
    PetscScalar mu,
    PetscMat    J,
    void*       ctx,
    ) except PETSC_ERR_PYTHON with gil:
    cdef NEP Nep  = ref_NEP(nep)
    cdef Mat Jmat = ref_Mat(J)
    (jacobian, args, kargs) = Nep.get_attr('__jacobian__')
    retv = jacobian(Nep, toScalar(mu), Jmat, *args, **kargs)
    cdef PetscMat Jtmp = NULL
    Jtmp = J; J = Jmat.mat; Jmat.mat = Jtmp
    return PETSC_SUCCESS

# -----------------------------------------------------------------------------

cdef PetscErrorCode NEP_Stopping(
    SlepcNEP                nep,
    PetscInt                its,
    PetscInt                max_it,
    PetscInt                nconv,
    PetscInt                nev,
    SlepcNEPConvergedReason *r,
    void                    *ctx,
    ) except PETSC_ERR_PYTHON with gil:
    cdef NEP Nep = ref_NEP(nep)
    (stopping, args, kargs) = Nep.get_attr('__stopping__')
    reason = stopping(Nep, toInt(its), toInt(max_it), toInt(nconv), toInt(nev), *args, **kargs)
    if   reason is None:  r[0] = NEP_CONVERGED_ITERATING
    elif reason is False: r[0] = NEP_CONVERGED_ITERATING
    elif reason is True:  r[0] = NEP_CONVERGED_USER
    else:                 r[0] = reason
    return PETSC_SUCCESS

# -----------------------------------------------------------------------------

cdef PetscErrorCode NEP_Monitor(
    SlepcNEP    nep,
    PetscInt    its,
    PetscInt    nconv,
    PetscScalar *eigr,
    PetscScalar *eigi,
    PetscReal   *errest,
    PetscInt    nest,
    void        *ctx,
    ) except PETSC_ERR_PYTHON with gil:
    cdef NEP Nep = ref_NEP(nep)
    cdef object monitorlist = Nep.get_attr('__monitor__')
    if monitorlist is None: return PETSC_SUCCESS
    cdef object eig = [toComplex(eigr[i], eigi[i]) for i in range(nest)]
    cdef object err = [toReal(errest[i]) for i in range(nest)]
    for (monitor, args, kargs) in monitorlist:
        monitor(Nep, toInt(its), toInt(nconv), eig, err, *args, **kargs)
    return PETSC_SUCCESS

# -----------------------------------------------------------------------------
