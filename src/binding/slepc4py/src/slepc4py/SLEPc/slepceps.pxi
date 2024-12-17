cdef extern from * nogil:

    ctypedef char* SlepcEPSType "const char*"
    SlepcEPSType EPSPOWER
    SlepcEPSType EPSSUBSPACE
    SlepcEPSType EPSARNOLDI
    SlepcEPSType EPSLANCZOS
    SlepcEPSType EPSKRYLOVSCHUR
    SlepcEPSType EPSGD
    SlepcEPSType EPSJD
    SlepcEPSType EPSRQCG
    SlepcEPSType EPSLOBPCG
    SlepcEPSType EPSCISS
    SlepcEPSType EPSLYAPII
    SlepcEPSType EPSLAPACK
    SlepcEPSType EPSARPACK
    SlepcEPSType EPSTRLAN
    SlepcEPSType EPSBLOPEX
    SlepcEPSType EPSPRIMME
    SlepcEPSType EPSFEAST
    SlepcEPSType EPSSCALAPACK
    SlepcEPSType EPSELPA
    SlepcEPSType EPSELEMENTAL
    SlepcEPSType EPSEVSL
    SlepcEPSType EPSCHASE

    ctypedef enum SlepcEPSProblemType "EPSProblemType":
        EPS_HEP
        EPS_GHEP
        EPS_NHEP
        EPS_GNHEP
        EPS_PGNHEP
        EPS_GHIEP
        EPS_BSE

    ctypedef enum SlepcEPSExtraction "EPSExtraction":
        EPS_RITZ
        EPS_HARMONIC
        EPS_HARMONIC_RELATIVE
        EPS_HARMONIC_RIGHT
        EPS_HARMONIC_LARGEST
        EPS_REFINED
        EPS_REFINED_HARMONIC

    ctypedef enum SlepcEPSWhich "EPSWhich":
        EPS_LARGEST_MAGNITUDE
        EPS_LARGEST_REAL
        EPS_LARGEST_IMAGINARY
        EPS_SMALLEST_MAGNITUDE
        EPS_SMALLEST_REAL
        EPS_SMALLEST_IMAGINARY
        EPS_TARGET_MAGNITUDE
        EPS_TARGET_REAL
        EPS_TARGET_IMAGINARY
        EPS_ALL
        EPS_WHICH_USER

    ctypedef enum SlepcEPSBalance "EPSBalance":
        EPS_BALANCE_NONE
        EPS_BALANCE_ONESIDE
        EPS_BALANCE_TWOSIDE
        EPS_BALANCE_USER

    ctypedef enum SlepcEPSErrorType "EPSErrorType":
        EPS_ERROR_ABSOLUTE
        EPS_ERROR_RELATIVE
        EPS_ERROR_BACKWARD

    ctypedef enum SlepcEPSConv "EPSConv":
        EPS_CONV_ABS
        EPS_CONV_REL
        EPS_CONV_NORM
        EPS_CONV_USER

    ctypedef enum SlepcEPSStop "EPSStop":
        EPS_STOP_BASIC
        EPS_STOP_USER
        EPS_STOP_THRESHOLD

    ctypedef enum SlepcEPSConvergedReason "EPSConvergedReason":
        EPS_CONVERGED_TOL
        EPS_CONVERGED_USER
        EPS_DIVERGED_ITS
        EPS_DIVERGED_BREAKDOWN
        EPS_DIVERGED_SYMMETRY_LOST
        EPS_CONVERGED_ITERATING

    ctypedef PetscErrorCode (*SlepcEPSCtxDel)(void*)
    ctypedef PetscErrorCode (*SlepcEPSStoppingFunction)(SlepcEPS,
                                             PetscInt,
                                             PetscInt,
                                             PetscInt,
                                             PetscInt,
                                             SlepcEPSConvergedReason*,
                                             void*) except PETSC_ERR_PYTHON
    ctypedef PetscErrorCode (*SlepcEPSMonitorFunction)(SlepcEPS,
                                            PetscInt,
                                            PetscInt,
                                            PetscScalar*,
                                            PetscScalar*,
                                            PetscReal*,
                                            PetscInt,
                                            void*) except PETSC_ERR_PYTHON
    ctypedef PetscErrorCode (*SlepcEPSArbitraryFunction)(PetscScalar,
                                             PetscScalar,
                                             PetscVec,
                                             PetscVec,
                                             PetscScalar*,
                                             PetscScalar*,
                                             void*) except PETSC_ERR_PYTHON
    ctypedef PetscErrorCode (*SlepcEPSComparisonFunction)(PetscScalar,
                                               PetscScalar,
                                               PetscScalar,
                                               PetscScalar,
                                               PetscInt*,
                                               void*) except PETSC_ERR_PYTHON

    PetscErrorCode EPSView(SlepcEPS,PetscViewer)
    PetscErrorCode EPSDestroy(SlepcEPS*)
    PetscErrorCode EPSReset(SlepcEPS)
    PetscErrorCode EPSCreate(MPI_Comm,SlepcEPS*)
    PetscErrorCode EPSSetType(SlepcEPS,SlepcEPSType)
    PetscErrorCode EPSGetType(SlepcEPS,SlepcEPSType*)
    PetscErrorCode EPSSetOptionsPrefix(SlepcEPS,char[])
    PetscErrorCode EPSAppendOptionsPrefix(SlepcEPS,char [])
    PetscErrorCode EPSGetOptionsPrefix(SlepcEPS,char*[])
    PetscErrorCode EPSSetFromOptions(SlepcEPS)

    PetscErrorCode EPSSetProblemType(SlepcEPS,SlepcEPSProblemType)
    PetscErrorCode EPSGetProblemType(SlepcEPS,SlepcEPSProblemType*)
    PetscErrorCode EPSIsGeneralized(SlepcEPS,PetscBool*)
    PetscErrorCode EPSIsHermitian(SlepcEPS,PetscBool*)
    PetscErrorCode EPSIsPositive(SlepcEPS,PetscBool*)
    PetscErrorCode EPSIsStructured(SlepcEPS,PetscBool*)
    PetscErrorCode EPSSetExtraction(SlepcEPS,SlepcEPSExtraction)
    PetscErrorCode EPSGetExtraction(SlepcEPS,SlepcEPSExtraction*)
    PetscErrorCode EPSSetBalance(SlepcEPS,SlepcEPSBalance,PetscInt,PetscReal)
    PetscErrorCode EPSGetBalance(SlepcEPS,SlepcEPSBalance*,PetscInt*,PetscReal*)
    PetscErrorCode EPSSetWhichEigenpairs(SlepcEPS,SlepcEPSWhich)
    PetscErrorCode EPSGetWhichEigenpairs(SlepcEPS,SlepcEPSWhich*)
    PetscErrorCode EPSSetThreshold(SlepcEPS,PetscReal,PetscBool)
    PetscErrorCode EPSGetThreshold(SlepcEPS,PetscReal*,PetscBool*)
    PetscErrorCode EPSSetTarget(SlepcEPS,PetscScalar)
    PetscErrorCode EPSGetTarget(SlepcEPS,PetscScalar*)
    PetscErrorCode EPSSetInterval(SlepcEPS,PetscReal,PetscReal)
    PetscErrorCode EPSGetInterval(SlepcEPS,PetscReal*,PetscReal*)
    PetscErrorCode EPSSetTolerances(SlepcEPS,PetscReal,PetscInt)
    PetscErrorCode EPSGetTolerances(SlepcEPS,PetscReal*,PetscInt*)
    PetscErrorCode EPSSetDimensions(SlepcEPS,PetscInt,PetscInt,PetscInt)
    PetscErrorCode EPSGetDimensions(SlepcEPS,PetscInt*,PetscInt*,PetscInt*)

    PetscErrorCode EPSSetBV(SlepcEPS,SlepcBV)
    PetscErrorCode EPSGetBV(SlepcEPS,SlepcBV*)
    PetscErrorCode EPSSetDS(SlepcEPS,SlepcDS)
    PetscErrorCode EPSGetDS(SlepcEPS,SlepcDS*)
    PetscErrorCode EPSSetST(SlepcEPS,SlepcST)
    PetscErrorCode EPSGetST(SlepcEPS,SlepcST*)
    PetscErrorCode EPSSetRG(SlepcEPS,SlepcRG)
    PetscErrorCode EPSGetRG(SlepcEPS,SlepcRG*)

    PetscErrorCode EPSSetOperators(SlepcEPS,PetscMat,PetscMat)
    PetscErrorCode EPSGetOperators(SlepcEPS,PetscMat*,PetscMat*)

    PetscErrorCode EPSSetTwoSided(SlepcEPS,PetscBool)
    PetscErrorCode EPSGetTwoSided(SlepcEPS,PetscBool*)
    PetscErrorCode EPSSetPurify(SlepcEPS,PetscBool)
    PetscErrorCode EPSGetPurify(SlepcEPS,PetscBool*)

    PetscErrorCode EPSSetConvergenceTest(SlepcEPS,SlepcEPSConv)
    PetscErrorCode EPSGetConvergenceTest(SlepcEPS,SlepcEPSConv*)

    PetscErrorCode EPSSetTrueResidual(SlepcEPS,PetscBool)
    PetscErrorCode EPSGetTrueResidual(SlepcEPS,PetscBool*)

    PetscErrorCode EPSSetTrackAll(SlepcEPS,PetscBool)
    PetscErrorCode EPSGetTrackAll(SlepcEPS,PetscBool*)

    PetscErrorCode EPSSetDeflationSpace(SlepcEPS,PetscInt,PetscVec*)
    PetscErrorCode EPSSetInitialSpace(SlepcEPS,PetscInt,PetscVec*)
    PetscErrorCode EPSSetLeftInitialSpace(SlepcEPS,PetscInt,PetscVec*)

    PetscErrorCode EPSSetUp(SlepcEPS)
    PetscErrorCode EPSSolve(SlepcEPS)

    PetscErrorCode EPSGetIterationNumber(SlepcEPS,PetscInt*)
    PetscErrorCode EPSGetConvergedReason(SlepcEPS,SlepcEPSConvergedReason*)
    PetscErrorCode EPSGetConverged(SlepcEPS,PetscInt*)
    PetscErrorCode EPSGetEigenvalue(SlepcEPS,PetscInt,PetscScalar*,PetscScalar*)
    PetscErrorCode EPSGetEigenvector(SlepcEPS,PetscInt,PetscVec,PetscVec)
    PetscErrorCode EPSGetLeftEigenvector(SlepcEPS,PetscInt,PetscVec,PetscVec)
    PetscErrorCode EPSGetEigenpair(SlepcEPS,PetscInt,PetscScalar*,PetscScalar*,PetscVec,PetscVec)
    PetscErrorCode EPSGetInvariantSubspace(SlepcEPS,PetscVec*)

    PetscErrorCode EPSSetStoppingTestFunction(SlepcEPS,SlepcEPSStoppingFunction,void*,SlepcEPSCtxDel)
    PetscErrorCode EPSStoppingBasic(SlepcEPS,PetscInt,PetscInt,PetscInt,PetscInt,SlepcEPSConvergedReason*,void*) except PETSC_ERR_PYTHON

    PetscErrorCode EPSSetArbitrarySelection(SlepcEPS,SlepcEPSArbitraryFunction,void*);
    PetscErrorCode EPSSetEigenvalueComparison(SlepcEPS,SlepcEPSComparisonFunction,void*);

    PetscErrorCode EPSGetErrorEstimate(SlepcEPS,PetscInt,PetscReal*)
    PetscErrorCode EPSComputeError(SlepcEPS,PetscInt,SlepcEPSErrorType,PetscReal*)
    PetscErrorCode EPSErrorView(SlepcEPS,SlepcEPSErrorType,PetscViewer)
    PetscErrorCode EPSValuesView(SlepcEPS,PetscViewer)
    PetscErrorCode EPSVectorsView(SlepcEPS,PetscViewer)

    PetscErrorCode EPSMonitorSet(SlepcEPS,SlepcEPSMonitorFunction,void*,SlepcEPSCtxDel)
    PetscErrorCode EPSMonitorCancel(SlepcEPS)

    ctypedef enum SlepcEPSPowerShiftType "EPSPowerShiftType":
        EPS_POWER_SHIFT_CONSTANT
        EPS_POWER_SHIFT_RAYLEIGH
        EPS_POWER_SHIFT_WILKINSON
    PetscErrorCode EPSPowerSetShiftType(SlepcEPS,SlepcEPSPowerShiftType)
    PetscErrorCode EPSPowerGetShiftType(SlepcEPS,SlepcEPSPowerShiftType*)

    PetscErrorCode EPSArnoldiSetDelayed(SlepcEPS,PetscBool)
    PetscErrorCode EPSArnoldiGetDelayed(SlepcEPS,PetscBool*)

    ctypedef enum SlepcEPSKrylovSchurBSEType "EPSKrylovSchurBSEType":
        EPS_KRYLOVSCHUR_BSE_SHAO
        EPS_KRYLOVSCHUR_BSE_GRUNING
        EPS_KRYLOVSCHUR_BSE_PROJECTEDBSE
    PetscErrorCode EPSKrylovSchurSetBSEType(SlepcEPS,SlepcEPSKrylovSchurBSEType)
    PetscErrorCode EPSKrylovSchurGetBSEType(SlepcEPS,SlepcEPSKrylovSchurBSEType*)
    PetscErrorCode EPSKrylovSchurSetRestart(SlepcEPS,PetscReal)
    PetscErrorCode EPSKrylovSchurGetRestart(SlepcEPS,PetscReal*)
    PetscErrorCode EPSKrylovSchurSetLocking(SlepcEPS,PetscBool)
    PetscErrorCode EPSKrylovSchurGetLocking(SlepcEPS,PetscBool*)
    PetscErrorCode EPSKrylovSchurSetPartitions(SlepcEPS,PetscInt)
    PetscErrorCode EPSKrylovSchurGetPartitions(SlepcEPS,PetscInt*)
    PetscErrorCode EPSKrylovSchurSetDetectZeros(SlepcEPS,PetscBool)
    PetscErrorCode EPSKrylovSchurGetDetectZeros(SlepcEPS,PetscBool*)
    PetscErrorCode EPSKrylovSchurSetDimensions(SlepcEPS,PetscInt,PetscInt,PetscInt)
    PetscErrorCode EPSKrylovSchurGetDimensions(SlepcEPS,PetscInt*,PetscInt*,PetscInt*)
    PetscErrorCode EPSKrylovSchurGetSubcommInfo(SlepcEPS,PetscInt*,PetscInt*,PetscVec*)
    PetscErrorCode EPSKrylovSchurGetSubcommPairs(SlepcEPS,PetscInt,PetscScalar*,PetscVec)
    PetscErrorCode EPSKrylovSchurGetSubcommMats(SlepcEPS,PetscMat*,PetscMat*)
    PetscErrorCode EPSKrylovSchurUpdateSubcommMats(SlepcEPS,PetscScalar,PetscScalar,PetscMat,PetscScalar,PetscScalar,PetscMat,PetscMatStructure,PetscBool)
    PetscErrorCode EPSKrylovSchurSetSubintervals(SlepcEPS,PetscReal*)
    PetscErrorCode EPSKrylovSchurGetSubintervals(SlepcEPS,PetscReal**)
    PetscErrorCode EPSKrylovSchurGetInertias(SlepcEPS,PetscInt*,PetscReal**,PetscInt**)
    PetscErrorCode EPSKrylovSchurGetKSP(SlepcEPS,PetscKSP*)

    ctypedef enum SlepcEPSLanczosReorthogType "EPSLanczosReorthogType":
        EPS_LANCZOS_REORTHOG_LOCAL
        EPS_LANCZOS_REORTHOG_FULL
        EPS_LANCZOS_REORTHOG_SELECTIVE
        EPS_LANCZOS_REORTHOG_PERIODIC
        EPS_LANCZOS_REORTHOG_PARTIAL
        EPS_LANCZOS_REORTHOG_DELAYED
    PetscErrorCode EPSLanczosSetReorthog(SlepcEPS,SlepcEPSLanczosReorthogType)
    PetscErrorCode EPSLanczosGetReorthog(SlepcEPS,SlepcEPSLanczosReorthogType*)

    PetscErrorCode EPSGDSetKrylovStart(SlepcEPS,PetscBool)
    PetscErrorCode EPSGDGetKrylovStart(SlepcEPS,PetscBool*)
    PetscErrorCode EPSGDSetBlockSize(SlepcEPS,PetscInt)
    PetscErrorCode EPSGDGetBlockSize(SlepcEPS,PetscInt*)
    PetscErrorCode EPSGDSetRestart(SlepcEPS,PetscInt,PetscInt)
    PetscErrorCode EPSGDGetRestart(SlepcEPS,PetscInt*,PetscInt*)
    PetscErrorCode EPSGDSetInitialSize(SlepcEPS,PetscInt)
    PetscErrorCode EPSGDGetInitialSize(SlepcEPS,PetscInt*)
    PetscErrorCode EPSGDSetBOrth(SlepcEPS,PetscBool)
    PetscErrorCode EPSGDGetBOrth(SlepcEPS,PetscBool*)
    PetscErrorCode EPSGDSetDoubleExpansion(SlepcEPS,PetscBool)
    PetscErrorCode EPSGDGetDoubleExpansion(SlepcEPS,PetscBool*)

    PetscErrorCode EPSJDSetKrylovStart(SlepcEPS,PetscBool)
    PetscErrorCode EPSJDGetKrylovStart(SlepcEPS,PetscBool*)
    PetscErrorCode EPSJDSetBlockSize(SlepcEPS,PetscInt)
    PetscErrorCode EPSJDGetBlockSize(SlepcEPS,PetscInt*)
    PetscErrorCode EPSJDSetRestart(SlepcEPS,PetscInt,PetscInt)
    PetscErrorCode EPSJDGetRestart(SlepcEPS,PetscInt*,PetscInt*)
    PetscErrorCode EPSJDSetInitialSize(SlepcEPS,PetscInt)
    PetscErrorCode EPSJDGetInitialSize(SlepcEPS,PetscInt*)
    PetscErrorCode EPSJDSetFix(SlepcEPS,PetscReal)
    PetscErrorCode EPSJDGetFix(SlepcEPS,PetscReal*)
    PetscErrorCode EPSJDSetConstCorrectionTol(SlepcEPS,PetscBool)
    PetscErrorCode EPSJDGetConstCorrectionTol(SlepcEPS,PetscBool*)
    PetscErrorCode EPSJDSetBOrth(SlepcEPS,PetscBool)
    PetscErrorCode EPSJDGetBOrth(SlepcEPS,PetscBool*)

    PetscErrorCode EPSRQCGSetReset(SlepcEPS,PetscInt)
    PetscErrorCode EPSRQCGGetReset(SlepcEPS,PetscInt*)

    PetscErrorCode EPSLOBPCGSetBlockSize(SlepcEPS,PetscInt)
    PetscErrorCode EPSLOBPCGGetBlockSize(SlepcEPS,PetscInt*)
    PetscErrorCode EPSLOBPCGSetRestart(SlepcEPS,PetscReal)
    PetscErrorCode EPSLOBPCGGetRestart(SlepcEPS,PetscReal*)
    PetscErrorCode EPSLOBPCGSetLocking(SlepcEPS,PetscBool)
    PetscErrorCode EPSLOBPCGGetLocking(SlepcEPS,PetscBool*)

    PetscErrorCode EPSLyapIISetRanks(SlepcEPS,PetscInt,PetscInt)
    PetscErrorCode EPSLyapIIGetRanks(SlepcEPS,PetscInt*,PetscInt*)

    ctypedef enum SlepcEPSCISSQuadRule "EPSCISSQuadRule":
        EPS_CISS_QUADRULE_TRAPEZOIDAL
        EPS_CISS_QUADRULE_CHEBYSHEV

    ctypedef enum SlepcEPSCISSExtraction "EPSCISSExtraction":
        EPS_CISS_EXTRACTION_RITZ
        EPS_CISS_EXTRACTION_HANKEL

    PetscErrorCode EPSCISSSetExtraction(SlepcEPS,SlepcEPSCISSExtraction)
    PetscErrorCode EPSCISSGetExtraction(SlepcEPS,SlepcEPSCISSExtraction*)
    PetscErrorCode EPSCISSSetQuadRule(SlepcEPS,SlepcEPSCISSQuadRule)
    PetscErrorCode EPSCISSGetQuadRule(SlepcEPS,SlepcEPSCISSQuadRule*)
    PetscErrorCode EPSCISSSetSizes(SlepcEPS,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscBool)
    PetscErrorCode EPSCISSGetSizes(SlepcEPS,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscBool*)
    PetscErrorCode EPSCISSSetThreshold(SlepcEPS,PetscReal,PetscReal)
    PetscErrorCode EPSCISSGetThreshold(SlepcEPS,PetscReal*,PetscReal*)
    PetscErrorCode EPSCISSSetRefinement(SlepcEPS,PetscInt,PetscInt)
    PetscErrorCode EPSCISSGetRefinement(SlepcEPS,PetscInt*,PetscInt*)
    PetscErrorCode EPSCISSSetUseST(SlepcEPS,PetscBool)
    PetscErrorCode EPSCISSGetUseST(SlepcEPS,PetscBool*)
    PetscErrorCode EPSCISSGetKSPs(SlepcEPS,PetscInt*,PetscKSP**)

cdef extern from * nogil:
    PetscErrorCode VecDuplicate(PetscVec,PetscVec*)
    PetscErrorCode MatCreateVecs(PetscMat,PetscVec*,PetscVec*)

# -----------------------------------------------------------------------------

cdef inline EPS ref_EPS(SlepcEPS eps):
    cdef EPS ob = <EPS> EPS()
    ob.eps = eps
    CHKERR( PetscINCREF(ob.obj) )
    return ob

# -----------------------------------------------------------------------------

cdef PetscErrorCode EPS_Stopping(
    SlepcEPS                eps,
    PetscInt                its,
    PetscInt                max_it,
    PetscInt                nconv,
    PetscInt                nev,
    SlepcEPSConvergedReason *r,
    void                    *ctx,
    ) except PETSC_ERR_PYTHON with gil:
    cdef EPS Eps = ref_EPS(eps)
    (stopping, args, kargs) = Eps.get_attr('__stopping__')
    reason = stopping(Eps, toInt(its), toInt(max_it), toInt(nconv), toInt(nev), *args, **kargs)
    if   reason is None:  r[0] = EPS_CONVERGED_ITERATING
    elif reason is False: r[0] = EPS_CONVERGED_ITERATING
    elif reason is True:  r[0] = EPS_CONVERGED_USER
    else:                 r[0] = reason

# -----------------------------------------------------------------------------

cdef PetscErrorCode EPS_Arbitrary(
    PetscScalar  er,
    PetscScalar  ei,
    PetscVec     xr,
    PetscVec     xi,
    PetscScalar* rr,
    PetscScalar* ri,
    void         *ctx,
    ) except PETSC_ERR_PYTHON with gil:
    (arbitrary, args, kargs) = <object>ctx
    cdef Vec Vr = ref_Vec(xr)
    cdef Vec Vi = ref_Vec(xi)
    r = arbitrary(toComplex(er, ei), Vr, Vi, args, **kargs)
    if sizeof(PetscScalar) == sizeof(PetscReal):
        rr[0] = asComplexReal(r)
        ri[0] = asComplexImag(r)
    else:
        rr[0] = asScalar(r)
        ri[0] = 0.0

# -----------------------------------------------------------------------------

cdef PetscErrorCode EPS_Comparison(
    PetscScalar  ar,
    PetscScalar  ai,
    PetscScalar  br,
    PetscScalar  bi,
    PetscInt*    res,
    void         *ctx,
    ) except PETSC_ERR_PYTHON with gil:
    (comparison, args, kargs) = <object>ctx
    r = comparison(toComplex(ar, ai), toComplex(br, bi), args, **kargs)
    res[0] = asInt(r)

# -----------------------------------------------------------------------------

cdef PetscErrorCode EPS_Monitor(
    SlepcEPS    eps,
    PetscInt    its,
    PetscInt    nconv,
    PetscScalar *eigr,
    PetscScalar *eigi,
    PetscReal   *errest,
    PetscInt    nest,
    void        *ctx,
    ) except PETSC_ERR_PYTHON with gil:
    cdef EPS Eps = ref_EPS(eps)
    cdef object monitorlist = Eps.get_attr('__monitor__')
    if monitorlist is None: return PETSC_SUCCESS
    cdef object eig = [toComplex(eigr[i], eigi[i]) for i in range(nest)]
    cdef object err = [toReal(errest[i]) for i in range(nest)]
    for (monitor, args, kargs) in monitorlist:
        monitor(Eps, toInt(its), toInt(nconv), eig, err, *args, **kargs)
    return PETSC_SUCCESS

# -----------------------------------------------------------------------------
