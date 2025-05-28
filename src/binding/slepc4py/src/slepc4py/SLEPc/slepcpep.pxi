cdef extern from * nogil:

    ctypedef char* SlepcPEPType "const char*"
    SlepcPEPType PEPLINEAR
    SlepcPEPType PEPQARNOLDI
    SlepcPEPType PEPTOAR
    SlepcPEPType PEPSTOAR
    SlepcPEPType PEPJD
    SlepcPEPType PEPCISS

    ctypedef enum SlepcPEPProblemType "PEPProblemType":
        PEP_GENERAL
        PEP_HERMITIAN
        PEP_HYPERBOLIC
        PEP_GYROSCOPIC

    ctypedef enum SlepcPEPRefine "PEPRefine":
        PEP_REFINE_NONE
        PEP_REFINE_SIMPLE
        PEP_REFINE_MULTIPLE

    ctypedef enum SlepcPEPExtract "PEPExtract":
        PEP_EXTRACT_NONE
        PEP_EXTRACT_NORM
        PEP_EXTRACT_RESIDUAL
        PEP_EXTRACT_STRUCTURED

    ctypedef enum SlepcPEPRefineScheme "PEPRefineScheme":
        PEP_REFINE_SCHEME_EXPLICIT
        PEP_REFINE_SCHEME_MBE
        PEP_REFINE_SCHEME_SCHUR

    ctypedef enum SlepcPEPErrorType "PEPErrorType":
        PEP_ERROR_ABSOLUTE
        PEP_ERROR_RELATIVE
        PEP_ERROR_BACKWARD

    ctypedef enum SlepcPEPWhich "PEPWhich":
        PEP_LARGEST_MAGNITUDE
        PEP_SMALLEST_MAGNITUDE
        PEP_LARGEST_REAL
        PEP_SMALLEST_REAL
        PEP_LARGEST_IMAGINARY
        PEP_SMALLEST_IMAGINARY
        PEP_TARGET_MAGNITUDE
        PEP_TARGET_REAL
        PEP_TARGET_IMAGINARY
        PEP_ALL
        PEP_WHICH_USER

    ctypedef enum SlepcPEPBasis "PEPBasis":
        PEP_BASIS_MONOMIAL
        PEP_BASIS_CHEBYSHEV1
        PEP_BASIS_CHEBYSHEV2
        PEP_BASIS_LEGENDRE
        PEP_BASIS_LAGUERRE
        PEP_BASIS_HERMITE

    ctypedef enum SlepcPEPScale "PEPScale":
        PEP_SCALE_NONE
        PEP_SCALE_SCALAR
        PEP_SCALE_DIAGONAL
        PEP_SCALE_BOTH

    ctypedef enum SlepcPEPConv "PEPConv":
        PEP_CONV_ABS
        PEP_CONV_REL
        PEP_CONV_NORM
        PEP_CONV_USER

    ctypedef enum SlepcPEPStop "PEPStop":
        PEP_STOP_BASIC
        PEP_STOP_USER

    ctypedef enum SlepcPEPConvergedReason "PEPConvergedReason":
        PEP_CONVERGED_TOL
        PEP_CONVERGED_USER
        PEP_DIVERGED_ITS
        PEP_DIVERGED_BREAKDOWN
        PEP_DIVERGED_SYMMETRY_LOST
        PEP_CONVERGED_ITERATING

    ctypedef PetscErrorCode (*SlepcPEPCtxDel)(void*)
    ctypedef PetscErrorCode (*SlepcPEPStoppingFunction)(SlepcPEP,
                                             PetscInt,
                                             PetscInt,
                                             PetscInt,
                                             PetscInt,
                                             SlepcPEPConvergedReason*,
                                             void*) except PETSC_ERR_PYTHON
    ctypedef PetscErrorCode (*SlepcPEPMonitorFunction)(SlepcPEP,
                                            PetscInt,
                                            PetscInt,
                                            PetscScalar*,
                                            PetscScalar*,
                                            PetscReal*,
                                            PetscInt,
                                            void*) except PETSC_ERR_PYTHON

    PetscErrorCode PEPCreate(MPI_Comm,SlepcPEP*)
    PetscErrorCode PEPDestroy(SlepcPEP*)
    PetscErrorCode PEPReset(SlepcPEP)
    PetscErrorCode PEPView(SlepcPEP,PetscViewer)

    PetscErrorCode PEPSetType(SlepcPEP,SlepcPEPType)
    PetscErrorCode PEPGetType(SlepcPEP,SlepcPEPType*)
    PetscErrorCode PEPSetBasis(SlepcPEP,SlepcPEPBasis)
    PetscErrorCode PEPGetBasis(SlepcPEP,SlepcPEPBasis*)
    PetscErrorCode PEPSetProblemType(SlepcPEP,SlepcPEPProblemType)
    PetscErrorCode PEPGetProblemType(SlepcPEP,SlepcPEPProblemType*)
    PetscErrorCode PEPSetOperators(SlepcPEP,PetscInt,PetscMat*)
    PetscErrorCode PEPGetOperators(SlepcPEP,PetscInt,PetscMat*)
    PetscErrorCode PEPGetNumMatrices(SlepcPEP,PetscInt*)
    PetscErrorCode PEPSetOptionsPrefix(SlepcPEP,char*)
    PetscErrorCode PEPGetOptionsPrefix(SlepcPEP,char*[])
    PetscErrorCode PEPSetFromOptions(SlepcPEP)
    PetscErrorCode PEPAppendOptionsPrefix(SlepcPEP,char*)
    PetscErrorCode PEPSetUp(SlepcPEP)
    PetscErrorCode PEPSolve(SlepcPEP)

    PetscErrorCode PEPSetBV(SlepcPEP,SlepcBV)
    PetscErrorCode PEPGetBV(SlepcPEP,SlepcBV*)
    PetscErrorCode PEPSetDS(SlepcPEP,SlepcDS)
    PetscErrorCode PEPGetDS(SlepcPEP,SlepcDS*)
    PetscErrorCode PEPSetST(SlepcPEP,SlepcST)
    PetscErrorCode PEPGetST(SlepcPEP,SlepcST*)
    PetscErrorCode PEPSetRG(SlepcPEP,SlepcRG)
    PetscErrorCode PEPGetRG(SlepcPEP,SlepcRG*)

    PetscErrorCode PEPSetTrackAll(SlepcPEP,PetscBool)
    PetscErrorCode PEPGetTrackAll(SlepcPEP,PetscBool*)

    PetscErrorCode PEPSetTolerances(SlepcPEP,PetscReal,PetscInt)
    PetscErrorCode PEPGetTolerances(SlepcPEP,PetscReal*,PetscInt*)
    PetscErrorCode PEPSetDimensions(SlepcPEP,PetscInt,PetscInt,PetscInt)
    PetscErrorCode PEPGetDimensions(SlepcPEP,PetscInt*,PetscInt*,PetscInt*)
    PetscErrorCode PEPSetScale(SlepcPEP,SlepcPEPScale,PetscReal,PetscVec,PetscVec,PetscInt,PetscReal)
    PetscErrorCode PEPGetScale(SlepcPEP,SlepcPEPScale*,PetscReal*,PetscVec*,PetscVec*,PetscInt*,PetscReal*)

    PetscErrorCode PEPGetConverged(SlepcPEP,PetscInt*)
    PetscErrorCode PEPGetEigenpair(SlepcPEP,PetscInt,PetscScalar*,PetscScalar*,PetscVec,PetscVec)
    PetscErrorCode PEPComputeError(SlepcPEP,PetscInt,SlepcPEPErrorType,PetscReal*)
    PetscErrorCode PEPErrorView(SlepcPEP,SlepcPEPErrorType,PetscViewer)
    PetscErrorCode PEPValuesView(SlepcPEP,PetscViewer)
    PetscErrorCode PEPVectorsView(SlepcPEP,PetscViewer)
    PetscErrorCode PEPGetErrorEstimate(SlepcPEP,PetscInt,PetscReal*)

    PetscErrorCode PEPSetStoppingTestFunction(SlepcPEP,SlepcPEPStoppingFunction,void*,SlepcPEPCtxDel)
    PetscErrorCode PEPStoppingBasic(SlepcPEP,PetscInt,PetscInt,PetscInt,PetscInt,SlepcPEPConvergedReason*,void*) except PETSC_ERR_PYTHON

    PetscErrorCode PEPSetConvergenceTest(SlepcPEP,SlepcPEPConv)
    PetscErrorCode PEPGetConvergenceTest(SlepcPEP,SlepcPEPConv*)
    PetscErrorCode PEPSetRefine(SlepcPEP,SlepcPEPRefine,PetscInt,PetscReal,PetscInt,SlepcPEPRefineScheme)
    PetscErrorCode PEPGetRefine(SlepcPEP,SlepcPEPRefine*,PetscInt*,PetscReal*,PetscInt*,SlepcPEPRefineScheme*)
    PetscErrorCode PEPRefineGetKSP(SlepcPEP,PetscKSP*)
    PetscErrorCode PEPSetExtract(SlepcPEP,SlepcPEPExtract);
    PetscErrorCode PEPGetExtract(SlepcPEP,SlepcPEPExtract*)

    PetscErrorCode PEPMonitorSet(SlepcPEP,SlepcPEPMonitorFunction,void*,SlepcPEPCtxDel)
    PetscErrorCode PEPMonitorCancel(SlepcPEP)
    PetscErrorCode PEPGetIterationNumber(SlepcPEP,PetscInt*)

    PetscErrorCode PEPSetInitialSpace(SlepcPEP,PetscInt,PetscVec*)
    PetscErrorCode PEPSetWhichEigenpairs(SlepcPEP,SlepcPEPWhich)
    PetscErrorCode PEPGetWhichEigenpairs(SlepcPEP,SlepcPEPWhich*)
    PetscErrorCode PEPSetTarget(SlepcPEP,PetscScalar)
    PetscErrorCode PEPGetTarget(SlepcPEP,PetscScalar*)
    PetscErrorCode PEPSetInterval(SlepcPEP,PetscReal,PetscReal)
    PetscErrorCode PEPGetInterval(SlepcPEP,PetscReal*,PetscReal*)
    PetscErrorCode PEPGetConvergedReason(SlepcPEP,SlepcPEPConvergedReason*)

    PetscErrorCode PEPLinearSetLinearization(SlepcPEP,PetscReal,PetscReal)
    PetscErrorCode PEPLinearGetLinearization(SlepcPEP,PetscReal*,PetscReal*)
    PetscErrorCode PEPLinearSetExplicitMatrix(SlepcPEP,PetscBool)
    PetscErrorCode PEPLinearGetExplicitMatrix(SlepcPEP,PetscBool*)
    PetscErrorCode PEPLinearSetEPS(SlepcPEP,SlepcEPS)
    PetscErrorCode PEPLinearGetEPS(SlepcPEP,SlepcEPS*)

    PetscErrorCode PEPQArnoldiSetRestart(SlepcPEP,PetscReal)
    PetscErrorCode PEPQArnoldiGetRestart(SlepcPEP,PetscReal*)
    PetscErrorCode PEPQArnoldiSetLocking(SlepcPEP,PetscBool)
    PetscErrorCode PEPQArnoldiGetLocking(SlepcPEP,PetscBool*)

    PetscErrorCode PEPTOARSetRestart(SlepcPEP,PetscReal)
    PetscErrorCode PEPTOARGetRestart(SlepcPEP,PetscReal*)
    PetscErrorCode PEPTOARSetLocking(SlepcPEP,PetscBool)
    PetscErrorCode PEPTOARGetLocking(SlepcPEP,PetscBool*)

    PetscErrorCode PEPSTOARSetLinearization(SlepcPEP,PetscReal,PetscReal)
    PetscErrorCode PEPSTOARGetLinearization(SlepcPEP,PetscReal*,PetscReal*)
    PetscErrorCode PEPSTOARSetLocking(SlepcPEP,PetscBool)
    PetscErrorCode PEPSTOARGetLocking(SlepcPEP,PetscBool*)
    PetscErrorCode PEPSTOARSetDetectZeros(SlepcPEP,PetscBool)
    PetscErrorCode PEPSTOARGetDetectZeros(SlepcPEP,PetscBool*)
    PetscErrorCode PEPSTOARSetDimensions(SlepcPEP,PetscInt,PetscInt,PetscInt)
    PetscErrorCode PEPSTOARGetDimensions(SlepcPEP,PetscInt*,PetscInt*,PetscInt*)
    PetscErrorCode PEPSTOARGetInertias(SlepcPEP,PetscInt*,PetscReal**,PetscInt**)
    PetscErrorCode PEPSTOARSetCheckEigenvalueType(SlepcPEP,PetscBool)
    PetscErrorCode PEPSTOARGetCheckEigenvalueType(SlepcPEP,PetscBool*)

    ctypedef enum SlepcPEPJDProjection "PEPJDProjection":
        PEP_JD_PROJECTION_HARMONIC
        PEP_JD_PROJECTION_ORTHOGONAL

    PetscErrorCode PEPJDSetRestart(SlepcPEP,PetscReal)
    PetscErrorCode PEPJDGetRestart(SlepcPEP,PetscReal*)
    PetscErrorCode PEPJDSetFix(SlepcPEP,PetscReal)
    PetscErrorCode PEPJDGetFix(SlepcPEP,PetscReal*)
    PetscErrorCode PEPJDSetReusePreconditioner(SlepcPEP,PetscBool)
    PetscErrorCode PEPJDGetReusePreconditioner(SlepcPEP,PetscBool*)
    PetscErrorCode PEPJDSetMinimalityIndex(SlepcPEP,PetscInt)
    PetscErrorCode PEPJDGetMinimalityIndex(SlepcPEP,PetscInt*)
    PetscErrorCode PEPJDSetProjection(SlepcPEP,SlepcPEPJDProjection)
    PetscErrorCode PEPJDGetProjection(SlepcPEP,SlepcPEPJDProjection*)

    ctypedef enum SlepcPEPCISSExtraction "PEPCISSExtraction":
        PEP_CISS_EXTRACTION_RITZ
        PEP_CISS_EXTRACTION_HANKEL
        PEP_CISS_EXTRACTION_CAA

    PetscErrorCode PEPCISSSetExtraction(SlepcPEP,SlepcPEPCISSExtraction)
    PetscErrorCode PEPCISSGetExtraction(SlepcPEP,SlepcPEPCISSExtraction*)
    PetscErrorCode PEPCISSSetSizes(SlepcPEP,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscBool)
    PetscErrorCode PEPCISSGetSizes(SlepcPEP,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscBool*)
    PetscErrorCode PEPCISSSetThreshold(SlepcPEP,PetscReal,PetscReal)
    PetscErrorCode PEPCISSGetThreshold(SlepcPEP,PetscReal*,PetscReal*)
    PetscErrorCode PEPCISSSetRefinement(SlepcPEP,PetscInt,PetscInt)
    PetscErrorCode PEPCISSGetRefinement(SlepcPEP,PetscInt*,PetscInt*)
    PetscErrorCode PEPCISSGetKSPs(SlepcPEP,PetscInt*,PetscKSP**)

# -----------------------------------------------------------------------------

cdef inline PEP ref_PEP(SlepcPEP pep):
    cdef PEP ob = <PEP> PEP()
    ob.pep = pep
    CHKERR( PetscINCREF(ob.obj) )
    return ob

# -----------------------------------------------------------------------------

cdef PetscErrorCode PEP_Stopping(
    SlepcPEP                pep,
    PetscInt                its,
    PetscInt                max_it,
    PetscInt                nconv,
    PetscInt                nev,
    SlepcPEPConvergedReason *r,
    void                    *ctx,
    ) except PETSC_ERR_PYTHON with gil:
    cdef PEP Pep = ref_PEP(pep)
    (stopping, args, kargs) = Pep.get_attr('__stopping__')
    reason = stopping(Pep, toInt(its), toInt(max_it), toInt(nconv), toInt(nev), *args, **kargs)
    if   reason is None:  r[0] = PEP_CONVERGED_ITERATING
    elif reason is False: r[0] = PEP_CONVERGED_ITERATING
    elif reason is True:  r[0] = PEP_CONVERGED_USER
    else:                 r[0] = reason
    return PETSC_SUCCESS

# -----------------------------------------------------------------------------

cdef PetscErrorCode PEP_Monitor(
    SlepcPEP    pep,
    PetscInt    its,
    PetscInt    nconv,
    PetscScalar *eigr,
    PetscScalar *eigi,
    PetscReal   *errest,
    PetscInt    nest,
    void        *ctx,
    ) except PETSC_ERR_PYTHON with gil:
    cdef PEP Pep = ref_PEP(pep)
    cdef object monitorlist = Pep.get_attr('__monitor__')
    if monitorlist is None: return PETSC_SUCCESS
    cdef object eig = [toComplex(eigr[i], eigi[i]) for i in range(nest)]
    cdef object err = [toReal(errest[i]) for i in range(nest)]
    for (monitor, args, kargs) in monitorlist:
        monitor(Pep, toInt(its), toInt(nconv), eig, err, *args, **kargs)
    return PETSC_SUCCESS

# -----------------------------------------------------------------------------
