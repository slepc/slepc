cdef extern from * nogil:

    ctypedef char* SlepcSVDType "const char*"
    SlepcSVDType SVDCROSS
    SlepcSVDType SVDCYCLIC
    SlepcSVDType SVDLAPACK
    SlepcSVDType SVDLANCZOS
    SlepcSVDType SVDTRLANCZOS
    SlepcSVDType SVDRANDOMIZED
    SlepcSVDType SVDSCALAPACK
    SlepcSVDType SVDKSVD
    SlepcSVDType SVDELEMENTAL
    SlepcSVDType SVDPRIMME

    ctypedef enum SlepcSVDProblemType "SVDProblemType":
        SVD_STANDARD
        SVD_GENERALIZED
        SVD_HYPERBOLIC

    ctypedef enum SlepcSVDWhich "SVDWhich":
        SVD_LARGEST
        SVD_SMALLEST

    ctypedef enum SlepcSVDErrorType "SVDErrorType":
        SVD_ERROR_ABSOLUTE
        SVD_ERROR_RELATIVE
        SVD_ERROR_NORM

    ctypedef enum SlepcSVDConv "SVDConv":
        SVD_CONV_ABS
        SVD_CONV_REL
        SVD_CONV_NORM
        SVD_CONV_MAXIT
        SVD_CONV_USER

    ctypedef enum SlepcSVDStop "SVDStop":
        SVD_STOP_BASIC
        SVD_STOP_USER
        SVD_STOP_THRESHOLD

    ctypedef enum SlepcSVDConvergedReason "SVDConvergedReason":
        SVD_CONVERGED_TOL
        SVD_CONVERGED_USER
        SVD_CONVERGED_MAXIT
        SVD_DIVERGED_ITS
        SVD_DIVERGED_BREAKDOWN
        SVD_DIVERGED_SYMMETRY_LOST
        SVD_CONVERGED_ITERATING

    ctypedef PetscErrorCode (*SlepcSVDCtxDel)(void*)
    ctypedef PetscErrorCode (*SlepcSVDStoppingFunction)(SlepcSVD,
                                             PetscInt,
                                             PetscInt,
                                             PetscInt,
                                             PetscInt,
                                             SlepcSVDConvergedReason*,
                                             void*) except PETSC_ERR_PYTHON
    ctypedef PetscErrorCode (*SlepcSVDMonitorFunction)(SlepcSVD,
                                            PetscInt,
                                            PetscInt,
                                            PetscReal*,
                                            PetscReal*,
                                            PetscInt,
                                            void*) except PETSC_ERR_PYTHON

    PetscErrorCode SVDCreate(MPI_Comm,SlepcSVD*)
    PetscErrorCode SVDView(SlepcSVD,PetscViewer)
    PetscErrorCode SVDDestroy(SlepcSVD*)
    PetscErrorCode SVDReset(SlepcSVD)
    PetscErrorCode SVDSetType(SlepcSVD,SlepcSVDType)
    PetscErrorCode SVDGetType(SlepcSVD,SlepcSVDType*)
    PetscErrorCode SVDSetOptionsPrefix(SlepcSVD,char[])
    PetscErrorCode SVDAppendOptionsPrefix(SlepcSVD,char[])
    PetscErrorCode SVDGetOptionsPrefix(SlepcSVD,char*[])
    PetscErrorCode SVDSetFromOptions(SlepcSVD)

    PetscErrorCode SVDSetProblemType(SlepcSVD,SlepcSVDProblemType)
    PetscErrorCode SVDGetProblemType(SlepcSVD,SlepcSVDProblemType*)
    PetscErrorCode SVDIsGeneralized(SlepcSVD,PetscBool*)
    PetscErrorCode SVDIsHyperbolic(SlepcSVD,PetscBool*)
    PetscErrorCode SVDSetBV(SlepcSVD,SlepcBV,SlepcBV)
    PetscErrorCode SVDGetBV(SlepcSVD,SlepcBV*,SlepcBV*)
    PetscErrorCode SVDSetDS(SlepcSVD,SlepcDS)
    PetscErrorCode SVDGetDS(SlepcSVD,SlepcDS*)

    PetscErrorCode SVDSetOperators(SlepcSVD,PetscMat,PetscMat)
    PetscErrorCode SVDGetOperators(SlepcSVD,PetscMat*,PetscMat*)
    PetscErrorCode SVDSetSignature(SlepcSVD,PetscVec)
    PetscErrorCode SVDGetSignature(SlepcSVD,PetscVec)

    PetscErrorCode SVDSetInitialSpaces(SlepcSVD,PetscInt,PetscVec*,PetscInt,PetscVec*)

    PetscErrorCode SVDSetImplicitTranspose(SlepcSVD,PetscBool)
    PetscErrorCode SVDGetImplicitTranspose(SlepcSVD,PetscBool*)
    PetscErrorCode SVDSetDimensions(SlepcSVD,PetscInt,PetscInt,PetscInt)
    PetscErrorCode SVDGetDimensions(SlepcSVD,PetscInt*,PetscInt*,PetscInt*)
    PetscErrorCode SVDSetTolerances(SlepcSVD,PetscReal,PetscInt)
    PetscErrorCode SVDGetTolerances(SlepcSVD,PetscReal*,PetscInt*)
    PetscErrorCode SVDSetWhichSingularTriplets(SlepcSVD,SlepcSVDWhich)
    PetscErrorCode SVDGetWhichSingularTriplets(SlepcSVD,SlepcSVDWhich*)
    PetscErrorCode SVDSetThreshold(SlepcSVD,PetscReal,PetscBool)
    PetscErrorCode SVDGetThreshold(SlepcSVD,PetscReal*,PetscBool*)
    PetscErrorCode SVDSetConvergenceTest(SlepcSVD,SlepcSVDConv)
    PetscErrorCode SVDGetConvergenceTest(SlepcSVD,SlepcSVDConv*)

    PetscErrorCode SVDMonitorSet(SlepcSVD,SlepcSVDMonitorFunction,void*,SlepcSVDCtxDel)
    PetscErrorCode SVDMonitorCancel(SlepcSVD)

    PetscErrorCode SVDSetStoppingTestFunction(SlepcSVD,SlepcSVDStoppingFunction,void*,SlepcSVDCtxDel)
    PetscErrorCode SVDStoppingBasic(SlepcSVD,PetscInt,PetscInt,PetscInt,PetscInt,SlepcSVDConvergedReason*,void*) except PETSC_ERR_PYTHON

    PetscErrorCode SVDSetTrackAll(SlepcSVD,PetscBool)
    PetscErrorCode SVDGetTrackAll(SlepcSVD,PetscBool*)

    PetscErrorCode SVDSetUp(SlepcSVD)
    PetscErrorCode SVDSolve(SlepcSVD)
    PetscErrorCode SVDGetIterationNumber(SlepcSVD,PetscInt*)
    PetscErrorCode SVDGetConvergedReason(SlepcSVD,SlepcSVDConvergedReason*)
    PetscErrorCode SVDGetConverged(SlepcSVD,PetscInt*)
    PetscErrorCode SVDGetSingularTriplet(SlepcSVD,PetscInt,PetscReal*,PetscVec,PetscVec)
    PetscErrorCode SVDComputeError(SlepcSVD,PetscInt,SlepcSVDErrorType,PetscReal*)
    PetscErrorCode SVDErrorView(SlepcSVD,SlepcSVDErrorType,PetscViewer)
    PetscErrorCode SVDValuesView(SlepcSVD,PetscViewer)
    PetscErrorCode SVDVectorsView(SlepcSVD,PetscViewer)

    PetscErrorCode SVDCrossSetExplicitMatrix(SlepcSVD,PetscBool)
    PetscErrorCode SVDCrossGetExplicitMatrix(SlepcSVD,PetscBool*)
    PetscErrorCode SVDCrossSetEPS(SlepcSVD,SlepcEPS)
    PetscErrorCode SVDCrossGetEPS(SlepcSVD,SlepcEPS*)

    PetscErrorCode SVDCyclicSetExplicitMatrix(SlepcSVD,PetscBool)
    PetscErrorCode SVDCyclicGetExplicitMatrix(SlepcSVD,PetscBool*)
    PetscErrorCode SVDCyclicSetEPS(SlepcSVD,SlepcEPS)
    PetscErrorCode SVDCyclicGetEPS(SlepcSVD,SlepcEPS*)

    PetscErrorCode SVDLanczosSetOneSide(SlepcSVD,PetscBool)
    PetscErrorCode SVDLanczosGetOneSide(SlepcSVD,PetscBool*)

    PetscErrorCode SVDTRLanczosSetOneSide(SlepcSVD,PetscBool)
    PetscErrorCode SVDTRLanczosGetOneSide(SlepcSVD,PetscBool*)
    PetscErrorCode SVDTRLanczosSetGBidiag(SlepcSVD,SlepcSVDTRLanczosGBidiag)
    PetscErrorCode SVDTRLanczosGetGBidiag(SlepcSVD,SlepcSVDTRLanczosGBidiag*)
    PetscErrorCode SVDTRLanczosSetKSP(SlepcSVD,PetscKSP)
    PetscErrorCode SVDTRLanczosGetKSP(SlepcSVD,PetscKSP*)
    PetscErrorCode SVDTRLanczosSetRestart(SlepcSVD,PetscReal)
    PetscErrorCode SVDTRLanczosGetRestart(SlepcSVD,PetscReal*)
    PetscErrorCode SVDTRLanczosSetLocking(SlepcSVD,PetscBool)
    PetscErrorCode SVDTRLanczosGetLocking(SlepcSVD,PetscBool*)
    PetscErrorCode SVDTRLanczosSetExplicitMatrix(SlepcSVD,PetscBool)
    PetscErrorCode SVDTRLanczosGetExplicitMatrix(SlepcSVD,PetscBool*)

    ctypedef enum SlepcSVDTRLanczosGBidiag "SVDTRLanczosGBidiag":
        SVD_TRLANCZOS_GBIDIAG_SINGLE
        SVD_TRLANCZOS_GBIDIAG_UPPER
        SVD_TRLANCZOS_GBIDIAG_LOWER

# -----------------------------------------------------------------------------

cdef inline SVD ref_SVD(SlepcSVD svd):
    cdef SVD ob = <SVD> SVD()
    ob.svd = svd
    CHKERR( PetscINCREF(ob.obj) )
    return ob

# -----------------------------------------------------------------------------

cdef PetscErrorCode SVD_Stopping(
    SlepcSVD                svd,
    PetscInt                its,
    PetscInt                max_it,
    PetscInt                nconv,
    PetscInt                nev,
    SlepcSVDConvergedReason *r,
    void                    *ctx,
    ) except PETSC_ERR_PYTHON with gil:
    cdef SVD Svd = ref_SVD(svd)
    (stopping, args, kargs) = Svd.get_attr('__stopping__')
    reason = stopping(Svd, toInt(its), toInt(max_it), toInt(nconv), toInt(nev), *args, **kargs)
    if   reason is None:  r[0] = SVD_CONVERGED_ITERATING
    elif reason is False: r[0] = SVD_CONVERGED_ITERATING
    elif reason is True:  r[0] = SVD_CONVERGED_USER
    else:                 r[0] = reason

# -----------------------------------------------------------------------------

cdef PetscErrorCode SVD_Monitor(
    SlepcSVD    svd,
    PetscInt    its,
    PetscInt    nconv,
    PetscReal   *sigma,
    PetscReal   *errest,
    PetscInt    nest,
    void        *ctx,
    ) except PETSC_ERR_PYTHON with gil:
    cdef SVD Svd = ref_SVD(svd)
    cdef object monitorlist = Svd.get_attr('__monitor__')
    if monitorlist is None: return PETSC_SUCCESS
    cdef object eig = [toReal(sigma[i]) for i in range(nest)]
    cdef object err = [toReal(errest[i]) for i in range(nest)]
    for (monitor, args, kargs) in monitorlist:
        monitor(Svd, toInt(its), toInt(nconv), eig, err, *args, **kargs)
    return PETSC_SUCCESS

# -----------------------------------------------------------------------------
