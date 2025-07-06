cdef extern from * nogil:

    ctypedef char* SlepcBVType "const char*"
    SlepcBVType BVMAT
    SlepcBVType BVSVEC
    SlepcBVType BVVECS
    SlepcBVType BVCONTIGUOUS
    SlepcBVType BVTENSOR

    ctypedef enum SlepcBVOrthogType "BVOrthogType":
        BV_ORTHOG_CGS
        BV_ORTHOG_MGS

    ctypedef enum SlepcBVOrthogRefineType "BVOrthogRefineType":
        BV_ORTHOG_REFINE_IFNEEDED
        BV_ORTHOG_REFINE_NEVER
        BV_ORTHOG_REFINE_ALWAYS

    ctypedef enum SlepcBVOrthogBlockType "BVOrthogBlockType":
        BV_ORTHOG_BLOCK_GS
        BV_ORTHOG_BLOCK_CHOL
        BV_ORTHOG_BLOCK_TSQR
        BV_ORTHOG_BLOCK_TSQRCHOL
        BV_ORTHOG_BLOCK_SVQB

    ctypedef enum SlepcBVMatMultType "BVMatMultType":
        BV_MATMULT_VECS
        BV_MATMULT_MAT

    ctypedef enum SlepcBVSVDMethod "BVSVDMethod":
        BV_SVD_METHOD_REFINE
        BV_SVD_METHOD_QR
        BV_SVD_METHOD_QR_CAA

    PetscErrorCode BVCreate(MPI_Comm,SlepcBV*)
    PetscErrorCode BVCreateMat(SlepcBV,PetscMat*)
    PetscErrorCode BVDuplicate(SlepcBV,SlepcBV*)
    PetscErrorCode BVDuplicateResize(SlepcBV,PetscInt,SlepcBV*)
    PetscErrorCode BVCopy(SlepcBV,SlepcBV)
    PetscErrorCode BVView(SlepcBV,PetscViewer)
    PetscErrorCode BVDestroy(SlepcBV*)
    PetscErrorCode BVSetType(SlepcBV,SlepcBVType)
    PetscErrorCode BVGetType(SlepcBV,SlepcBVType*)
    PetscErrorCode BVSetSizes(SlepcBV,PetscInt,PetscInt,PetscInt)
    PetscErrorCode BVSetSizesFromVec(SlepcBV,PetscVec,PetscInt)
    PetscErrorCode BVGetSizes(SlepcBV,PetscInt*,PetscInt*,PetscInt*)
    PetscErrorCode BVResize(SlepcBV,PetscInt,PetscBool)
    PetscErrorCode BVSetLeadingDimension(SlepcBV,PetscInt)
    PetscErrorCode BVGetLeadingDimension(SlepcBV,PetscInt*)

    PetscErrorCode BVSetOptionsPrefix(SlepcBV,char[])
    PetscErrorCode BVGetOptionsPrefix(SlepcBV,char*[])
    PetscErrorCode BVAppendOptionsPrefix(SlepcBV,char[])
    PetscErrorCode BVSetFromOptions(SlepcBV)

    PetscErrorCode BVSetOrthogonalization(SlepcBV,SlepcBVOrthogType,SlepcBVOrthogRefineType,PetscReal,SlepcBVOrthogBlockType)
    PetscErrorCode BVGetOrthogonalization(SlepcBV,SlepcBVOrthogType*,SlepcBVOrthogRefineType*,PetscReal*,SlepcBVOrthogBlockType*)
    PetscErrorCode BVSetMatMultMethod(SlepcBV,SlepcBVMatMultType)
    PetscErrorCode BVGetMatMultMethod(SlepcBV,SlepcBVMatMultType*)

    PetscErrorCode BVSetRandom(SlepcBV)
    PetscErrorCode BVSetRandomNormal(SlepcBV)
    PetscErrorCode BVSetRandomSign(SlepcBV)
    PetscErrorCode BVSetRandomColumn(SlepcBV,PetscInt)
    PetscErrorCode BVSetRandomCond(SlepcBV,PetscReal)
    PetscErrorCode BVSetRandomContext(SlepcBV,PetscRandom)
    PetscErrorCode BVGetRandomContext(SlepcBV,PetscRandom*)

    PetscErrorCode BVSetMatrix(SlepcBV,PetscMat,PetscBool)
    PetscErrorCode BVGetMatrix(SlepcBV,PetscMat*,PetscBool*)
    PetscErrorCode BVApplyMatrix(SlepcBV,PetscVec,PetscVec)

    PetscErrorCode BVSetActiveColumns(SlepcBV,PetscInt,PetscInt)
    PetscErrorCode BVGetActiveColumns(SlepcBV,PetscInt*,PetscInt*)
    PetscErrorCode BVSetDefiniteTolerance(SlepcBV,PetscReal)
    PetscErrorCode BVGetDefiniteTolerance(SlepcBV,PetscReal*)

    PetscErrorCode BVCreateVec(SlepcBV,PetscVec*)
    PetscErrorCode BVSetVecType(SlepcBV,PetscVecType)
    PetscErrorCode BVGetVecType(SlepcBV,PetscVecType*)
    PetscErrorCode BVInsertVec(SlepcBV,PetscInt,PetscVec)
    PetscErrorCode BVInsertVecs(SlepcBV,PetscInt,PetscInt*,PetscVec*,PetscBool)
    PetscErrorCode BVInsertConstraints(SlepcBV,PetscInt*,PetscVec*)
    PetscErrorCode BVSetNumConstraints(SlepcBV,PetscInt)
    PetscErrorCode BVGetNumConstraints(SlepcBV,PetscInt*)
    PetscErrorCode BVGetColumn(SlepcBV,PetscInt,PetscVec*)
    PetscErrorCode BVRestoreColumn(SlepcBV,PetscInt,PetscVec*)
    PetscErrorCode BVCopyVec(SlepcBV,PetscInt,PetscVec)
    PetscErrorCode BVCopyColumn(SlepcBV,PetscInt,PetscInt)

    PetscErrorCode BVDot(SlepcBV,SlepcBV,PetscMat)
    PetscErrorCode BVDotVec(SlepcBV,PetscVec,PetscScalar*)
    PetscErrorCode BVDotColumn(SlepcBV,PetscInt,PetscScalar*)

    PetscErrorCode BVMatProject(SlepcBV,PetscMat,SlepcBV,PetscMat)
    PetscErrorCode BVMatMult(SlepcBV,PetscMat,SlepcBV)
    PetscErrorCode BVMatMultHermitianTranspose(SlepcBV,PetscMat,SlepcBV)
    PetscErrorCode BVMatMultColumn(SlepcBV,PetscMat,PetscInt)
    PetscErrorCode BVMatMultTransposeColumn(SlepcBV,PetscMat,PetscInt)
    PetscErrorCode BVMatMultHermitianTransposeColumn(SlepcBV,PetscMat,PetscInt)

    PetscErrorCode BVMult(SlepcBV,PetscScalar,PetscScalar,SlepcBV,PetscMat)
    PetscErrorCode BVMultColumn(SlepcBV,PetscScalar,PetscScalar,PetscInt,PetscScalar*)
    PetscErrorCode BVMultInPlace(SlepcBV,PetscMat,PetscInt,PetscInt)
    PetscErrorCode BVMultVec(SlepcBV,PetscScalar,PetscScalar,PetscVec,PetscScalar*)

    PetscErrorCode BVScaleColumn(SlepcBV,PetscInt,PetscScalar)
    PetscErrorCode BVScale(SlepcBV,PetscScalar)

    PetscErrorCode BVNormColumn(SlepcBV,PetscInt,PetscNormType,PetscReal*)
    PetscErrorCode BVNorm(SlepcBV,PetscNormType,PetscReal*)

    PetscErrorCode BVOrthogonalizeVec(SlepcBV,PetscVec,PetscScalar*,PetscReal*,PetscBool*)
    PetscErrorCode BVOrthogonalizeColumn(SlepcBV,PetscInt,PetscScalar*,PetscReal*,PetscBool*)
    PetscErrorCode BVOrthonormalizeColumn(SlepcBV,PetscInt,PetscBool,PetscReal*,PetscBool*)
    PetscErrorCode BVOrthogonalize(SlepcBV,PetscMat)

    PetscErrorCode BVCreateFromMat(PetscMat,SlepcBV*)
    PetscErrorCode BVGetMat(SlepcBV,PetscMat*)
    PetscErrorCode BVRestoreMat(SlepcBV,PetscMat*)

cdef inline PetscErrorCode BV_Sizes(
    object size,
    PetscInt *_n,
    PetscInt *_N,
    ) except PETSC_ERR_PYTHON:
    # unpack and get local and global sizes
    cdef PetscInt n=PETSC_DECIDE, N=PETSC_DECIDE
    cdef object on, oN
    try:
        on, oN = size
    except (TypeError, ValueError):
        on = None; oN = size
    if on is not None: n = asInt(on)
    if oN is not None: N = asInt(oN)
    if n==PETSC_DECIDE and N==PETSC_DECIDE: raise ValueError(
        "local and global sizes cannot be both 'DECIDE'")
    # return result to the caller
    if _n != NULL: _n[0] = n
    if _N != NULL: _N[0] = N
    return PETSC_SUCCESS

# --------------------------------------------------------------------

# unary operations

cdef BV bv_pos(BV self):
    cdef BV bv = type(self)()
    CHKERR(BVDuplicate(self.bv, &bv.bv))
    CHKERR(BVCopy(self.bv, bv.bv))
    return bv

cdef BV bv_neg(BV self):
    cdef BV bv = <BV> bv_pos(self)
    CHKERR(BVScale(bv.bv, -1))
    return bv

# inplace binary operations

cdef BV bv_iadd(BV self, other):
    cdef PetscScalar alpha = 1, beta = 1
    cdef BV bv
    if isinstance(other, BV):
        alpha = 1; bv = other
        CHKERR(BVMult(self.bv, alpha, beta, bv.bv, NULL))
        return self
    elif isinstance(other, (tuple, list)):
        other, bv = other
        alpha = asScalar(other)
        CHKERR(BVMult(self.bv, alpha, beta, bv.bv, NULL))
        return self
    return NotImplemented

cdef BV bv_isub(BV self, other):
    cdef PetscScalar alpha = 1, beta = 1
    cdef BV bv
    if isinstance(other, BV):
        alpha = 1; bv = other
        CHKERR(BVMult(self.bv, -alpha, beta, bv.bv, NULL))
        return self
    elif isinstance(other, (tuple, list)):
        other, bv = other
        alpha = asScalar(other)
        CHKERR(BVMult(self.bv, -alpha, beta, bv.bv, NULL))
        return self
    return NotImplemented

cdef BV bv_imul(BV self, other):
    self.scale(other)
    return self

cdef BV bv_idiv(BV self, other):
    other = 1/other
    self.scale(other)
    return self

# binary operations

cdef BV bv_add(BV self, other):
    return bv_iadd(bv_pos(self), other)

cdef BV bv_sub(BV self, other):
    return bv_isub(bv_pos(self), other)

cdef BV bv_mul(BV self, other):
    return bv_imul(bv_pos(self), other)

cdef BV bv_div(BV self, other):
    return bv_idiv(bv_pos(self), other)

# reflected binary operations

cdef BV bv_radd(BV self, other):
    return bv_add(self, other)

cdef BV bv_rsub(BV self, other):
    cdef BV bv = <BV> bv_sub(self, other)
    CHKERR(BVScale(bv.bv, -1))
    return bv

cdef BV bv_rmul(BV self, other):
    return bv_mul(self, other)

cdef BV bv_rdiv(BV self, other):
    <void>self; <void>other # unused
    return NotImplemented

