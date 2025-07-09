cdef extern from * nogil:

    ctypedef char* SlepcFNType "const char*"
    SlepcFNType FNCOMBINE
    SlepcFNType FNRATIONAL
    SlepcFNType FNEXP
    SlepcFNType FNLOG
    SlepcFNType FNPHI
    SlepcFNType FNSQRT
    SlepcFNType FNINVSQRT

    ctypedef enum SlepcFNCombineType "FNCombineType":
        FN_COMBINE_ADD
        FN_COMBINE_MULTIPLY
        FN_COMBINE_DIVIDE
        FN_COMBINE_COMPOSE

    ctypedef enum SlepcFNParallelType "FNParallelType":
        FN_PARALLEL_REDUNDANT
        FN_PARALLEL_SYNCHRONIZED

    PetscErrorCode FNCreate(MPI_Comm,SlepcFN*)
    PetscErrorCode FNView(SlepcFN,PetscViewer)
    PetscErrorCode FNDestroy(SlepcFN*)
    PetscErrorCode FNReset(SlepcFN)
    PetscErrorCode FNSetType(SlepcFN,SlepcFNType)
    PetscErrorCode FNGetType(SlepcFN,SlepcFNType*)

    PetscErrorCode FNSetOptionsPrefix(SlepcFN,char[])
    PetscErrorCode FNGetOptionsPrefix(SlepcFN,char*[])
    PetscErrorCode FNAppendOptionsPrefix(SlepcFN,char[])
    PetscErrorCode FNSetFromOptions(SlepcFN)
    PetscErrorCode FNDuplicate(SlepcFN,MPI_Comm,SlepcFN*)

    PetscErrorCode FNSetScale(SlepcFN,PetscScalar,PetscScalar)
    PetscErrorCode FNGetScale(SlepcFN,PetscScalar*,PetscScalar*)
    PetscErrorCode FNSetMethod(SlepcFN,PetscInt)
    PetscErrorCode FNGetMethod(SlepcFN,PetscInt*)
    PetscErrorCode FNSetParallel(SlepcFN,SlepcFNParallelType)
    PetscErrorCode FNGetParallel(SlepcFN,SlepcFNParallelType*)
    PetscErrorCode FNEvaluateFunction(SlepcFN,PetscScalar,PetscScalar*)
    PetscErrorCode FNEvaluateDerivative(SlepcFN,PetscScalar,PetscScalar*)
    PetscErrorCode FNEvaluateFunctionMat(SlepcFN,PetscMat,PetscMat)
    PetscErrorCode FNEvaluateFunctionMatVec(SlepcFN,PetscMat,PetscVec)

    PetscErrorCode FNRationalSetNumerator(SlepcFN,PetscInt,PetscScalar[])
    PetscErrorCode FNRationalGetNumerator(SlepcFN,PetscInt*,PetscScalar*[])
    PetscErrorCode FNRationalSetDenominator(SlepcFN,PetscInt,PetscScalar[])
    PetscErrorCode FNRationalGetDenominator(SlepcFN,PetscInt*,PetscScalar*[])

    PetscErrorCode FNCombineSetChildren(SlepcFN,SlepcFNCombineType,SlepcFN,SlepcFN)
    PetscErrorCode FNCombineGetChildren(SlepcFN,SlepcFNCombineType*,SlepcFN*,SlepcFN*)

    PetscErrorCode FNPhiSetIndex(SlepcFN,PetscInt)
    PetscErrorCode FNPhiGetIndex(SlepcFN,PetscInt*)

# --------------------------------------------------------------------

# unary operations

cdef FN fn_pos(FN self):
    cdef FN fn = type(self)()
    cdef MPI_Comm comm = def_Comm(None, PetscObjectComm(<PetscObject>self.fn))
    CHKERR(FNDuplicate(self.fn, comm, &fn.fn))
    return fn

cdef FN fn_neg(FN self):
    cdef PetscScalar alpha = 1, beta = 1
    cdef FN fn = <FN> fn_pos(self)
    CHKERR(FNGetScale(fn.fn, &alpha, &beta))
    CHKERR(FNSetScale(fn.fn, alpha, -beta))
    return fn

# inplace binary operations

cdef FN fn_iadd(FN self, other):
    cdef PetscScalar alpha = 1
    cdef FN comb = type(self)()
    cdef SlepcFN fn = NULL
    cdef MPI_Comm comm = def_Comm(None, PetscObjectComm(<PetscObject>self.fn))
    CHKERR(FNCreate(comm, &comb.fn))
    CHKERR(FNSetType(comb.fn, FNCOMBINE))
    CHKERR(FNSetFromOptions(comb.fn))
    if isinstance(other, FN):
        fn = (<FN>other).fn
        CHKERR(PetscObjectReference(<PetscObject>fn))
    else:
        alpha = asScalar(other)
        CHKERR(FNCreate(comm, &fn))
        CHKERR(FNSetType(fn, FNRATIONAL))
        CHKERR(FNSetFromOptions(fn))
        CHKERR(FNRationalSetNumerator(fn, 1, &alpha));
    CHKERR(FNCombineSetChildren(comb.fn, FN_COMBINE_ADD, self.fn, fn))
    CHKERR(FNDestroy(&fn))
    return comb

cdef FN fn_isub(FN self, other):
    cdef PetscScalar alpha = 1, beta = 1
    cdef FN comb = type(self)()
    cdef SlepcFN fn = NULL
    cdef MPI_Comm comm = def_Comm(None, PetscObjectComm(<PetscObject>self.fn))
    CHKERR(FNCreate(comm, &comb.fn))
    CHKERR(FNSetType(comb.fn, FNCOMBINE))
    CHKERR(FNSetFromOptions(comb.fn))
    if isinstance(other, FN):
        CHKERR(FNDuplicate((<FN>other).fn, comm, &fn))
        CHKERR(FNGetScale(fn, &alpha, &beta))
        CHKERR(FNSetScale(fn, alpha, -beta))
    else:
        alpha = -asScalar(other)
        CHKERR(FNCreate(comm, &fn))
        CHKERR(FNSetType(fn, FNRATIONAL))
        CHKERR(FNSetFromOptions(fn))
        CHKERR(FNRationalSetNumerator(fn, 1, &alpha));
    CHKERR(FNCombineSetChildren(comb.fn, FN_COMBINE_ADD, self.fn, fn))
    CHKERR(FNDestroy(&fn))
    return comb

cdef FN fn_imul(FN self, other):
    cdef PetscScalar alpha = 1, beta = 1, gamma = 1
    cdef FN comb
    cdef MPI_Comm comm = def_Comm(None, PetscObjectComm(<PetscObject>self.fn))
    if isinstance(other, FN):
        comb = type(self)()
        CHKERR(FNCreate(comm, &comb.fn))
        CHKERR(FNSetType(comb.fn, FNCOMBINE))
        CHKERR(FNSetFromOptions(comb.fn))
        CHKERR(FNCombineSetChildren(comb.fn, FN_COMBINE_MULTIPLY, self.fn, (<FN>other).fn))
        return comb
    else:
        gamma = asScalar(other)
        CHKERR(FNGetScale(self.fn, &alpha, &beta))
        CHKERR(FNSetScale(self.fn, alpha, beta*gamma))
        return self

cdef FN fn_idiv(FN self, other):
    cdef PetscScalar alpha = 1, beta = 1, gamma = 1
    cdef FN comb
    cdef MPI_Comm comm = def_Comm(None, PetscObjectComm(<PetscObject>self.fn))
    if isinstance(other, FN):
        comb = type(self)()
        CHKERR(FNCreate(comm, &comb.fn))
        CHKERR(FNSetType(comb.fn, FNCOMBINE))
        CHKERR(FNSetFromOptions(comb.fn))
        CHKERR(FNCombineSetChildren(comb.fn, FN_COMBINE_DIVIDE, self.fn, (<FN>other).fn))
        return comb
    else:
        gamma = asScalar(other)
        CHKERR(FNGetScale(self.fn, &alpha, &beta))
        CHKERR(FNSetScale(self.fn, alpha, beta/gamma))
        return self

# binary operations

cdef FN fn_add(FN self, other):
    return fn_iadd(fn_pos(self), other)

cdef FN fn_sub(FN self, other):
    return fn_isub(fn_pos(self), other)

cdef FN fn_mul(FN self, other):
    return fn_imul(fn_pos(self), other)

cdef FN fn_div(FN self, other):
    return fn_idiv(fn_pos(self), other)

# reflected binary operations

cdef FN fn_radd(FN self, other):
    return fn_add(self, other)

cdef FN fn_rsub(FN self, other):
    cdef FN fn = <FN> fn_sub(self, other)
    CHKERR(FNSetScale(fn.fn, 1, -1))
    return fn

cdef FN fn_rmul(FN self, other):
    return fn_mul(self, other)

cdef FN fn_rdiv(FN self, other):
    <void>self; <void>other # unused
    return NotImplemented

# composition

cdef FN fn_matmul(FN self, other):
    cdef FN fn, comb
    cdef MPI_Comm comm = def_Comm(None, PetscObjectComm(<PetscObject>self.fn))
    if isinstance(other, FN):
        comb = type(self)()
        CHKERR(FNCreate(comm, &comb.fn))
        CHKERR(FNSetType(comb.fn, FNCOMBINE))
        CHKERR(FNSetFromOptions(comb.fn))
        CHKERR(FNCombineSetChildren(comb.fn, FN_COMBINE_COMPOSE, (<FN>other).fn, self.fn))
        return comb
    else:
        return NotImplemented

