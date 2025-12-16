# -----------------------------------------------------------------------------

class FNType(object):
    """
    FN type.

    - `COMBINE`: A math function defined by combining two functions.
    - `RATIONAL`: A rational function :math:`f(x)=p(x)/q(x)`.
    - `EXP`: The exponential function :math:`f(x)=e^x`.
    - `LOG`: The logarithm function :math:`f(x)=\log{x}`.
    - `PHI`: One of the Phi_k functions with index k.
    - `SQRT`: The square root function :math:`f(x)=\sqrt{x}`.
    - `INVSQRT`: The inverse square root function.

    See Also
    --------
    slepc.FNType
    """
    COMBINE  = S_(FNCOMBINE)
    RATIONAL = S_(FNRATIONAL)
    EXP      = S_(FNEXP)
    LOG      = S_(FNLOG)
    PHI      = S_(FNPHI)
    SQRT     = S_(FNSQRT)
    INVSQRT  = S_(FNINVSQRT)

class FNCombineType(object):
    """
    FN type of combination of child functions.

    - `ADD`:       Addition       :math:`f(x) = f_1(x)+f_2(x)`
    - `MULTIPLY`:  Multiplication :math:`f(x) = f_1(x)f_2(x)`
    - `DIVIDE`:    Division       :math:`f(x) = f_1(x)/f_2(x)`
    - `COMPOSE`:   Composition    :math:`f(x) = f_2(f_1(x))`

    See Also
    --------
    slepc.FNCombineType
    """
    ADD      = FN_COMBINE_ADD
    MULTIPLY = FN_COMBINE_MULTIPLY
    DIVIDE   = FN_COMBINE_DIVIDE
    COMPOSE  = FN_COMBINE_COMPOSE

class FNParallelType(object):
    """
    FN parallel types.

    - `REDUNDANT`:    Every process performs the computation redundantly.
    - `SYNCHRONIZED`: The first process sends the result to the rest.

    See Also
    --------
    slepc.FNParallelType
    """
    REDUNDANT    = FN_PARALLEL_REDUNDANT
    SYNCHRONIZED = FN_PARALLEL_SYNCHRONIZED

# -----------------------------------------------------------------------------

cdef class FN(Object):

    """
    Mathematical Function.

    The `FN` package provides the functionality to represent a simple
    mathematical function such as an exponential, a polynomial or a rational
    function. This is used as a building block for defining the function
    associated to the nonlinear eigenproblem, as well as for specifying which
    function to use when computing the action of a matrix function on a vector.
    """

    Type         = FNType
    CombineType  = FNCombineType
    ParallelType = FNParallelType

    def __cinit__(self):
        self.obj = <PetscObject*> &self.fn
        self.fn = NULL

    # unary operations

    def __pos__(self):
        return fn_pos(self)

    def __neg__(self):
        return fn_neg(self)

    # inplace binary operations

    def __iadd__(self, other):
        return fn_iadd(self, other)

    def __isub__(self, other):
        return fn_isub(self, other)

    def __imul__(self, other):
        return fn_imul(self, other)

    def __idiv__(self, other):
        return fn_idiv(self, other)

    def __itruediv__(self, other):
        return fn_idiv(self, other)

    # binary operations

    def __add__(self, other):
        return fn_add(self, other)

    def __radd__(self, other):
        return fn_radd(self, other)

    def __sub__(self, other):
        return fn_sub(self, other)

    def __rsub__(self, other):
        return fn_rsub(self, other)

    def __mul__(self, other):
        return fn_mul(self, other)

    def __rmul__(self, other):
        return fn_rmul(self, other)

    def __div__(self, other):
        return fn_div(self, other)

    def __rdiv__(self, other):
        return fn_rdiv(self, other)

    def __truediv__(self, other):
        return fn_div(self, other)

    def __rtruediv__(self, other):
        return fn_rdiv(self, other)

    def __matmul__(self, other):
        return fn_matmul(self, other)

    def __call__(self, arg):
        if isinstance(arg, Mat):
            return self.evaluateFunctionMat(arg)
        else:
            return self.evaluateFunction(arg)

    #

    def view(self, Viewer viewer=None) -> None:
        """
        Print the FN data structure.

        Collective.

        Parameters
        ----------
        viewer
            Visualization context; if not provided, the standard
            output is used.

        See Also
        --------
        slepc.FNView
        """
        cdef PetscViewer vwr = def_Viewer(viewer)
        CHKERR( FNView(self.fn, vwr) )

    def destroy(self) -> Self:
        """
        Destroy the FN object.

        Collective.

        See Also
        --------
        slepc.FNDestroy
        """
        CHKERR( FNDestroy(&self.fn) )
        self.fn = NULL
        return self

    def create(self, comm: Comm | None = None) -> Self:
        """
        Create the FN object.

        Collective.

        Parameters
        ----------
        comm
            MPI communicator; if not provided, it defaults to all processes.

        See Also
        --------
        slepc.FNCreate
        """
        cdef MPI_Comm ccomm = def_Comm(comm, SLEPC_COMM_DEFAULT())
        cdef SlepcFN newfn = NULL
        CHKERR( FNCreate(ccomm, &newfn) )
        CHKERR( SlepcCLEAR(self.obj) ); self.fn = newfn
        return self

    def setType(self, fn_type: Type | str) -> None:
        """
        Set the type for the FN object.

        Logically collective.

        Parameters
        ----------
        fn_type
            The math function type to be used.

        See Also
        --------
        getType, slepc.FNSetType
        """
        cdef SlepcFNType cval = NULL
        fn_type = str2bytes(fn_type, &cval)
        CHKERR( FNSetType(self.fn, cval) )

    def getType(self) -> str:
        """
        Get the FN type of this object.

        Not collective.

        Returns
        -------
        str
            The math function type currently being used.

        See Also
        --------
        setType, slepc.FNGetType
        """
        cdef SlepcFNType fn_type = NULL
        CHKERR( FNGetType(self.fn, &fn_type) )
        return bytes2str(fn_type)

    def setOptionsPrefix(self, prefix: str | None = None) -> None:
        """
        Set the prefix used for searching for all FN options in the database.

        Logically collective.

        Parameters
        ----------
        prefix
            The prefix string to prepend to all FN option requests.

        Notes
        -----
        A hyphen (``-``) must NOT be given at the beginning of the
        prefix name.  The first character of all runtime options is
        AUTOMATICALLY the hyphen.

        See Also
        --------
        appendOptionsPrefix, getOptionsPrefix, slepc.FNGetOptionsPrefix
        """
        cdef const char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( FNSetOptionsPrefix(self.fn, cval) )

    def appendOptionsPrefix(self, prefix: str | None = None) -> None:
        """
        Append to the prefix used for searching for all FN options in the database.

        Logically collective.

        Parameters
        ----------
        prefix
            The prefix string to prepend to all FN option requests.

        See Also
        --------
        setOptionsPrefix, getOptionsPrefix, slepc.FNAppendOptionsPrefix
        """
        cdef const char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( FNAppendOptionsPrefix(self.fn, cval) )

    def getOptionsPrefix(self) -> str:
        """
        Get the prefix used for searching for all FN options in the database.

        Not collective.

        Returns
        -------
        str
            The prefix string set for this FN object.

        See Also
        --------
        setOptionsPrefix, appendOptionsPrefix, slepc.FNGetOptionsPrefix
        """
        cdef const char *prefix = NULL
        CHKERR( FNGetOptionsPrefix(self.fn, &prefix) )
        return bytes2str(prefix)

    def setFromOptions(self) -> None:
        """
        Set FN options from the options database.

        Collective.

        Notes
        -----
        To see all options, run your program with the ``-help``
        option.

        See Also
        --------
        setOptionsPrefix, slepc.FNSetFromOptions
        """
        CHKERR( FNSetFromOptions(self.fn) )

    def duplicate(self, comm: Comm | None = None) -> FN:
        """
        Duplicate the FN object copying all parameters.

        Collective.

        Duplicate the FN object copying all parameters, possibly with a
        different communicator.

        Parameters
        ----------
        comm
            MPI communicator; if not provided, it defaults to the
            object's communicator.

        Returns
        -------
        FN
            The new object.

        See Also
        --------
        create, slepc.FNDuplicate
        """
        cdef MPI_Comm ccomm = def_Comm(comm, PetscObjectComm(<PetscObject>self.fn))
        cdef FN fn = type(self)()
        CHKERR( FNDuplicate(self.fn, ccomm, &fn.fn) )
        return fn

    #

    def evaluateFunction(self, x: Scalar) -> Scalar:
        """
        Compute the value of the function :math:`f(x)` for a given x.

        Not collective.

        Parameters
        ----------
        x
            Value where the function must be evaluated.

        Returns
        -------
        Scalar
            The result of :math:`f(x)`.

        Notes
        -----
        Scaling factors are taken into account, so the actual function
        evaluation will return :math:`b f(a x)`.

        See Also
        --------
        evaluateDerivative, evaluateFunctionMat, setScale, slepc.FNEvaluateFunction
        """
        cdef PetscScalar sval = 0
        cdef PetscScalar sarg = asScalar(x)
        CHKERR( FNEvaluateFunction(self.fn, sarg, &sval) )
        return toScalar(sval)

    def evaluateDerivative(self, x: Scalar) -> Scalar:
        """
        Compute the value of the derivative :math:`f'(x)` for a given x.

        Not collective.

        Parameters
        ----------
        x
            Value where the derivative must be evaluated.

        Returns
        -------
        Scalar
            The result of :math:`f'(x)`.

        Notes
        -----
        Scaling factors are taken into account, so the actual derivative
        evaluation will return :math:`ab f'(a x)`.

        See Also
        --------
        evaluateFunction, setScale, slepc.FNEvaluateDerivative
        """
        cdef PetscScalar sval = 0
        cdef PetscScalar sarg = asScalar(x)
        CHKERR( FNEvaluateDerivative(self.fn, sarg, &sval) )
        return toScalar(sval)

    def evaluateFunctionMat(self, Mat A, Mat B: Mat | None = None) -> Mat:
        """
        Compute the value of the function :math:`f(A)` for a given matrix A.

        Logically collective.

        Parameters
        ----------
        A
            Matrix on which the function must be evaluated.
        B
            Placeholder for the result.

        Returns
        -------
        petsc4py.PETSc.Mat
            The result of :math:`f(A)`.

        Notes
        -----
        Scaling factors are taken into account, so the actual function
        evaluation will return :math:`b f(a A)`.

        See Also
        --------
        evaluateFunction, evaluateFunctionMatVec, slepc.FNEvaluateFunctionMat
        """
        if B is None: B = A.duplicate()
        CHKERR( FNEvaluateFunctionMat(self.fn, A.mat, B.mat) )
        return B

    def evaluateFunctionMatVec(self, Mat A, Vec v: Vec | None = None) -> Vec:
        """
        Compute the first column of the matrix :math:`f(A)`.

        Logically collective.

        Parameters
        ----------
        A
            Matrix on which the function must be evaluated.

        Returns
        -------
        petsc4py.PETSc.Vec
            The first column of the result :math:`f(A)`.

        Notes
        -----
        This operation is similar to `evaluateFunctionMat()` but returns only
        the first column of :math:`f(A)`, hence saving computations in most
        cases.

        See Also
        --------
        evaluateFunctionMat, slepc.FNEvaluateFunctionMatVec
        """
        if v is None: v = A.createVecs('left')
        CHKERR( FNEvaluateFunctionMatVec(self.fn, A.mat, v.vec) )
        return v

    def setScale(self, alpha: Scalar | None = None, beta: Scalar | None = None) -> None:
        """
        Set the scaling parameters that define the matematical function.

        Logically collective.

        Parameters
        ----------
        alpha
            Inner scaling (argument), default is 1.0.
        beta
            Outer scaling (result), default is 1.0.

        See Also
        --------
        getScale, evaluateFunction, slepc.FNSetScale
        """
        cdef PetscScalar aval = 1.0
        cdef PetscScalar bval = 1.0
        if alpha is not None: aval = asScalar(alpha)
        if beta  is not None: bval = asScalar(beta)
        CHKERR( FNSetScale(self.fn, aval, bval) )

    def getScale(self) -> tuple[Scalar, Scalar]:
        """
        Get the scaling parameters that define the matematical function.

        Not collective.

        Returns
        -------
        alpha: Scalar
            Inner scaling (argument).
        beta: Scalar
            Outer scaling (result).

        See Also
        --------
        setScale, slepc.FNGetScale
        """
        cdef PetscScalar aval = 0, bval = 0
        CHKERR( FNGetScale(self.fn, &aval, &bval) )
        return (toScalar(aval), toScalar(bval))

    def setMethod(self, meth: int) -> None:
        """
        Set the method to be used to evaluate functions of matrices.

        Logically collective.

        Parameters
        ----------
        meth
            An index identifying the method.

        Notes
        -----
        In some `FN` types there are more than one algorithms available
        for computing matrix functions. In that case, this function allows
        choosing the wanted method.

        If ``meth`` is currently set to 0 and the input argument of
        `FN.evaluateFunctionMat()` is a symmetric/Hermitian matrix, then
        the computation is done via the eigendecomposition, rather than
        with the general algorithm.

        See Also
        --------
        getMethod, slepc.FNSetMethod
        """
        cdef PetscInt val = asInt(meth)
        CHKERR( FNSetMethod(self.fn, val) )

    def getMethod(self) -> int:
        """
        Get the method currently used for matrix functions.

        Not collective.

        Returns
        -------
        int
            An index identifying the method.

        See Also
        --------
        setMethod, slepc.FNGetMethod
        """
        cdef PetscInt val = 0
        CHKERR( FNGetMethod(self.fn, &val) )
        return toInt(val)

    def setParallel(self, pmode: ParallelType) -> None:
        """
        Set the mode of operation in parallel runs.

        Logically collective.

        Parameters
        ----------
        pmode
            The parallel mode.

        Notes
        -----
        This is relevant only when the function is evaluated on a matrix, with
        either `evaluateFunctionMat()` or `evaluateFunctionMatVec()`.

        See Also
        --------
        evaluateFunctionMat, getParallel, slepc.FNSetParallel
        """
        cdef SlepcFNParallelType val = pmode
        CHKERR( FNSetParallel(self.fn, val) )

    def getParallel(self) -> ParallelType:
        """
        Get the mode of operation in parallel runs.

        Not collective.

        Returns
        -------
        ParallelType
            The parallel mode.

        See Also
        --------
        setParallel, slepc.FNGetParallel
        """
        cdef SlepcFNParallelType val = FN_PARALLEL_REDUNDANT
        CHKERR( FNGetParallel(self.fn, &val) )
        return val

    #

    def setRationalNumerator(self, alpha: Sequence[Scalar]) -> None:
        """
        Set the coefficients of the numerator of the rational function.

        Logically collective.

        Parameters
        ----------
        alpha
            Coefficients.

        See Also
        --------
        setRationalDenominator, slepc.FNRationalSetNumerator
        """
        cdef PetscInt na = 0
        cdef PetscScalar *a = NULL
        cdef object tmp1 = iarray_s(alpha, &na, &a)
        CHKERR( FNRationalSetNumerator(self.fn, na, a) )

    def getRationalNumerator(self) -> ArrayScalar:
        """
        Get the coefficients of the numerator of the rational function.

        Not collective.

        Returns
        -------
        ArrayScalar
            Coefficients.

        See Also
        --------
        setRationalNumerator, slepc.FNRationalGetNumerator
        """
        cdef PetscInt np = 0
        cdef PetscScalar *coeff = NULL
        CHKERR( FNRationalGetNumerator(self.fn, &np, &coeff) )
        cdef object ocoeff = None
        try:
            ocoeff = array_s(np, coeff)
        finally:
            CHKERR( PetscFree(coeff) )
        return ocoeff

    def setRationalDenominator(self, alpha: Sequence[Scalar]) -> None:
        """
        Set the coefficients of the denominator of the rational function.

        Logically collective.

        Parameters
        ----------
        alpha
            Coefficients.

        See Also
        --------
        setRationalNumerator, slepc.FNRationalSetDenominator
        """
        cdef PetscInt na = 0
        cdef PetscScalar *a = NULL
        cdef object tmp1 = iarray_s(alpha, &na, &a)
        CHKERR( FNRationalSetDenominator(self.fn, na, a) )

    def getRationalDenominator(self) -> ArrayScalar:
        """
        Get the coefficients of the denominator of the rational function.

        Not collective.

        Returns
        -------
        ArrayScalar
            Coefficients.

        See Also
        --------
        setRationalDenominator, slepc.FNRationalGetDenominator
        """
        cdef PetscInt np = 0
        cdef PetscScalar *coeff = NULL
        CHKERR( FNRationalGetDenominator(self.fn, &np, &coeff) )
        cdef object ocoeff = None
        try:
            ocoeff = array_s(np, coeff)
        finally:
            CHKERR( PetscFree(coeff) )
        return ocoeff

    def setCombineChildren(self, comb: CombineType, FN f1, FN f2) -> None:
        """
        Set the two child functions that constitute this combined function.

        Logically collective.

        Set the two child functions that constitute this combined function,
        and the way they must be combined.

        Parameters
        ----------
        comb
            How to combine the functions (addition, multiplication, division,
            composition).
        f1
            First function.
        f2
            Second function.

        See Also
        --------
        getCombineChildren, slepc.FNCombineSetChildren
        """
        cdef SlepcFNCombineType val = comb
        CHKERR( FNCombineSetChildren(self.fn, val, f1.fn, f2.fn) )

    def getCombineChildren(self) -> tuple[CombineType, FN, FN]:
        """
        Get the two child functions that constitute this combined function.

        Not collective.

        Get the two child functions that constitute this combined
        function, and the way they must be combined.

        Returns
        -------
        comb: CombineType
            How to combine the functions (addition, multiplication, division,
            composition).
        f1: FN
            First function.
        f2: FN
            Second function.

        See Also
        --------
        setCombineChildren, slepc.FNCombineGetChildren
        """
        cdef SlepcFNCombineType comb
        cdef FN f1 = FN()
        cdef FN f2 = FN()
        CHKERR( FNCombineGetChildren(self.fn, &comb, &f1.fn, &f2.fn) )
        CHKERR( PetscINCREF(f1.obj) )
        CHKERR( PetscINCREF(f2.obj) )
        return (comb, f1, f2)

    def setPhiIndex(self, k: int) -> None:
        """
        Set the index of the phi-function.

        Logically collective.

        Parameters
        ----------
        k
            The index.

        Notes
        -----
        If not set, the default index is 1.

        See Also
        --------
        getPhiIndex, slepc.FNPhiSetIndex
        """
        cdef PetscInt val = asInt(k)
        CHKERR( FNPhiSetIndex(self.fn, val) )

    def getPhiIndex(self) -> int:
        """
        Get the index of the phi-function.

        Not collective.

        Returns
        -------
        int
            The index.

        See Also
        --------
        setPhiIndex, slepc.FNPhiGetIndex
        """
        cdef PetscInt val = 0
        CHKERR( FNPhiGetIndex(self.fn, &val) )
        return toInt(val)

    #

    property method:
        """The method to be used to evaluate functions of matrices."""
        def __get__(self) -> int:
            return self.getMethod()
        def __set__(self, value):
            self.setMethod(value)

    property parallel:
        """The mode of operation in parallel runs."""
        def __get__(self) -> FNParallelType:
            return self.getParallel()
        def __set__(self, value):
            self.setParallel(value)

# -----------------------------------------------------------------------------

del FNType
del FNCombineType
del FNParallelType

# -----------------------------------------------------------------------------
