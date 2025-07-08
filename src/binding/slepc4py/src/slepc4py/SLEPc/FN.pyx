# -----------------------------------------------------------------------------

class FNType(object):
    """
    FN type
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
    FN type of combination of child functions

    - `ADD`:       Addition         f(x) = f1(x)+f2(x)
    - `MULTIPLY`:  Multiplication   f(x) = f1(x)*f2(x)
    - `DIVIDE`:    Division         f(x) = f1(x)/f2(x)
    - `COMPOSE`:   Composition      f(x) = f2(f1(x))
    """
    ADD      = FN_COMBINE_ADD
    MULTIPLY = FN_COMBINE_MULTIPLY
    DIVIDE   = FN_COMBINE_DIVIDE
    COMPOSE  = FN_COMBINE_COMPOSE

class FNParallelType(object):
    """
    FN parallel types

    - `REDUNDANT`:    Every process performs the computation redundantly.
    - `SYNCHRONIZED`: The first process sends the result to the rest.
    """
    REDUNDANT    = FN_PARALLEL_REDUNDANT
    SYNCHRONIZED = FN_PARALLEL_SYNCHRONIZED

# -----------------------------------------------------------------------------

cdef class FN(Object):

    """
    FN
    """

    Type         = FNType
    CombineType  = FNCombineType
    ParallelType = FNParallelType

    def __cinit__(self):
        self.obj = <PetscObject*> &self.fn
        self.fn = NULL

    def view(self, Viewer viewer=None):
        """
        Prints the FN data structure.

        Parameters
        ----------
        viewer: Viewer, optional
                Visualization context; if not provided, the standard
                output is used.
        """
        cdef PetscViewer vwr = def_Viewer(viewer)
        CHKERR( FNView(self.fn, vwr) )

    def destroy(self):
        """
        Destroys the FN object.
        """
        CHKERR( FNDestroy(&self.fn) )
        self.fn = NULL
        return self

    def create(self, comm=None):
        """
        Creates the FN object.

        Parameters
        ----------
        comm: Comm, optional
              MPI communicator; if not provided, it defaults to all
              processes.
        """
        cdef MPI_Comm ccomm = def_Comm(comm, SLEPC_COMM_DEFAULT())
        cdef SlepcFN newfn = NULL
        CHKERR( FNCreate(ccomm, &newfn) )
        CHKERR( SlepcCLEAR(self.obj) ); self.fn = newfn
        return self

    def setType(self, fn_type):
        """
        Selects the type for the FN object.

        Parameters
        ----------
        fn_type: `FN.Type` enumerate
                  The inner product type to be used.
        """
        cdef SlepcFNType cval = NULL
        fn_type = str2bytes(fn_type, &cval)
        CHKERR( FNSetType(self.fn, cval) )

    def getType(self):
        """
        Gets the FN type of this object.

        Returns
        -------
        type: `FN.Type` enumerate
              The inner product type currently being used.
        """
        cdef SlepcFNType fn_type = NULL
        CHKERR( FNGetType(self.fn, &fn_type) )
        return bytes2str(fn_type)

    def setOptionsPrefix(self, prefix):
        """
        Sets the prefix used for searching for all FN options in the
        database.

        Parameters
        ----------
        prefix: string
                The prefix string to prepend to all FN option
                requests.

        Notes
        -----
        A hyphen (``-``) must NOT be given at the beginning of the
        prefix name.  The first character of all runtime options is
        AUTOMATICALLY the hyphen.
        """
        cdef const char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( FNSetOptionsPrefix(self.fn, cval) )

    def getOptionsPrefix(self):
        """
        Gets the prefix used for searching for all FN options in the
        database.

        Returns
        -------
        prefix: string
                The prefix string set for this FN object.
        """
        cdef const char *prefix = NULL
        CHKERR( FNGetOptionsPrefix(self.fn, &prefix) )
        return bytes2str(prefix)

    def setFromOptions(self):
        """
        Sets FN options from the options database.

        Notes
        -----
        To see all options, run your program with the ``-help``
        option.
        """
        CHKERR( FNSetFromOptions(self.fn) )

    def duplicate(self, comm=None):
        """
        Duplicate the FN object copying all parameters, possibly with a
        different communicator.

        Parameters
        ----------
        comm: Comm, optional
              MPI communicator; if not provided, it defaults to the
              object's communicator.
        """
        cdef MPI_Comm ccomm = def_Comm(comm, PetscObjectComm(<PetscObject>self.fn))
        cdef FN fn = type(self)()
        CHKERR( FNDuplicate(self.fn, ccomm, &fn.fn) )
        return fn

    #

    def evaluateFunction(self, x):
        """
        Computes the value of the function f(x) for a given x.

        Parameters
        ----------
        x: scalar
            Value where the function must be evaluated.

        Returns
        -------
        y: scalar
            The result of f(x).
        """
        cdef PetscScalar sval = 0
        cdef PetscScalar sarg = asScalar(x)
        CHKERR( FNEvaluateFunction(self.fn, sarg, &sval) )
        return toScalar(sval)

    def evaluateDerivative(self, x):
        """
        Computes the value of the derivative f'(x) for a given x.

        Parameters
        ----------
        x: scalar
            Value where the derivative must be evaluated.

        Returns
        -------
        y: scalar
            The result of f'(x).
        """
        cdef PetscScalar sval = 0
        cdef PetscScalar sarg = asScalar(x)
        CHKERR( FNEvaluateDerivative(self.fn, sarg, &sval) )
        return toScalar(sval)

    def evaluateFunctionMat(self, Mat A, Mat B=None):
        """
        Computes the value of the function f(A) for a given matrix A.

        Parameters
        ----------
        A: Mat
           Matrix on which the function must be evaluated.
        B: Mat, optional
           Placeholder for the result.

        Returns
        -------
        B: Mat
           The result of f(A).
        """
        if B is None: B = A.duplicate()
        CHKERR( FNEvaluateFunctionMat(self.fn, A.mat, B.mat) )
        return B

    def evaluateFunctionMatVec(self, Mat A, Vec v=None):
        """
        Computes the first column of the matrix f(A) for a given matrix A.

        Parameters
        ----------
        A: Mat
           Matrix on which the function must be evaluated.

        Returns
        -------
        v: Vec
           The first column of the result f(A).
        """
        if v is None: v = A.createVecs('left')
        CHKERR( FNEvaluateFunctionMatVec(self.fn, A.mat, v.vec) )
        return v

    def setScale(self, alpha=None, beta=None):
        """
        Sets the scaling parameters that define the matematical function.

        Parameters
        ----------
        alpha: scalar (possibly complex), optional
               Inner scaling (argument), default is 1.0.
        beta: scalar (possibly complex), optional
               Outer scaling (result), default is 1.0.
        """
        cdef PetscScalar aval = 1.0
        cdef PetscScalar bval = 1.0
        if alpha is not None: aval = asScalar(alpha)
        if beta  is not None: bval = asScalar(beta)
        CHKERR( FNSetScale(self.fn, aval, bval) )

    def getScale(self):
        """
        Gets the scaling parameters that define the matematical function.

        Returns
        -------
        alpha: scalar (possibly complex)
               Inner scaling (argument).
        beta: scalar (possibly complex)
               Outer scaling (result).
        """
        cdef PetscScalar aval = 0, bval = 0
        CHKERR( FNGetScale(self.fn, &aval, &bval) )
        return (toScalar(aval), toScalar(bval))

    def setMethod(self, meth):
        """
        Selects the method to be used to evaluate functions of matrices.

        Parameters
        ----------
        meth: int
              An index identifying the method.

        Notes
        -----
        In some `FN` types there are more than one algorithms available
        for computing matrix functions. In that case, this function allows
        choosing the wanted method.

        If `meth` is currently set to 0 and the input argument of
        `FN.evaluateFunctionMat()` is a symmetric/Hermitian matrix, then
        the computation is done via the eigendecomposition, rather than
        with the general algorithm.
        """
        cdef PetscInt val = asInt(meth)
        CHKERR( FNSetMethod(self.fn, val) )

    def getMethod(self):
        """
        Gets the method currently used for matrix functions.

        Returns
        -------
        meth: int
              An index identifying the method.
        """
        cdef PetscInt val = 0
        CHKERR( FNGetMethod(self.fn, &val) )
        return toInt(val)

    def setParallel(self, pmode):
        """
        Selects the mode of operation in parallel runs.

        Parameters
        ----------
        pmode: `FN.ParallelType` enumerate
               The parallel mode.
        """
        cdef SlepcFNParallelType val = pmode
        CHKERR( FNSetParallel(self.fn, val) )

    def getParallel(self):
        """
        Gets the mode of operation in parallel runs.

        Returns
        -------
        pmode: `FN.ParallelType` enumerate
               The parallel mode.
        """
        cdef SlepcFNParallelType val = FN_PARALLEL_REDUNDANT
        CHKERR( FNGetParallel(self.fn, &val) )
        return val

    #

    def setRationalNumerator(self, alpha):
        """
        Sets the coefficients of the numerator of the rational function.

        Parameters
        ----------
        alpha: array of scalars
            Coefficients.
        """
        cdef PetscInt na = 0
        cdef PetscScalar *a = NULL
        cdef object tmp1 = iarray_s(alpha, &na, &a)
        CHKERR( FNRationalSetNumerator(self.fn, na, a) )

    def getRationalNumerator(self):
        """
        Gets the coefficients of the numerator of the rational function.

        Returns
        -------
        alpha: array of scalars
            Coefficients.
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

    def setRationalDenominator(self, alpha):
        """
        Sets the coefficients of the denominator of the rational function.

        Parameters
        ----------
        alpha: array of scalars
            Coefficients.
        """
        cdef PetscInt na = 0
        cdef PetscScalar *a = NULL
        cdef object tmp1 = iarray_s(alpha, &na, &a)
        CHKERR( FNRationalSetDenominator(self.fn, na, a) )

    def getRationalDenominator(self):
        """
        Gets the coefficients of the denominator of the rational function.

        Returns
        -------
        alpha: array of scalars
            Coefficients.
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

    def setCombineChildren(self, comb, FN f1, FN f2):
        """
        Sets the two child functions that constitute this combined
        function, and the way they must be combined.

        Parameters
        ----------
        comb: `FN.CombineType` enumerate
            How to combine the functions (addition, multiplication, division, composition).
        f1: FN
            First function.
        f2: FN
            Second function.
        """
        cdef SlepcFNCombineType val = comb
        CHKERR( FNCombineSetChildren(self.fn, val, f1.fn, f2.fn) )

    def getCombineChildren(self):
        """
        Gets the two child functions that constitute this combined
        function, and the way they must be combined.

        Returns
        -------
        comb: `FN.CombineType` enumerate
            How to combine the functions (addition, multiplication, division, composition).
        f1: FN
            First function.
        f2: FN
            Second function.
        """
        cdef SlepcFNCombineType comb
        cdef FN f1 = FN()
        cdef FN f2 = FN()
        CHKERR( FNCombineGetChildren(self.fn, &comb, &f1.fn, &f2.fn) )
        CHKERR( PetscINCREF(f1.obj) )
        CHKERR( PetscINCREF(f2.obj) )
        return (comb, f1, f2)

    def setPhiIndex(self, k):
        """
        Sets the index of the phi-function.

        Parameters
        ----------
        k: int
           The index.
        """
        cdef PetscInt val = asInt(k)
        CHKERR( FNPhiSetIndex(self.fn, val) )

    def getPhiIndex(self):
        """
        Gets the index of the phi-function.

        Returns
        -------
        k: int
           The index.
        """
        cdef PetscInt val = 0
        CHKERR( FNPhiGetIndex(self.fn, &val) )
        return toInt(val)

    #

    property method:
        def __get__(self):
            return self.getMethod()
        def __set__(self, value):
            self.setMethod(value)

    property parallel:
        def __get__(self):
            return self.getParallel()
        def __set__(self, value):
            self.setParallel(value)

# -----------------------------------------------------------------------------

del FNType
del FNCombineType
del FNParallelType

# -----------------------------------------------------------------------------
