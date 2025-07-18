# -----------------------------------------------------------------------------

class MFNType(object):
    """
    MFN type

    Action of a matrix function on a vector.

    - `KRYLOV`:  Restarted Krylov solver.
    - `EXPOKIT`: Implementation of the method in Expokit.
    """
    KRYLOV   = S_(MFNKRYLOV)
    EXPOKIT  = S_(MFNEXPOKIT)

class MFNConvergedReason(object):
    CONVERGED_TOL       = MFN_CONVERGED_TOL
    CONVERGED_ITS       = MFN_CONVERGED_ITS
    DIVERGED_ITS        = MFN_DIVERGED_ITS
    DIVERGED_BREAKDOWN  = MFN_DIVERGED_BREAKDOWN
    CONVERGED_ITERATING = MFN_CONVERGED_ITERATING
    ITERATING           = MFN_CONVERGED_ITERATING

# -----------------------------------------------------------------------------

cdef class MFN(Object):

    """
    MFN
    """

    Type            = MFNType
    ConvergedReason = MFNConvergedReason

    def __cinit__(self):
        self.obj = <PetscObject*> &self.mfn
        self.mfn = NULL

    def view(self, Viewer viewer=None) -> None:
        """
        Prints the MFN data structure.

        Parameters
        ----------
        viewer
            Visualization context; if not provided, the standard
            output is used.
        """
        cdef PetscViewer vwr = def_Viewer(viewer)
        CHKERR( MFNView(self.mfn, vwr) )

    def destroy(self) -> Self:
        """
        Destroys the MFN object.
        """
        CHKERR( MFNDestroy(&self.mfn) )
        self.mfn = NULL
        return self

    def reset(self) -> None:
        """
        Resets the MFN object.
        """
        CHKERR( MFNReset(self.mfn) )

    def create(self, comm: Comm | None = None) -> Self:
        """
        Creates the MFN object.

        Parameters
        ----------
        comm
            MPI communicator. If not provided, it defaults to all processes.
        """
        cdef MPI_Comm ccomm = def_Comm(comm, SLEPC_COMM_DEFAULT())
        cdef SlepcMFN newmfn = NULL
        CHKERR( MFNCreate(ccomm, &newmfn) )
        CHKERR( SlepcCLEAR(self.obj) ); self.mfn = newmfn
        return self

    def setType(self, mfn_type: Type | str) -> None:
        """
        Selects the particular solver to be used in the MFN object.

        Parameters
        ----------
        mfn_type
            The solver to be used.
        """
        cdef SlepcMFNType cval = NULL
        mfn_type = str2bytes(mfn_type, &cval)
        CHKERR( MFNSetType(self.mfn, cval) )

    def getType(self) -> str:
        """
        Gets the MFN type of this object.

        Returns
        -------
        str
            The solver currently being used.
        """
        cdef SlepcMFNType mfn_type = NULL
        CHKERR( MFNGetType(self.mfn, &mfn_type) )
        return bytes2str(mfn_type)

    def getOptionsPrefix(self) -> str:
        """
        Gets the prefix used for searching for all MFN options in the
        database.

        Returns
        -------
        str
            The prefix string set for this MFN object.
        """
        cdef const char *prefix = NULL
        CHKERR( MFNGetOptionsPrefix(self.mfn, &prefix) )
        return bytes2str(prefix)

    def setOptionsPrefix(self, prefix: str | None = None) -> None:
        """
        Sets the prefix used for searching for all MFN options in the
        database.

        Parameters
        ----------
        prefix
            The prefix string to prepend to all MFN option requests.
        """
        cdef const char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( MFNSetOptionsPrefix(self.mfn, cval) )

    def appendOptionsPrefix(self, prefix: str | None = None) -> None:
        """
        Appends to the prefix used for searching for all MFN options
        in the database.

        Parameters
        ----------
        prefix
            The prefix string to prepend to all MFN option requests.
        """
        cdef const char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( MFNAppendOptionsPrefix(self.mfn, cval) )

    def setFromOptions(self) -> None:
        """
        Sets MFN options from the options database. This routine must
        be called before `setUp()` if the user is to be allowed to set
        the solver type.
        """
        CHKERR( MFNSetFromOptions(self.mfn) )

    def getTolerances(self) -> tuple[float, int]:
        """
        Gets the tolerance and maximum iteration count used by the
        default MFN convergence tests.

        Returns
        -------
        tol: float
            The convergence tolerance.
        max_it: int
            The maximum number of iterations
        """
        cdef PetscReal rval = 0
        cdef PetscInt  ival = 0
        CHKERR( MFNGetTolerances(self.mfn, &rval, &ival) )
        return (toReal(rval), toInt(ival))

    def setTolerances(self, tol: float | None = None, max_it: int | None = None) -> None:
        """
        Sets the tolerance and maximum iteration count used by the
        default MFN convergence tests.

        Parameters
        ----------
        tol
            The convergence tolerance.
        max_it
            The maximum number of iterations
        """
        cdef PetscReal rval = PETSC_CURRENT
        cdef PetscInt  ival = PETSC_CURRENT
        if tol    is not None: rval = asReal(tol)
        if max_it is not None: ival = asInt(max_it)
        CHKERR( MFNSetTolerances(self.mfn, rval, ival) )

    def getDimensions(self) -> int:
        """
        Gets the dimension of the subspace used by the solver.

        Returns
        -------
        int
            Maximum dimension of the subspace to be used by the solver.
        """
        cdef PetscInt ival = 0
        CHKERR( MFNGetDimensions(self.mfn, &ival) )
        return toInt(ival)

    def setDimensions(self, ncv: int) -> None:
        """
        Sets the dimension of the subspace to be used by the solver.

        Parameters
        ----------
        ncv
            Maximum dimension of the subspace to be used by the solver.
        """
        cdef PetscInt ival = asInt(ncv)
        CHKERR( MFNSetDimensions(self.mfn, ival) )

    def getFN(self) -> FN:
        """
        Obtain the math function object associated to the MFN object.

        Returns
        -------
        FN
            The math function context.
        """
        cdef FN fn = FN()
        CHKERR( MFNGetFN(self.mfn, &fn.fn) )
        CHKERR( PetscINCREF(fn.obj) )
        return fn

    def setFN(self, FN fn) -> None:
        """
        Associates a math function object to the MFN object.

        Parameters
        ----------
        fn
            The math function context.
        """
        CHKERR( MFNSetFN(self.mfn, fn.fn) )

    def getBV(self) -> BV:
        """
        Obtain the basis vector object associated to the MFN object.

        Returns
        -------
        BV
            The basis vectors context.
        """
        cdef BV bv = BV()
        CHKERR( MFNGetBV(self.mfn, &bv.bv) )
        CHKERR( PetscINCREF(bv.obj) )
        return bv

    def setBV(self, BV bv) -> None:
        """
        Associates a basis vector object to the MFN object.

        Parameters
        ----------
        bv
            The basis vectors context.
        """
        CHKERR( MFNSetBV(self.mfn, bv.bv) )

    def getOperator(self) -> Mat:
        """
        Gets the matrix associated with the MFN object.

        Returns
        -------
        Mat
            The matrix for which the matrix function is to be computed.
        """
        cdef Mat A = Mat()
        CHKERR( MFNGetOperator(self.mfn, &A.mat) )
        CHKERR( PetscINCREF(A.obj) )
        return A

    def setOperator(self, Mat A) -> None:
        """
        Sets the matrix associated with the MFN object.

        Parameters
        ----------
        A
            The problem matrix.
        """
        CHKERR( MFNSetOperator(self.mfn, A.mat) )

    #

    def setMonitor(
        self,
        monitor: MFNMonitorFunction | None,
        args: tuple[Any, ...] | None = None,
        kargs: dict[str, Any] | None = None,
    ) -> None:
        """
        Appends a monitor function to the list of monitors.
        """
        if monitor is None: return
        cdef object monitorlist = self.get_attr('__monitor__')
        if monitorlist is None:
            monitorlist = []
            self.set_attr('__monitor__', monitorlist)
            CHKERR( MFNMonitorSet(self.mfn, MFN_Monitor, NULL, NULL) )
        if args is None: args = ()
        if kargs is None: kargs = {}
        monitorlist.append((monitor, args, kargs))

    def getMonitor(self) -> MFNMonitorFunction:
        """
        Gets the list of monitor functions.
        """
        return self.get_attr('__monitor__')

    def cancelMonitor(self) -> None:
        """
        Clears all monitors for an `MFN` object.
        """
        CHKERR( MFNMonitorCancel(self.mfn) )
        self.set_attr('__monitor__', None)

    #

    def setUp(self) -> None:
        """
        Sets up all the internal data structures necessary for the
        execution of the eigensolver.
        """
        CHKERR( MFNSetUp(self.mfn) )

    def solve(self, Vec b, Vec x) -> None:
        """
        Solves the matrix function problem. Given a vector b, the
        vector x = f(A)*b is returned.

        Parameters
        ----------
        b
            The right hand side vector.
        x
            The solution.
        """
        CHKERR( MFNSolve(self.mfn, b.vec, x.vec) )

    def solveTranspose(self, Vec b, Vec x) -> None:
        """
        Solves the transpose matrix function problem. Given a vector b, the
        vector x = f(A^T)*b is returned.

        Parameters
        ----------
        b
            The right hand side vector.
        x
            The solution.
        """
        CHKERR( MFNSolveTranspose(self.mfn, b.vec, x.vec) )

    def getIterationNumber(self) -> int:
        """
        Gets the current iteration number. If the call to `solve()` is
        complete, then it returns the number of iterations carried out
        by the solution method.

        Returns
        -------
        int
            Iteration number.
        """
        cdef PetscInt ival = 0
        CHKERR( MFNGetIterationNumber(self.mfn, &ival) )
        return toInt(ival)

    def getConvergedReason(self) -> ConvergedReason:
        """
        Gets the reason why the `solve()` iteration was stopped.

        Returns
        -------
        ConvergedReason
            Negative value indicates diverged, positive value converged.
        """
        cdef SlepcMFNConvergedReason val = MFN_CONVERGED_ITERATING
        CHKERR( MFNGetConvergedReason(self.mfn, &val) )
        return val

    def setErrorIfNotConverged(self, flg: bool = True) -> None:
        """
        Causes `solve()` to generate an error if the solver has not converged.

        Parameters
        ----------
        flg
            True indicates you want the error generated.
        """
        cdef PetscBool tval = flg
        CHKERR( MFNSetErrorIfNotConverged(self.mfn, tval) )

    def getErrorIfNotConverged(self) -> bool:
        """
        Return a flag indicating whether `solve()` will generate an
        error if the solver does not converge.

        Returns
        -------
        bool
            True indicates you want the error generated.
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( MFNGetErrorIfNotConverged(self.mfn, &tval) )
        return toBool(tval)

    #

    property tol:
        def __get__(self):
            return self.getTolerances()[0]
        def __set__(self, value):
            self.setTolerances(tol=value)

    property max_it:
        def __get__(self):
            return self.getTolerances()[1]
        def __set__(self, value):
            self.setTolerances(max_it=value)

    property fn:
        def __get__(self):
            return self.getFN()
        def __set__(self, value):
            self.setBV(value)

    property bv:
        def __get__(self):
            return self.getFN()
        def __set__(self, value):
            self.setBV(value)

# -----------------------------------------------------------------------------

del MFNType
del MFNConvergedReason

# -----------------------------------------------------------------------------
