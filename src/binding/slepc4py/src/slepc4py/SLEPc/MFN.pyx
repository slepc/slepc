# -----------------------------------------------------------------------------

class MFNType(object):
    """
    MFN type.

    Action of a matrix function on a vector.

    - `KRYLOV`:  Restarted Krylov solver.
    - `EXPOKIT`: Implementation of the method in Expokit.
    """
    KRYLOV   = S_(MFNKRYLOV)
    EXPOKIT  = S_(MFNEXPOKIT)

class MFNConvergedReason(object):
    """
    MFN convergence reasons.

    - 'MFN_CONVERGED_TOL': All eigenpairs converged to requested tolerance.
    - 'MFN_CONVERGED_ITS': Solver completed the requested number of steps.
    - 'MFN_DIVERGED_ITS': Maximum number of iterations exceeded.
    - 'MFN_DIVERGED_BREAKDOWN': Generic breakdown in method.
    """
    CONVERGED_TOL       = MFN_CONVERGED_TOL
    CONVERGED_ITS       = MFN_CONVERGED_ITS
    DIVERGED_ITS        = MFN_DIVERGED_ITS
    DIVERGED_BREAKDOWN  = MFN_DIVERGED_BREAKDOWN
    CONVERGED_ITERATING = MFN_CONVERGED_ITERATING
    ITERATING           = MFN_CONVERGED_ITERATING

# -----------------------------------------------------------------------------

cdef class MFN(Object):

    """MFN."""

    Type            = MFNType
    ConvergedReason = MFNConvergedReason

    def __cinit__(self):
        self.obj = <PetscObject*> &self.mfn
        self.mfn = NULL

    def view(self, Viewer viewer=None) -> None:
        """
        Print the MFN data structure.

        Collective.

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
        Destroy the MFN object.

        Logically collective.
        """
        CHKERR( MFNDestroy(&self.mfn) )
        self.mfn = NULL
        return self

    def reset(self) -> None:
        """
        Reset the MFN object.

        Collective.
        """
        CHKERR( MFNReset(self.mfn) )

    def create(self, comm: Comm | None = None) -> Self:
        """
        Create the MFN object.

        Collective.

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
        Set the particular solver to be used in the MFN object.

        Logically collective.

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
        Get the MFN type of this object.

        Not collective.

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
        Get the prefix used for searching for all MFN options in the database.

        Not collective.

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
        Set the prefix used for searching for all MFN options in the database.

        Logically collective.

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
        Append to the prefix used for searching for all MFN options in the database.

        Logically collective.

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
        Set MFN options from the options database.

        Collective.

        Set MFN options from the options database. This routine must
        be called before `setUp()` if the user is to be allowed to set
        the solver type.
        """
        CHKERR( MFNSetFromOptions(self.mfn) )

    def getTolerances(self) -> tuple[float, int]:
        """
        Get the tolerance and maximum iteration count.

        Not collective.

        Get the tolerance and maximum iteration count used by the
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
        Set the tolerance and maximum iteration count.

        Logically collective.

        Set the tolerance and maximum iteration count used by the
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
        Get the dimension of the subspace used by the solver.

        Not collective.

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
        Set the dimension of the subspace to be used by the solver.

        Logically collective.

        Parameters
        ----------
        ncv
            Maximum dimension of the subspace to be used by the solver.
        """
        cdef PetscInt ival = asInt(ncv)
        CHKERR( MFNSetDimensions(self.mfn, ival) )

    def getFN(self) -> FN:
        """
        Get the math function object associated to the MFN object.

        Not collective.

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
        Set a math function object associated to the MFN object.

        Collective.

        Parameters
        ----------
        fn
            The math function context.
        """
        CHKERR( MFNSetFN(self.mfn, fn.fn) )

    def getBV(self) -> BV:
        """
        Get the basis vector object associated to the MFN object.

        Not collective.

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
        Set a basis vector object associated to the MFN object.

        Collective.

        Parameters
        ----------
        bv
            The basis vectors context.
        """
        CHKERR( MFNSetBV(self.mfn, bv.bv) )

    def getOperator(self) -> petsc4py.PETSc.Mat:
        """
        Get the matrix associated with the MFN object.

        Collective.

        Returns
        -------
        petsc4py.PETSc.Mat
            The matrix for which the matrix function is to be computed.
        """
        cdef Mat A = Mat()
        CHKERR( MFNGetOperator(self.mfn, &A.mat) )
        CHKERR( PetscINCREF(A.obj) )
        return A

    def setOperator(self, Mat A) -> None:
        """
        Set the matrix associated with the MFN object.

        Collective.

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
        Append a monitor function to the list of monitors.

        Logically collective.
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
        """Get the list of monitor functions."""
        return self.get_attr('__monitor__')

    def cancelMonitor(self) -> None:
        """
        Clear all monitors for an `MFN` object.

        Logically collective.
        """
        CHKERR( MFNMonitorCancel(self.mfn) )
        self.set_attr('__monitor__', None)

    #

    def setUp(self) -> None:
        """
        Set up all the necessary internal data structures.

        Collective.

        Set up all the internal data structures necessary for the execution
        of the eigensolver.
        """
        CHKERR( MFNSetUp(self.mfn) )

    def solve(self, Vec b, Vec x) -> None:
        """
        Solve the matrix function problem.

        Collective.

        Given a vector :math:`b`, the vector :math:`x = f(A) b` is
        returned.

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
        Solve the transpose matrix function problem.

        Collective.

        Given a vector :math:`b`, the vector :math:`x = f(A^T) b` is
        returned.

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
        Get the current iteration number.

        Not collective.

        Get the current iteration number. If the call to `solve()` is
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
        Get the reason why the `solve()` iteration was stopped.

        Not collective.

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
        Set `solve()` to generate an error if the solver does not converge.

        Logically collective.

        Parameters
        ----------
        flg
            True indicates you want the error generated.
        """
        cdef PetscBool tval = flg
        CHKERR( MFNSetErrorIfNotConverged(self.mfn, tval) )

    def getErrorIfNotConverged(self) -> bool:
        """
        Get if `solve()` generates an error if the solver does not converge.

        Not collective.

        Get a flag indicating whether `solve()` will generate an error if the
        solver does not converge.

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
        """The tolerance count used by the MFN convergence tests."""
        def __get__(self) -> float:
            return self.getTolerances()[0]
        def __set__(self, value):
            self.setTolerances(tol=value)

    property max_it:
        """The maximum iteration count used by the MFN convergence tests."""
        def __get__(self) -> int:
            return self.getTolerances()[1]
        def __set__(self, value):
            self.setTolerances(max_it=value)

    property fn:
        """The math function (FN) object associated to the MFN object."""
        def __get__(self) -> FN:
            return self.getFN()
        def __set__(self, value):
            self.setBV(value)

    property bv:
        """The basis vectors (BV) object associated to the MFN object."""
        def __get__(self) -> BV:
            return self.getFN()
        def __set__(self, value):
            self.setBV(value)

# -----------------------------------------------------------------------------

del MFNType
del MFNConvergedReason

# -----------------------------------------------------------------------------
