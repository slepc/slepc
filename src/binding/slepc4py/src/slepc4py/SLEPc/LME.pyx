# -----------------------------------------------------------------------------

class LMEType(object):
    """
    LME type.

    - `KRYLOV`:  Restarted Krylov solver.
    """
    KRYLOV   = S_(LMEKRYLOV)

class LMEConvergedReason(object):
    """
    LME convergence reasons.

    - `CONVERGED_TOL`:       All eigenpairs converged to requested tolerance.
    - `DIVERGED_ITS`:        Maximum number of iterations exceeded.
    - `DIVERGED_BREAKDOWN`:  Solver failed due to breakdown.
    - `CONVERGED_ITERATING`: Iteration not finished yet.
    """
    CONVERGED_TOL       = LME_CONVERGED_TOL
    DIVERGED_ITS        = LME_DIVERGED_ITS
    DIVERGED_BREAKDOWN  = LME_DIVERGED_BREAKDOWN
    CONVERGED_ITERATING = LME_CONVERGED_ITERATING
    ITERATING           = LME_CONVERGED_ITERATING

class LMEProblemType(object):
    """
    LME problem type.

    - `LYAPUNOV`:      Continuous-time Lyapunov.
    - `SYLVESTER`:     Continuous-time Sylvester.
    - `GEN_LYAPUNOV`:  Generalized Lyapunov.
    - `GEN_SYLVESTER`: Generalized Sylvester.
    - `DT_LYAPUNOV`:   Discrete-time Lyapunov.
    - `STEIN`:         Stein.
    """
    LYAPUNOV        =  LME_LYAPUNOV
    SYLVESTER       =  LME_SYLVESTER
    GEN_LYAPUNOV    =  LME_GEN_LYAPUNOV
    GEN_SYLVESTER   =  LME_GEN_SYLVESTER
    DT_LYAPUNOV     =  LME_DT_LYAPUNOV
    STEIN           =  LME_STEIN

# -----------------------------------------------------------------------------

cdef class LME(Object):

    """LME."""

    Type            = LMEType
    ProblemType     = LMEProblemType
    ConvergedReason = LMEConvergedReason

    def __cinit__(self):
        self.obj = <PetscObject*> &self.lme
        self.lme = NULL

    def view(self, Viewer viewer=None) -> None:
        """
        Print the LME data structure.

        Collective.

        Parameters
        ----------
        viewer
            Visualization context; if not provided, the standard
            output is used.
        """
        cdef PetscViewer vwr = def_Viewer(viewer)
        CHKERR( LMEView(self.lme, vwr) )

    def destroy(self) -> Self:
        """
        Destroy the LME object.

        Collective.
        """
        CHKERR( LMEDestroy(&self.lme) )
        self.lme = NULL
        return self

    def reset(self) -> None:
        """
        Reset the LME object.

        Collective.
        """
        CHKERR( LMEReset(self.lme) )

    def create(self, comm: Comm | None = None) -> Self:
        """
        Create the LME object.

        Collective.

        Parameters
        ----------
        comm
            MPI communicator. If not provided, it defaults to all processes.
        """
        cdef MPI_Comm ccomm = def_Comm(comm, SLEPC_COMM_DEFAULT())
        cdef SlepcLME newlme = NULL
        CHKERR( LMECreate(ccomm, &newlme) )
        CHKERR( SlepcCLEAR(self.obj) ); self.lme = newlme
        return self

    def setType(self, lme_type: Type | str) -> None:
        """
        Set the particular solver to be used in the LME object.

        Logically collective.

        Parameters
        ----------
        lme_type
            The solver to be used.
        """
        cdef SlepcLMEType cval = NULL
        lme_type = str2bytes(lme_type, &cval)
        CHKERR( LMESetType(self.lme, cval) )

    def getType(self) -> str:
        """
        Get the LME type of this object.

        Not collective.

        Returns
        -------
        str
            The solver currently being used.
        """
        cdef SlepcLMEType lme_type = NULL
        CHKERR( LMEGetType(self.lme, &lme_type) )
        return bytes2str(lme_type)

    def setProblemType(self, lme_problem_type: ProblemType | str) -> None:
        """
        Set the LME problem type of this object.

        Logically collective.

        Parameters
        ----------
        lme_problem_type
            The problem type to be used.
        """
        cdef SlepcLMEProblemType val = lme_problem_type
        CHKERR( LMESetProblemType(self.lme, val) )

    def getProblemType(self) -> ProblemType:
        """
        Get the LME problem type of this object.

        Not collective.

        Returns
        -------
        ProblemType
            The problem type currently being used.
        """
        cdef SlepcLMEProblemType val = LME_LYAPUNOV
        CHKERR( LMEGetProblemType(self.lme, &val))
        return val

    def setCoefficients(self, Mat A, Mat B = None, Mat D = None, Mat E = None) -> None:
        """
        Set the coefficient matrices.

        Collective.

        Set the coefficient matrices that define the linear matrix equation
        to be solved.

        Parameters
        ----------
        A
            First coefficient matrix
        B
            Second coefficient matrix, optional
        D
            Third coefficient matrix, optional
        E
            Fourth coefficient matrix, optional
        """
        cdef PetscMat Amat = A.mat
        cdef PetscMat Bmat = B.mat if B is not None else <PetscMat>NULL
        cdef PetscMat Dmat = D.mat if D is not None else <PetscMat>NULL
        cdef PetscMat Emat = E.mat if E is not None else <PetscMat>NULL
        CHKERR( LMESetCoefficients(self.lme, Amat, Bmat, Dmat, Emat))

    def getCoefficients(self) -> tuple[Mat, Mat | None, Mat | None, Mat | None]:
        """
        Get the coefficient matrices of the matrix equation.

        Collective.

        Returns
        -------
        ``A``
            First coefficient matrix
        ``B``
            Second coefficient matrix, if available
        ``D``
            Third coefficient matrix, if available
        ``E``
            Fourth coefficient matrix, if available
        """
        cdef PetscMat Amat, Bmat, Dmat, Emat
        cdef Mat A = Mat(), B = None, D = None, E = None
        CHKERR( LMEGetCoefficients(self.lme, &Amat, &Bmat, &Dmat, &Emat) )
        A.mat = Amat
        CHKERR( PetscINCREF(A.obj) )
        if Bmat:
            B = Mat(); B.mat = Bmat; CHKERR( PetscINCREF(B.obj) )
        if Dmat:
            D = Mat(); D.mat = Dmat; CHKERR( PetscINCREF(D.obj) )
        if Emat:
            E = Mat(); E.mat = Emat; CHKERR( PetscINCREF(E.obj) )
        return (A, B, D, E)

    def setRHS(self, Mat C) -> None:
        """
        Set the right-hand side of the matrix equation.

        Collective.

        Set the right-hand side of the matrix equation, as a low-rank
        matrix.

        Parameters
        ----------
        C
            The right-hand side matrix
        """
        CHKERR( LMESetRHS(self.lme, C.mat) )

    def getRHS(self) -> Mat:
        """
        Get the right-hand side of the matrix equation.

        Collective.

        Returns
        -------
        C
            The low-rank matrix
        """
        cdef Mat C = Mat()
        CHKERR( LMEGetRHS(self.lme, &C.mat) )
        CHKERR( PetscINCREF(C.obj) )
        return C

    def setSolution(self, Mat X) -> None:
        """
        Set the placeholder for the solution of the matrix equation.

        Collective.

        Set the placeholder for the solution of the matrix
        equation, as a low-rank matrix.

        Parameters
        ----------
        X
            The solution matrix
        """
        CHKERR( LMESetSolution(self.lme, X.mat) )

    def getSolution(self) -> Mat:
        """
        Get the solution of the matrix equation.

        Collective.

        Returns
        -------
        X
            The low-rank matrix
        """
        cdef Mat X = Mat()
        CHKERR( LMEGetSolution(self.lme, &X.mat) )
        CHKERR( PetscINCREF(X.obj) )
        return X

    def getErrorEstimate(self) -> float:
        """
        Get the error estimate obtained during solve.

        Not collective.

        Returns
        -------
        float
            The error estimate
        """
        cdef PetscReal rval = 0
        CHKERR( LMEGetErrorEstimate(self.lme, &rval) )
        return toReal(rval)

    def computeError(self) -> float:
        """
        Compute the error associated with the last equation solved.

        Collective.

        Computes the error (based on the residual norm) associated with the
        last equation solved.

        Returns
        -------
        float
            The error
        """
        cdef PetscReal rval = 0
        CHKERR( LMEComputeError(self.lme, &rval) )
        return toReal(rval)

    def getOptionsPrefix(self) -> str:
        """
        Get the prefix used for searching for all LME options in the database.

        Not collective.

        Returns
        -------
        str
            The prefix string set for this LME object.
        """
        cdef const char *prefix = NULL
        CHKERR( LMEGetOptionsPrefix(self.lme, &prefix) )
        return bytes2str(prefix)

    def setOptionsPrefix(self, prefix: str | None = None) -> None:
        """
        Set the prefix used for searching for all LME options in the database.

        Logically collective.

        Parameters
        ----------
        prefix
            The prefix string to prepend to all LME option requests.
        """
        cdef const char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( LMESetOptionsPrefix(self.lme, cval) )

    def appendOptionsPrefix(self, prefix: str | None = None) -> None:
        """
        Append to the prefix used for searching in the database.

        Logically collective.

        Append to the prefix used for searching for all LME options in the
        database.

        Parameters
        ----------
        prefix
            The prefix string to prepend to all LME option requests.
        """
        cdef const char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( LMEAppendOptionsPrefix(self.lme, cval) )

    def setFromOptions(self) -> None:
        """
        Set LME options from the options database.

        Collective.

        Sets LME options from the options database. This routine must
        be called before `setUp()` if the user is to be allowed to set
        the solver type.
        """
        CHKERR( LMESetFromOptions(self.lme) )

    def getTolerances(self) -> tuple[float, int]:
        """
        Get the tolerance and maximum iteration count.

        Not collective.

        Get the tolerance and maximum iteration count used by the
        default LME convergence tests.

        Returns
        -------
        tol: float
            The convergence tolerance.
        max_it: int
            The maximum number of iterations
        """
        cdef PetscReal rval = 0
        cdef PetscInt  ival = 0
        CHKERR( LMEGetTolerances(self.lme, &rval, &ival) )
        return (toReal(rval), toInt(ival))

    def setTolerances(self, tol: float | None = None, max_it: int | None = None) -> None:
        """
        Set the tolerance and maximum iteration count.

        Logically collective.

        Set the tolerance and maximum iteration count used by the
        default LME convergence tests.

        Parameters
        ----------
        tol
            The convergence tolerance.
        max_it
            The maximum number of iterations
        """
        cdef PetscReal rval = PETSC_DEFAULT
        cdef PetscInt  ival = PETSC_DEFAULT
        if tol    is not None: rval = asReal(tol)
        if max_it is not None: ival = asInt(max_it)
        CHKERR( LMESetTolerances(self.lme, rval, ival) )

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
        CHKERR( LMEGetDimensions(self.lme, &ival) )
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
        CHKERR( LMESetDimensions(self.lme, ival) )

    def getBV(self) -> BV:
        """
        Get the basis vector object associated to the LME object.

        Not collective.

        Returns
        -------
        BV
            The basis vectors context.
        """
        cdef BV bv = BV()
        CHKERR( LMEGetBV(self.lme, &bv.bv) )
        CHKERR( PetscINCREF(bv.obj) )
        return bv

    def setBV(self, BV bv) -> None:
        """
        Set a basis vector object to the LME object.

        Collective.

        Parameters
        ----------
        bv
            The basis vectors context.
        """
        CHKERR( LMESetBV(self.lme, bv.bv) )

    def setMonitor(
        self,
        monitor: LMEMonitorFunction | None,
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
            CHKERR( LMEMonitorSet(self.lme, LME_Monitor, NULL, NULL) )
        if args is None: args = ()
        if kargs is None: kargs = {}
        monitorlist.append((monitor, args, kargs))

    def getMonitor(self) -> LMEMonitorFunction:
        """
        Get the list of monitor functions.
        """
        return self.get_attr('__monitor__')

    def cancelMonitor(self) -> None:
        """
        Clear all monitors for an `LME` object.
        """
        CHKERR( LMEMonitorCancel(self.lme) )
        self.set_attr('__monitor__', None)

    def setUp(self) -> None:
        """
        Set up all the internal necessary data structures.

        Collective.

        Set up all the internal data structures necessary for the
        execution of the eigensolver.
        """
        CHKERR( LMESetUp(self.lme) )

    def solve(self) -> None:
        """
        Solve the linear matrix equation.

        Collective.
        """
        CHKERR( LMESolve(self.lme) )

    def getIterationNumber(self) -> int:
        """
        Get the current iteration number.

        Not collective.

        If the call to `solve()` is complete, then it returns the number of
        iterations carried out by the solution method.

        Returns
        -------
        int
            Iteration number.
        """
        cdef PetscInt ival = 0
        CHKERR( LMEGetIterationNumber(self.lme, &ival) )
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
        cdef SlepcLMEConvergedReason val = LME_CONVERGED_ITERATING
        CHKERR( LMEGetConvergedReason(self.lme, &val) )
        return val

    def setErrorIfNotConverged(self, flg: bool = True) -> None:
        """
        Set `solve()` to generate an error if the solver has not converged.

        Logically collective.

        Parameters
        ----------
        flg
            True indicates you want the error generated.
        """
        cdef PetscBool tval = flg
        CHKERR( LMESetErrorIfNotConverged(self.lme, tval) )

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
        CHKERR( LMEGetErrorIfNotConverged(self.lme, &tval) )
        return toBool(tval)

    #

    property tol:
        """The tolerance value used by the LME convergence tests."""
        def __get__(self) -> float:
            return self.getTolerances()[0]
        def __set__(self, value):
            self.setTolerances(tol=value)

    property max_it:
        """The maximum iteration count used by the LME convergence tests."""
        def __get__(self) -> int:
            return self.getTolerances()[1]
        def __set__(self, value):
            self.setTolerances(max_it=value)

    property fn:
        """The math function (FN) object associated to the LME object."""
        def __get__(self) -> FN:
            return self.getFN()
        def __set__(self, value):
            self.setFN(value)

    property bv:
        """The basis vectors (BV) object associated to the LME object."""
        def __get__(self) -> BV:
            return self.getBV()
        def __set__(self, value):
            self.setBV(value)

# -----------------------------------------------------------------------------

del LMEType
del LMEConvergedReason
del LMEProblemType

# -----------------------------------------------------------------------------
