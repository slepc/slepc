# -----------------------------------------------------------------------------

class LMEType(object):
    """
    LME type.

    - `KRYLOV`:  Restarted Krylov solver.

    See Also
    --------
    slepc.LMEType
    """
    KRYLOV   = S_(LMEKRYLOV)

class LMEConvergedReason(object):
    """
    LME convergence reasons.

    - `CONVERGED_TOL`:       All eigenpairs converged to requested tolerance.
    - `DIVERGED_ITS`:        Maximum number of iterations exceeded.
    - `DIVERGED_BREAKDOWN`:  Solver failed due to breakdown.
    - `CONVERGED_ITERATING`: Iteration not finished yet.

    See Also
    --------
    slepc.LMEConvergedReason
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

    See Also
    --------
    slepc.LMEProblemType
    """
    LYAPUNOV        =  LME_LYAPUNOV
    SYLVESTER       =  LME_SYLVESTER
    GEN_LYAPUNOV    =  LME_GEN_LYAPUNOV
    GEN_SYLVESTER   =  LME_GEN_SYLVESTER
    DT_LYAPUNOV     =  LME_DT_LYAPUNOV
    STEIN           =  LME_STEIN

# -----------------------------------------------------------------------------

cdef class LME(Object):

    """
    Linear Matrix Equation.

    Linear Matrix Equation (`LME`) is the object provided by slepc4py for
    solving linear matrix equations such as Lyapunov or Sylvester where
    the solution has low rank.
    """

    Type            = LMEType
    ConvergedReason = LMEConvergedReason
    ProblemType     = LMEProblemType

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

        See Also
        --------
        slepc.LMEView
        """
        cdef PetscViewer vwr = def_Viewer(viewer)
        CHKERR( LMEView(self.lme, vwr) )

    def destroy(self) -> Self:
        """
        Destroy the LME object.

        Collective.

        See Also
        --------
        slepc.LMEDestroy
        """
        CHKERR( LMEDestroy(&self.lme) )
        self.lme = NULL
        return self

    def reset(self) -> None:
        """
        Reset the LME object.

        Collective.

        See Also
        --------
        slepc.LMEReset
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

        See Also
        --------
        slepc.LMECreate
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

        Notes
        -----
        The default is ``KRYLOV``. Normally, it is best to use
        `setFromOptions()` and then set the LME type from the options
        database rather than by using this routine. Using the options
        database provides the user with maximum flexibility in
        evaluating the different available methods.

        See Also
        --------
        getType, slepc.LMESetType
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

        See Also
        --------
        setType, slepc.LMEGetType
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

        See Also
        --------
        getProblemType, slepc.LMESetProblemType
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

        See Also
        --------
        setProblemType, slepc.LMEGetProblemType
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

        Notes
        -----
        The matrix equation takes the general form :math:`AXE+DXB=C`, where
        matrix :math:`C` is not provided here but with `setRHS()`. Not all
        four matrices must be passed.

        This must be called before `setUp()`. If called again after `setUp()`
        then the `LME` object is reset.

        See Also
        --------
        getCoefficients, solve, setRHS, setUp, slepc.LMESetCoefficients
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
        A: petsc4py.PETSc.Mat
            First coefficient matrix.
        B: petsc4py.PETSc.Mat
            Second coefficient matrix, if available.
        D: petsc4py.PETSc.Mat
            Third coefficient matrix, if available.
        E: petsc4py.PETSc.Mat
            Fourth coefficient matrix, if available.

        See Also
        --------
        setCoefficients, slepc.LMEGetCoefficients
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

        Parameters
        ----------
        C
            The right-hand side matrix

        Notes
        -----
        The matrix equation takes the general form :math:`AXE+DXB=C`, where
        matrix :math:`C` is given with this function. It must be a low-rank
        matrix of type `petsc4py.PETSc.Mat.Type.LRC`, that is,
        :math:`C = UDV^*` where :math:`D` is diagonal of order :math:`k`,
        and :math:`U,V` are dense tall-skinny matrices with :math:`k` columns.
        No sparse matrix must be provided when creating the ``LRC`` matrix.

        In equation types that require :math:`C` to be symmetric, such as
        Lyapunov, ``C`` must be created with :math:`V=U`.

        See Also
        --------
        getRHS, setSolution, slepc.LMESetRHS
        """
        CHKERR( LMESetRHS(self.lme, C.mat) )

    def getRHS(self) -> Mat:
        """
        Get the right-hand side of the matrix equation.

        Collective.

        Returns
        -------
        C: petsc4py.PETSc.Mat
            The low-rank matrix.

        See Also
        --------
        setRHS, slepc.LMEGetRHS
        """
        cdef Mat C = Mat()
        CHKERR( LMEGetRHS(self.lme, &C.mat) )
        CHKERR( PetscINCREF(C.obj) )
        return C

    def setSolution(self, Mat X = None) -> None:
        """
        Set the placeholder for the solution of the matrix equation.

        Collective.

        Parameters
        ----------
        X
            The solution matrix

        Notes
        -----
        The matrix equation takes the general form :math:`AXE+DXB=C`, where
        the solution matrix is of low rank and is written in factored form
        :math:`X = UDV^*`. This function provides a matrix object of type
        `petsc4py.PETSc.Mat.Type.LRC` that stores :math:`U,V` and
        (optionally) :math:`D`. These factors will be computed during `solve()`.

        In equation types whose solution :math:`X` is symmetric, such as
        Lyapunov, ``X`` must be created with :math:`V=U`.

        If the user provides ``X`` with this function, then the solver will
        return a solution with rank at most the number of columns of :math:`U`.
        Alternatively, it is possible to let the solver choose the rank of the
        solution, by passing ``None`` and then calling `getSolution()` after
        `solve()`.

        See Also
        --------
        solve, setRHS, getSolution, slepc.LMESetSolution
        """
        cdef PetscMat Xmat = X.mat if X is not None else <PetscMat>NULL
        CHKERR( LMESetSolution(self.lme, Xmat) )

    def getSolution(self) -> Mat:
        """
        Get the solution of the matrix equation.

        Collective.

        Returns
        -------
        X: petsc4py.PETSc.Mat
            The low-rank matrix.

        Notes
        -----
        If called after `solve()`, ``X`` will contain the solution of the
        equation.

        The matrix ``X`` may have been passed by the user via `setSolution()`,
        although this is not required.

        See Also
        --------
        solve, setSolution, slepc.LMEGetSolution
        """
        cdef Mat X = Mat()
        CHKERR( LMEGetSolution(self.lme, &X.mat) )
        CHKERR( PetscINCREF(X.obj) )
        return X

    def getErrorEstimate(self) -> float:
        """
        Get the error estimate obtained during the solve.

        Not collective.

        Returns
        -------
        float
            The error estimate.

        Notes
        -----
        This is the error estimated internally by the solver. The actual
        error bound can be computed with `computeError()`. Note that some
        solvers may not be able to provide an error estimate.

        See Also
        --------
        computeError, slepc.LMEGetErrorEstimate
        """
        cdef PetscReal rval = 0
        CHKERR( LMEGetErrorEstimate(self.lme, &rval) )
        return toReal(rval)

    def computeError(self) -> float:
        """
        Compute the error associated with the last equation solved.

        Collective.

        Returns
        -------
        float
            The error.

        Notes
        -----
        The error is based on the residual norm.

        This function is not scalable (in terms of memory or parallel
        communication), so it should not be called except in the case of
        small problem size. For large equations, use `getErrorEstimate()`.

        See Also
        --------
        getErrorEstimate, slepc.LMEComputeError
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

        See Also
        --------
        setOptionsPrefix, appendOptionsPrefix, slepc.LMEGetOptionsPrefix
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

        Notes
        -----
        A hyphen (-) must NOT be given at the beginning of the prefix
        name.  The first character of all runtime options is
        AUTOMATICALLY the hyphen.

        For example, to distinguish between the runtime options for
        two different LME contexts, one could call::

            L1.setOptionsPrefix("lme1_")
            L2.setOptionsPrefix("lme2_")

        See Also
        --------
        appendOptionsPrefix, getOptionsPrefix, slepc.LMEGetOptionsPrefix
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

        See Also
        --------
        setOptionsPrefix, getOptionsPrefix, slepc.LMEAppendOptionsPrefix
        """
        cdef const char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( LMEAppendOptionsPrefix(self.lme, cval) )

    def setFromOptions(self) -> None:
        """
        Set LME options from the options database.

        Collective.

        Notes
        -----
        To see all options, run your program with the ``-help`` option.

        This routine must be called before `setUp()` if the user is to be
        allowed to set the solver type.

        See Also
        --------
        setOptionsPrefix, slepc.LMESetFromOptions
        """
        CHKERR( LMESetFromOptions(self.lme) )

    def getTolerances(self) -> tuple[float, int]:
        """
        Get the tolerance and maximum iteration count.

        Not collective.

        Returns
        -------
        tol: float
            The convergence tolerance.
        max_it: int
            The maximum number of iterations.

        See Also
        --------
        setTolerances, slepc.LMEGetTolerances
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
            The maximum number of iterations.

        See Also
        --------
        getTolerances, slepc.LMESetTolerances
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

        See Also
        --------
        setDimensions, slepc.LMEGetDimensions
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

        See Also
        --------
        getDimensions, slepc.LMESetDimensions
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

        See Also
        --------
        setBV, slepc.LMEGetBV
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

        See Also
        --------
        getBV, slepc.LMESetBV
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

        See Also
        --------
        getMonitor, cancelMonitor, slepc.LMEMonitorSet
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

        Not collective.

        Returns
        -------
        LMEMonitorFunction
            The list of monitor functions.
        """
        return self.get_attr('__monitor__')

    def cancelMonitor(self) -> None:
        """
        Clear all monitors for an `LME` object.

        Logically collective.

        See Also
        --------
        slepc.LMEMonitorCancel
        """
        CHKERR( LMEMonitorCancel(self.lme) )
        self.set_attr('__monitor__', None)

    def setUp(self) -> None:
        """
        Set up all the internal necessary data structures.

        Collective.

        Set up all the internal data structures necessary for the
        execution of the eigensolver.

        See Also
        --------
        solve, slepc.LMESetUp
        """
        CHKERR( LMESetUp(self.lme) )

    def solve(self) -> None:
        """
        Solve the linear matrix equation.

        Collective.

        Notes
        -----
        The matrix coefficients are specified with `setCoefficients()`.
        The right-hand side is specified with `setRHS()`. The placeholder
        for the solution is specified with `setSolution()`.
        See Also
        --------
        setCoefficients, setRHS, setSolution, slepc.LMESolve
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

        See Also
        --------
        getConvergedReason, slepc.LMEGetIterationNumber
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

        See Also
        --------
        setTolerances, solve, setErrorIfNotConverged, slepc.LMEGetConvergedReason
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
            ``True`` indicates you want the error generated.

        Notes
        -----
        Normally SLEPc continues if the solver fails to converge, you can
        call `getConvergedReason()` after a `solve()` to determine if it
        has converged.

        See Also
        --------
        getConvergedReason, solve, slepc.LMESetErrorIfNotConverged
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
            ``True`` indicates you want the error generated.

        See Also
        --------
        setErrorIfNotConverged, slepc.LMEGetErrorIfNotConverged
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
        """The math function (`FN`) object associated to the LME object."""
        def __get__(self) -> FN:
            return self.getFN()
        def __set__(self, value):
            self.setFN(value)

    property bv:
        """The basis vectors (`BV`) object associated to the LME object."""
        def __get__(self) -> BV:
            return self.getBV()
        def __set__(self, value):
            self.setBV(value)

# -----------------------------------------------------------------------------

del LMEType
del LMEConvergedReason
del LMEProblemType

# -----------------------------------------------------------------------------
