# -----------------------------------------------------------------------------

class NEPType(object):
    """
    NEP type.

    - `RII`:      Residual inverse iteration.
    - `SLP`:      Successive linear problems.
    - `NARNOLDI`: Nonlinear Arnoldi.
    - `NLEIGS`:   Fully rational Krylov method for nonlinear eigenproblems.
    - `CISS`:     Contour integral spectrum slice.
    - `INTERPOL`: Polynomial interpolation.

    See Also
    --------
    slepc.NEPType
    """
    RII      = S_(NEPRII)
    SLP      = S_(NEPSLP)
    NARNOLDI = S_(NEPNARNOLDI)
    NLEIGS   = S_(NEPNLEIGS)
    CISS     = S_(NEPCISS)
    INTERPOL = S_(NEPINTERPOL)

class NEPProblemType(object):
    """
    NEP problem type.

    - `GENERAL`:  General nonlinear eigenproblem.
    - `RATIONAL`: NEP defined in split form with all :math:`f_i` rational.

    See Also
    --------
    slepc.NEPProblemType
    """
    GENERAL  = NEP_GENERAL
    RATIONAL = NEP_RATIONAL

class NEPErrorType(object):
    """
    NEP error type to assess accuracy of computed solutions.

    - `ABSOLUTE`: Absolute error.
    - `RELATIVE`: Relative error.
    - `BACKWARD`: Backward error.

    See Also
    --------
    slepc.NEPErrorType
    """
    ABSOLUTE = NEP_ERROR_ABSOLUTE
    RELATIVE = NEP_ERROR_RELATIVE
    BACKWARD = NEP_ERROR_BACKWARD

class NEPWhich(object):
    """
    NEP desired part of spectrum.

    - `LARGEST_MAGNITUDE`:  Largest magnitude (default).
    - `SMALLEST_MAGNITUDE`: Smallest magnitude.
    - `LARGEST_REAL`:       Largest real parts.
    - `SMALLEST_REAL`:      Smallest real parts.
    - `LARGEST_IMAGINARY`:  Largest imaginary parts in magnitude.
    - `SMALLEST_IMAGINARY`: Smallest imaginary parts in magnitude.
    - `TARGET_MAGNITUDE`:   Closest to target (in magnitude).
    - `TARGET_REAL`:        Real part closest to target.
    - `TARGET_IMAGINARY`:   Imaginary part closest to target.
    - `ALL`:                All eigenvalues in a region.
    - `USER`:               User defined selection.

    See Also
    --------
    slepc.NEPWhich
    """
    LARGEST_MAGNITUDE  = NEP_LARGEST_MAGNITUDE
    SMALLEST_MAGNITUDE = NEP_SMALLEST_MAGNITUDE
    LARGEST_REAL       = NEP_LARGEST_REAL
    SMALLEST_REAL      = NEP_SMALLEST_REAL
    LARGEST_IMAGINARY  = NEP_LARGEST_IMAGINARY
    SMALLEST_IMAGINARY = NEP_SMALLEST_IMAGINARY
    TARGET_MAGNITUDE   = NEP_TARGET_MAGNITUDE
    TARGET_REAL        = NEP_TARGET_REAL
    TARGET_IMAGINARY   = NEP_TARGET_IMAGINARY
    ALL                = NEP_ALL
    USER               = NEP_WHICH_USER

class NEPConvergedReason(object):
    """
    NEP convergence reasons.

    - `CONVERGED_TOL`: All eigenpairs converged to requested tolerance.
    - `CONVERGED_USER`: User-defined convergence criterion satisfied.
    - `DIVERGED_ITS`: Maximum number of iterations exceeded.
    - `DIVERGED_BREAKDOWN`: Solver failed due to breakdown.
    - `DIVERGED_LINEAR_SOLVE`: Inner linear solve failed.
    - `DIVERGED_SUBSPACE_EXHAUSTED`: Run out of space for the basis in an
      unrestarted solver.
    - `CONVERGED_ITERATING`: Iteration not finished yet.

    See Also
    --------
    slepc.NEPConvergedReason
    """
    CONVERGED_TOL               = NEP_CONVERGED_TOL
    CONVERGED_USER              = NEP_CONVERGED_USER
    DIVERGED_ITS                = NEP_DIVERGED_ITS
    DIVERGED_BREAKDOWN          = NEP_DIVERGED_BREAKDOWN
    DIVERGED_LINEAR_SOLVE       = NEP_DIVERGED_LINEAR_SOLVE
    DIVERGED_SUBSPACE_EXHAUSTED = NEP_DIVERGED_SUBSPACE_EXHAUSTED
    CONVERGED_ITERATING         = NEP_CONVERGED_ITERATING
    ITERATING                   = NEP_CONVERGED_ITERATING

class NEPRefine(object):
    """
    NEP refinement strategy.

    - `NONE`:     No refinement.
    - `SIMPLE`:   Refine eigenpairs one by one.
    - `MULTIPLE`: Refine all eigenpairs simultaneously (invariant pair).

    See Also
    --------
    slepc.NEPRefine
    """
    NONE     = NEP_REFINE_NONE
    SIMPLE   = NEP_REFINE_SIMPLE
    MULTIPLE = NEP_REFINE_MULTIPLE

class NEPRefineScheme(object):
    """
    NEP scheme for solving linear systems during iterative refinement.

    - `SCHUR`:    Schur complement.
    - `MBE`:      Mixed block elimination.
    - `EXPLICIT`: Build the explicit matrix.

    See Also
    --------
    slepc.NEPRefineScheme
    """
    SCHUR    = NEP_REFINE_SCHEME_SCHUR
    MBE      = NEP_REFINE_SCHEME_MBE
    EXPLICIT = NEP_REFINE_SCHEME_EXPLICIT

class NEPConv(object):
    """
    NEP convergence test.

    - `ABS`:  Absolute convergence test.
    - `REL`:  Convergence test relative to the eigenvalue.
    - `NORM`: Convergence test relative to the matrix norms.
    - `USER`: User-defined convergence test.

    See Also
    --------
    slepc.NEPConv
    """
    ABS  = NEP_CONV_ABS
    REL  = NEP_CONV_REL
    NORM = NEP_CONV_NORM
    USER = NEP_CONV_USER

class NEPStop(object):
    """
    NEP stopping test.

    - `BASIC`: Default stopping test.
    - `USER`:  User-defined stopping test.

    See Also
    --------
    slepc.NEPStop
    """
    BASIC = NEP_STOP_BASIC
    USER  = NEP_STOP_USER

class NEPCISSExtraction(object):
    """
    NEP CISS extraction technique.

    - `RITZ`:   Ritz extraction.
    - `HANKEL`: Extraction via Hankel eigenproblem.
    - `CAA`:    Communication-avoiding Arnoldi.

    See Also
    --------
    slepc.NEPCISSExtraction
    """
    RITZ   = NEP_CISS_EXTRACTION_RITZ
    HANKEL = NEP_CISS_EXTRACTION_HANKEL
    CAA    = NEP_CISS_EXTRACTION_CAA

# -----------------------------------------------------------------------------

cdef class NEP(Object):

    """
    Nonlinear Eigenvalue Problem Solver.

    The Nonlinear Eigenvalue Problem (`NEP`) solver is the object provided
    by slepc4py for specifying an eigenvalue problem that is nonlinear with
    respect to the eigenvalue (not the eigenvector). This is intended for
    general nonlinear problems (rather than polynomial eigenproblems)
    described as :math:`T(\lambda) x=0`.
    """

    Type            = NEPType
    ProblemType     = NEPProblemType
    ErrorType       = NEPErrorType
    Which           = NEPWhich
    ConvergedReason = NEPConvergedReason
    Refine          = NEPRefine
    RefineScheme    = NEPRefineScheme
    Conv            = NEPConv
    Stop            = NEPStop

    CISSExtraction  = NEPCISSExtraction

    def __cinit__(self):
        self.obj = <PetscObject*> &self.nep
        self.nep = NULL

    def view(self, Viewer viewer=None) -> None:
        """
        Print the NEP data structure.

        Collective.

        Parameters
        ----------
        viewer
            Visualization context; if not provided, the standard
            output is used.

        See Also
        --------
        slepc.NEPView
        """
        cdef PetscViewer vwr = def_Viewer(viewer)
        CHKERR( NEPView(self.nep, vwr) )

    def destroy(self) -> Self:
        """
        Destroy the NEP object.

        Collective.

        See Also
        --------
        slepc.NEPDestroy
        """
        CHKERR( NEPDestroy(&self.nep) )
        self.nep = NULL
        return self

    def reset(self) -> None:
        """
        Reset the NEP object.

        Collective.

        See Also
        --------
        slepc.NEPReset
        """
        CHKERR( NEPReset(self.nep) )

    def create(self, comm: Comm | None = None) -> Self:
        """
        Create the NEP object.

        Collective.

        Parameters
        ----------
        comm
            MPI communicator. If not provided, it defaults to all processes.

        See Also
        --------
        slepc.NEPCreate
        """
        cdef MPI_Comm ccomm = def_Comm(comm, SLEPC_COMM_DEFAULT())
        cdef SlepcNEP newnep = NULL
        CHKERR( NEPCreate(ccomm, &newnep) )
        CHKERR( SlepcCLEAR(self.obj) ); self.nep = newnep
        return self

    def setType(self, nep_type: Type | str) -> None:
        """
        Set the particular solver to be used in the NEP object.

        Logically collective.

        Parameters
        ----------
        nep_type
            The solver to be used.

        Notes
        -----
        The default is `RII`. Normally, it is best to use
        `setFromOptions()` and then set the NEP type from the options
        database rather than by using this routine. Using the options
        database provides the user with maximum flexibility in
        evaluating the different available methods.

        See Also
        --------
        getType, slepc.NEPSetType
        """
        cdef SlepcNEPType cval = NULL
        nep_type = str2bytes(nep_type, &cval)
        CHKERR( NEPSetType(self.nep, cval) )

    def getType(self) -> str:
        """
        Get the NEP type of this object.

        Not collective.

        Returns
        -------
        str
            The solver currently being used.

        See Also
        --------
        setType, slepc.NEPGetType
        """
        cdef SlepcNEPType nep_type = NULL
        CHKERR( NEPGetType(self.nep, &nep_type) )
        return bytes2str(nep_type)

    def getOptionsPrefix(self) -> str:
        """
        Get the prefix used for searching for all NEP options in the database.

        Not collective.

        Returns
        -------
        str
            The prefix string set for this NEP object.

        See Also
        --------
        setOptionsPrefix, appendOptionsPrefix, slepc.NEPGetOptionsPrefix
        """
        cdef const char *prefix = NULL
        CHKERR( NEPGetOptionsPrefix(self.nep, &prefix) )
        return bytes2str(prefix)

    def setOptionsPrefix(self, prefix: str | None = None) -> None:
        """
        Set the prefix used for searching for all NEP options in the database.

        Logically collective.

        Parameters
        ----------
        prefix
            The prefix string to prepend to all NEP option requests.

        Notes
        -----
        A hyphen (-) must NOT be given at the beginning of the prefix
        name.  The first character of all runtime options is
        AUTOMATICALLY the hyphen.

        For example, to distinguish between the runtime options for
        two different NEP contexts, one could call::

            N1.setOptionsPrefix("nep1_")
            N2.setOptionsPrefix("nep2_")

        See Also
        --------
        appendOptionsPrefix, getOptionsPrefix, slepc.NEPGetOptionsPrefix
        """
        cdef const char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( NEPSetOptionsPrefix(self.nep, cval) )

    def appendOptionsPrefix(self, prefix: str | None = None) -> None:
        """
        Append to the prefix used for searching for all NEP options in the database.

        Logically collective.

        Parameters
        ----------
        prefix
            The prefix string to prepend to all NEP option requests.

        See Also
        --------
        setOptionsPrefix, getOptionsPrefix, slepc.NEPAppendOptionsPrefix
        """
        cdef const char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( NEPAppendOptionsPrefix(self.nep, cval) )

    def setFromOptions(self) -> None:
        """
        Set NEP options from the options database.

        Collective.

        Notes
        -----
        To see all options, run your program with the ``-help`` option.

        This routine must be called before `setUp()` if the user is to be
        allowed to set the solver type.

        See Also
        --------
        setOptionsPrefix, slepc.NEPSetFromOptions
        """
        CHKERR( NEPSetFromOptions(self.nep) )

    def getProblemType(self) -> ProblemType:
        """
        Get the problem type from the `NEP` object.

        Not collective.

        Returns
        -------
        ProblemType
            The problem type that was previously set.

        See Also
        --------
        setProblemType, slepc.NEPGetProblemType
        """
        cdef SlepcNEPProblemType val = NEP_GENERAL
        CHKERR( NEPGetProblemType(self.nep, &val) )
        return val

    def setProblemType(self, problem_type: ProblemType) -> None:
        """
        Set the type of the eigenvalue problem.

        Logically collective.

        Parameters
        ----------
        problem_type
            The problem type to be set.

        Notes
        -----
        This function is used to provide a hint to the `NEP` solver
        to exploit certain properties of the nonlinear eigenproblem.
        This hint may be used or not, depending on the solver. By
        default, no particular structure is assumed.

        See Also
        --------
        getProblemType, slepc.NEPSetProblemType
        """
        cdef SlepcNEPProblemType val = problem_type
        CHKERR( NEPSetProblemType(self.nep, val) )

    def getWhichEigenpairs(self) -> Which:
        """
        Get which portion of the spectrum is to be sought.

        Not collective.

        Returns
        -------
        Which
            The portion of the spectrum to be sought by the solver.

        See Also
        --------
        setWhichEigenpairs, slepc.NEPGetWhichEigenpairs
        """
        cdef SlepcNEPWhich val = NEP_LARGEST_MAGNITUDE
        CHKERR( NEPGetWhichEigenpairs(self.nep, &val) )
        return val

    def setWhichEigenpairs(self, which: Which) -> None:
        """
        Set which portion of the spectrum is to be sought.

        Logically collective.

        Parameters
        ----------
        which
            The portion of the spectrum to be sought by the solver.

        Notes
        -----
        Not all eigensolvers implemented in NEP account for all the
        possible values. Also, some values make sense only for certain
        types of problems. If SLEPc is compiled for real numbers
        `NEP.Which.LARGEST_IMAGINARY` and
        `NEP.Which.SMALLEST_IMAGINARY` use the absolute value of the
        imaginary part for eigenvalue selection.

        The target is a scalar value provided with `setTarget()`.

        The criterion `NEP.Which.TARGET_IMAGINARY` is available only
        in case PETSc and SLEPc have been built with complex scalars.

        `NEP.Which.ALL` is intended for use in the context of the
        `PEP.Type.CISS` solver for computing all eigenvalues in a region.

        See Also
        --------
        getWhichEigenpairs, setTarget, slepc.PEPSetWhichEigenpairs
        """
        cdef SlepcNEPWhich val = which
        CHKERR( NEPSetWhichEigenpairs(self.nep, val) )

    def getTarget(self) -> Scalar:
        """
        Get the value of the target.

        Not collective.

        Returns
        -------
        Scalar
            The value of the target.

        Notes
        -----
        If the target was not set by the user, then zero is returned.

        See Also
        --------
        setTarget, slepc.NEPGetTarget
        """
        cdef PetscScalar sval = 0
        CHKERR( NEPGetTarget(self.nep, &sval) )
        return toScalar(sval)

    def setTarget(self, target: Scalar) -> None:
        """
        Set the value of the target.

        Logically collective.

        Parameters
        ----------
        target
            The value of the target.

        Notes
        -----
        The target is a scalar value used to determine the portion of
        the spectrum of interest. It is used in combination with
        `setWhichEigenpairs()`.

        When PETSc is built with real scalars, it is not possible to
        specify a complex target.

        See Also
        --------
        getTarget, setWhichEigenpairs, slepc.NEPSetTarget
        """
        cdef PetscScalar sval = asScalar(target)
        CHKERR( NEPSetTarget(self.nep, sval) )

    def getTolerances(self) -> tuple[float, int]:
        """
        Get the tolerance and maximum iteration count.

        Not collective.

        Get the tolerance and maximum iteration count used by the default
        NEP convergence tests.

        Returns
        -------
        tol: float
            The convergence tolerance.
        maxit: int
            The maximum number of iterations.

        See Also
        --------
        setTolerances, slepc.NEPGetTolerances
        """
        cdef PetscReal rval = 0
        cdef PetscInt  ival = 0
        CHKERR( NEPGetTolerances(self.nep, &rval, &ival) )
        return (toReal(rval), toInt(ival))

    def setTolerances(self, tol: float | None = None, maxit: int | None = None) -> None:
        """
        Set the tolerance and max. iteration count used in convergence tests.

        Logically collective.

        Parameters
        ----------
        tol
            The convergence tolerance.
        maxit
            The maximum number of iterations.

        Notes
        -----
        Use `DETERMINE` for ``max_it`` to assign a reasonably good value,
        which is dependent on the solution method.

        See Also
        --------
        getTolerances, slepc.NEPSetTolerances
        """
        cdef PetscReal rval = PETSC_CURRENT
        cdef PetscInt  ival = PETSC_CURRENT
        if tol   is not None: rval = asReal(tol)
        if maxit is not None: ival = asInt(maxit)
        CHKERR( NEPSetTolerances(self.nep, rval, ival) )

    def getConvergenceTest(self) -> Conv:
        """
        Get the method used to compute the error estimate used in the convergence test.

        Not collective.

        Returns
        -------
        Conv
            The method used to compute the error estimate
            used in the convergence test.

        See Also
        --------
        setConvergenceTest, slepc.NEPGetConvergenceTest
        """
        cdef SlepcNEPConv conv = NEP_CONV_REL
        CHKERR( NEPGetConvergenceTest(self.nep, &conv) )
        return conv

    def setConvergenceTest(self, conv: Conv) -> None:
        """
        Set how to compute the error estimate used in the convergence test.

        Logically collective.

        Parameters
        ----------
        conv
            The method used to compute the error estimate
            used in the convergence test.

        See Also
        --------
        getConvergenceTest, slepc.NEPSetConvergenceTest
        """
        cdef SlepcNEPConv tconv = conv
        CHKERR( NEPSetConvergenceTest(self.nep, tconv) )

    def getRefine(self) -> tuple[Refine, int, float, int, RefineScheme]:
        """
        Get the refinement strategy used by the NEP object.

        Not collective.

        Returns
        -------
        ref: Refine
            The refinement type.
        npart: int
            The number of partitions of the communicator.
        tol: float
            The convergence tolerance.
        its: int
            The maximum number of refinement iterations.
        scheme: RefineScheme
            Scheme for solving linear systems.

        See Also
        --------
        setRefine, slepc.NEPGetRefine
        """
        cdef SlepcNEPRefine ref = NEP_REFINE_NONE
        cdef PetscInt npart = 1
        cdef PetscReal tol = PETSC_DEFAULT
        cdef PetscInt its = PETSC_DEFAULT
        cdef SlepcNEPRefineScheme scheme = NEP_REFINE_SCHEME_MBE
        CHKERR( NEPGetRefine(self.nep, &ref, &npart, &tol, &its, &scheme) )
        return (ref, toInt(npart), toReal(tol), toInt(its), scheme)

    def setRefine(
        self,
        ref: Refine,
        npart: int | None = None,
        tol: float | None = None,
        its: int | None = None,
        scheme: RefineScheme | None = None,
    ) -> None:
        """
        Set the refinement strategy used by the NEP object.

        Logically collective.

        Set the refinement strategy used by the NEP object and the associated
        parameters.

        Parameters
        ----------
        ref
            The refinement type.
        npart
            The number of partitions of the communicator.
        tol
            The convergence tolerance.
        its
            The maximum number of refinement iterations.
        scheme
            Scheme for solving linear systems.

        See Also
        --------
        getRefine, slepc.NEPSetRefine
        """
        cdef SlepcNEPRefine tref = ref
        cdef PetscInt tnpart = PETSC_CURRENT
        cdef PetscReal ttol = PETSC_CURRENT
        cdef PetscInt tits = PETSC_CURRENT
        cdef SlepcNEPRefineScheme tscheme = NEP_REFINE_SCHEME_MBE
        if npart is not None: tnpart = asInt(npart)
        if tol is not None: ttol = asReal(tol)
        if its is not None: tits = asInt(its)
        if scheme is not None: tscheme = scheme
        CHKERR( NEPSetRefine(self.nep, tref, tnpart, ttol, tits, tscheme) )

    def getRefineKSP(self) -> KSP:
        """
        Get the ``KSP`` object used by the eigensolver in the refinement phase.

        Collective.

        Returns
        -------
        `petsc4py.PETSc.KSP`
            The linear solver object.

        See Also
        --------
        setRefine, slepc.NEPRefineGetKSP
        """
        cdef KSP ksp = KSP()
        CHKERR( NEPRefineGetKSP(self.nep, &ksp.ksp) )
        CHKERR( PetscINCREF(ksp.obj) )
        return ksp

    def getTrackAll(self) -> bool:
        """
        Get the flag indicating whether all residual norms must be computed.

        Not collective.

        Returns
        -------
        bool
            Whether the solver computes all residuals or not.

        See Also
        --------
        setTrackAll, slepc.NEPGetTrackAll
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( NEPGetTrackAll(self.nep, &tval) )
        return toBool(tval)

    def setTrackAll(self, trackall: bool) -> None:
        """
        Set if the solver must compute the residual of all approximate eigenpairs.

        Logically collective.

        Parameters
        ----------
        trackall
            Whether to compute all residuals or not.

        See Also
        --------
        getTrackAll, slepc.NEPSetTrackAll
        """
        cdef PetscBool tval = trackall
        CHKERR( NEPSetTrackAll(self.nep, tval) )

    def getDimensions(self) -> tuple[int, int, int]:
        """
        Get the number of eigenvalues to compute.

        Not collective.

        Get the number of eigenvalues to compute, and the dimension of the
        subspace.

        Returns
        -------
        nev: int
            Number of eigenvalues to compute.
        ncv: int
            Maximum dimension of the subspace to be used by the solver.
        mpd: int
            Maximum dimension allowed for the projected problem.

        See Also
        --------
        setDimensions, slepc.NEPGetDimensions
        """
        cdef PetscInt ival1 = 0
        cdef PetscInt ival2 = 0
        cdef PetscInt ival3 = 0
        CHKERR( NEPGetDimensions(self.nep, &ival1, &ival2, &ival3) )
        return (toInt(ival1), toInt(ival2), toInt(ival3))

    def setDimensions(
        self,
        nev: int | None = None,
        ncv: int | None = None,
        mpd: int | None = None,
    ) -> None:
        """
        Set the number of eigenvalues to compute.

        Logically collective.

        Set the number of eigenvalues to compute and the dimension of the
        subspace.

        Parameters
        ----------
        nev
            Number of eigenvalues to compute.
        ncv
            Maximum dimension of the subspace to be used by the solver.
        mpd
            Maximum dimension allowed for the projected problem.

        Notes
        -----
        Use `DETERMINE` for ``ncv`` and ``mpd`` to assign a reasonably good
        value, which is dependent on the solution method.

        The parameters ``ncv`` and ``mpd`` are intimately related, so that
        the user is advised to set one of them at most. Normal usage
        is the following:

        + In cases where ``nev`` is small, the user sets ``ncv``
          (a reasonable default is 2 * ``nev``).

        + In cases where ``nev`` is large, the user sets ``mpd``.

        The value of ``ncv`` should always be between ``nev`` and (``nev`` +
        ``mpd``), typically ``ncv`` = ``nev`` + ``mpd``. If ``nev`` is not too
        large, ``mpd`` = ``nev`` is a reasonable choice, otherwise a
        smaller value should be used.

        See Also
        --------
        getDimensions, slepc.NEPSetDimensions
        """
        cdef PetscInt ival1 = PETSC_CURRENT
        cdef PetscInt ival2 = PETSC_CURRENT
        cdef PetscInt ival3 = PETSC_CURRENT
        if nev is not None: ival1 = asInt(nev)
        if ncv is not None: ival2 = asInt(ncv)
        if mpd is not None: ival3 = asInt(mpd)
        CHKERR( NEPSetDimensions(self.nep, ival1, ival2, ival3) )

    def getBV(self) -> BV:
        """
        Get the basis vectors object associated to the eigensolver.

        Not collective.

        Returns
        -------
        BV
            The basis vectors context.

        See Also
        --------
        setBV, slepc.NEPGetBV
        """
        cdef BV bv = BV()
        CHKERR( NEPGetBV(self.nep, &bv.bv) )
        CHKERR( PetscINCREF(bv.obj) )
        return bv

    def setBV(self, BV bv) -> None:
        """
        Set the basis vectors object associated to the eigensolver.

        Collective.

        Parameters
        ----------
        bv
            The basis vectors context.

        See Also
        --------
        getBV, slepc.NEPSetBV
        """
        CHKERR( NEPSetBV(self.nep, bv.bv) )

    def getRG(self) -> RG:
        """
        Get the region object associated to the eigensolver.

        Not collective.

        Returns
        -------
        RG
            The region context.

        See Also
        --------
        setRG, slepc.NEPGetRG
        """
        cdef RG rg = RG()
        CHKERR( NEPGetRG(self.nep, &rg.rg) )
        CHKERR( PetscINCREF(rg.obj) )
        return rg

    def setRG(self, RG rg) -> None:
        """
        Set a region object associated to the eigensolver.

        Collective.

        Parameters
        ----------
        rg
            The region context.

        See Also
        --------
        getRG, slepc.NEPSetRG
        """
        CHKERR( NEPSetRG(self.nep, rg.rg) )

    def getDS(self) -> DS:
        """
        Get the direct solver associated to the eigensolver.

        Not collective.

        Returns
        -------
        DS
            The direct solver context.

        See Also
        --------
        setDS, slepc.NEPGetDS
        """
        cdef DS ds = DS()
        CHKERR( NEPGetDS(self.nep, &ds.ds) )
        CHKERR( PetscINCREF(ds.obj) )
        return ds

    def setDS(self, DS ds) -> None:
        """
        Set a direct solver object associated to the eigensolver.

        Collective.

        Parameters
        ----------
        ds
            The direct solver context.

        See Also
        --------
        getDS, slepc.NEPSetDS
        """
        CHKERR( NEPSetDS(self.nep, ds.ds) )

    #

    def setInitialSpace(self, space: Vec or list[Vec]) -> None:
        """
        Set the initial space from which the eigensolver starts to iterate.

        Collective.

        Parameters
        ----------
        space
            The initial space.

        Notes
        -----
        Some solvers start to iterate on a single vector (initial vector).
        In that case, only the first vector is taken into account and the
        other vectors are ignored.

        These vectors do not persist from one `solve()` call to the other,
        so the initial space should be set every time.

        The vectors do not need to be mutually orthonormal, since they are
        explicitly orthonormalized internally.

        Common usage of this function is when the user can provide a rough
        approximation of the wanted eigenspace. Then, convergence may be faster.

        See Also
        --------
        setUp, slepc.NEPSetInitialSpace
        """
        if isinstance(space, Vec): space = [space]
        cdef PetscVec *vs = NULL
        cdef Py_ssize_t i = 0, ns = len(space)
        cdef tmp = allocate(<size_t>ns*sizeof(PetscVec),<void**>&vs)
        for i in range(ns): vs[i] = (<Vec?>space[i]).vec
        CHKERR( NEPSetInitialSpace(self.nep, <PetscInt>ns, vs) )

    #

    def setStoppingTest(
        self,
        stopping: NEPStoppingFunction | None,
        args: tuple[Any, ...] | None = None,
        kargs: dict[str, Any] | None = None,
    ) -> None:
        """
        Set a function to decide when to stop the outer iteration of the eigensolver.

        Logically collective.

        See Also
        --------
        getStoppingTest, slepc.NEPSetStoppingTestFunction
        """
        if stopping is not None:
            if args is None: args = ()
            if kargs is None: kargs = {}
            self.set_attr('__stopping__', (stopping, args, kargs))
            CHKERR( NEPSetStoppingTestFunction(self.nep, NEP_Stopping, NULL, NULL) )
        else:
            self.set_attr('__stopping__', None)
            CHKERR( NEPSetStoppingTestFunction(self.nep, NEPStoppingBasic, NULL, NULL) )

    def getStoppingTest(self) -> NEPStoppingFunction:
        """
        Get the stopping test function.

        Not collective.

        Returns
        -------
        NEPStoppingFunction
            The stopping test function.

        See Also
        --------
        setStoppingTest
        """
        return self.get_attr('__stopping__')

    def setEigenvalueComparison(
        self,
        comparison: NEPEigenvalueComparison | None,
        args: tuple[Any, ...] | None = None,
        kargs: dict[str, Any] | None = None,
    ) -> None:
        """
        Set an eigenvalue comparison function.

        Logically collective.

        Notes
        -----
        This eigenvalue comparison function is used when `setWhichEigenpairs()`
        is set to `NEP.Which.USER`.

        See Also
        --------
        getEigenvalueComparison, slepc.NEPSetEigenvalueComparison
        """
        if comparison is not None:
            if args is None: args = ()
            if kargs is None: kargs = {}
            self.set_attr('__comparison__', (comparison, args, kargs))
            ctx = self.get_attr('__comparison__')
            CHKERR( NEPSetEigenvalueComparison(self.nep, NEP_Comparison, <void*>ctx) )
        else:
            self.set_attr('__comparison__', None)
            CHKERR( NEPSetEigenvalueComparison(self.nep, NULL, NULL) )

    def getEigenvalueComparison(self) -> NEPEigenvalueComparison:
        """
        Get the eigenvalue comparison function.

        Not collective.

        Returns
        -------
        NEPEigenvalueComparison
            The eigenvalue comparison function.

        See Also
        --------
        setEigenvalueComparison
        """
        return self.get_attr('__comparison__')

    def setMonitor(
        self,
        monitor: NEPMonitorFunction | None,
        args: tuple[Any, ...] | None = None,
        kargs: dict[str, Any] | None = None,
    ) -> None:
        """
        Append a monitor function to the list of monitors.

        Logically collective.

        See Also
        --------
        getMonitor, cancelMonitor, slepc.NEPMonitorSet
        """
        if monitor is None: return
        cdef object monitorlist = self.get_attr('__monitor__')
        if monitorlist is None:
            monitorlist = []
            self.set_attr('__monitor__', monitorlist)
            CHKERR( NEPMonitorSet(self.nep, NEP_Monitor, NULL, NULL) )
        if args is None: args = ()
        if kargs is None: kargs = {}
        monitorlist.append((monitor, args, kargs))

    def getMonitor(self) -> NEPMonitorFunction:
        """
        Get the list of monitor functions.

        Not collective.

        Returns
        -------
        NEPMonitorFunction
            The list of monitor functions.

        See Also
        --------
        setMonitor
        """
        return self.get_attr('__monitor__')

    def cancelMonitor(self) -> None:
        """
        Clear all monitors for a `NEP` object.

        Logically collective.

        See Also
        --------
        slepc.NEPMonitorCancel
        """
        CHKERR( NEPMonitorCancel(self.nep) )
        self.set_attr('__monitor__', None)

    #

    def setUp(self) -> None:
        """
        Set up all the internal data structures.

        Collective.

        Notes
        -----
        Sets up all the internal data structures necessary for the execution
        of the eigensolver.

        This function need not be called explicitly in most cases,
        since `solve()` calls it. It can be useful when one wants to
        measure the set-up time separately from the solve time.

        See Also
        --------
        solve, slepc.NEPSetUp
        """
        CHKERR( NEPSetUp(self.nep) )

    def solve(self) -> None:
        """
        Solve the nonlinear eigenproblem.

        Collective.

        Notes
        -----
        `solve()` will return without generating an error regardless of
        whether all requested solutions were computed or not. Call
        `getConverged()` to get the actual number of computed solutions,
        and `getConvergedReason()` to determine if the solver converged
        or failed and why.

        See Also
        --------
        setUp, getConverged, getConvergedReason, slepc.NEPSolve
        """
        CHKERR( NEPSolve(self.nep) )

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
        getConvergedReason, setTolerances, slepc.NEPGetIterationNumber
        """
        cdef PetscInt ival = 0
        CHKERR( NEPGetIterationNumber(self.nep, &ival) )
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
        setTolerances, solve, slepc.NEPGetConvergedReason
        """
        cdef SlepcNEPConvergedReason val = NEP_CONVERGED_ITERATING
        CHKERR( NEPGetConvergedReason(self.nep, &val) )
        return val

    def getConverged(self) -> int:
        """
        Get the number of converged eigenpairs.

        Not collective.

        Returns
        -------
        nconv: int
            Number of converged eigenpairs.

        Notes
        -----
        This function should be called after `solve()` has finished.

        The value ``nconv`` may be different from the number of requested
        solutions ``nev``, but not larger than ``ncv``, see `setDimensions()`.

        See Also
        --------
        setDimensions, solve, getEigenpair, slepc.NEPGetConverged
        """
        cdef PetscInt ival = 0
        CHKERR( NEPGetConverged(self.nep, &ival) )
        return toInt(ival)

    def getEigenpair(self, i: int, Vec Vr = None, Vec Vi = None) -> None:
        """
        Get the i-th solution of the eigenproblem as computed by `solve()`.

        Collective.

        The solution consists of both the eigenvalue and the eigenvector.

        Parameters
        ----------
        i
            Index of the solution to be obtained.
        Vr
            Placeholder for the returned eigenvector (real part).
        Vi
            Placeholder for the returned eigenvector (imaginary part).

        Returns
        -------
        complex
            The computed eigenvalue.

        Notes
        -----
        The index ``i`` should be a value between ``0`` and ``nconv-1`` (see
        `getConverged()`). Eigenpairs are indexed according to the ordering
        criterion established with `setWhichEigenpairs()`.

        The eigenvector is normalized to have unit norm.

        See Also
        --------
        solve, getConverged, setWhichEigenpairs, slepc.NEPGetEigenpair
        """
        cdef PetscScalar sval1 = 0
        cdef PetscScalar sval2 = 0
        cdef PetscVec vecr = Vr.vec if Vr is not None else <PetscVec>NULL
        cdef PetscVec veci = Vi.vec if Vi is not None else <PetscVec>NULL
        CHKERR( NEPGetEigenpair(self.nep, i, &sval1, &sval2, vecr, veci) )
        return toComplex(sval1, sval2)

    def getLeftEigenvector(self, i: int, Vec Wr, Vec Wi=None) -> None:
        """
        Get the i-th left eigenvector as computed by `solve()`.

        Collective.

        Parameters
        ----------
        i
            Index of the solution to be obtained.
        Wr
            Placeholder for the returned eigenvector (real part).
        Wi
            Placeholder for the returned eigenvector (imaginary part).

        Notes
        -----
        The index ``i`` should be a value between ``0`` and
        ``nconv-1`` (see `getConverged()`). Eigensolutions are indexed
        according to the ordering criterion established with
        `setWhichEigenpairs()`.

        Left eigenvectors are available only if the ``twosided`` flag was
        set with `setTwoSided()`.

        See Also
        --------
        getEigenpair, getConverged, setTwoSided, slepc.NEPGetLeftEigenvector
        """
        cdef PetscVec vecr = Wr.vec if Wr is not None else <PetscVec>NULL
        cdef PetscVec veci = Wi.vec if Wi is not None else <PetscVec>NULL
        CHKERR( NEPGetLeftEigenvector(self.nep, i, vecr, veci) )

    def getErrorEstimate(self, i: int) -> float:
        """
        Get the error estimate associated to the i-th computed eigenpair.

        Not collective.

        Parameters
        ----------
        i
            Index of the solution to be considered.

        Returns
        -------
        float
            Error estimate.

        Notes
        -----
        This is the error estimate used internally by the eigensolver.
        The actual error bound can be computed with `computeError()`.

        See Also
        --------
        computeError, slepc.NEPGetErrorEstimate
        """
        cdef PetscReal rval = 0
        CHKERR( NEPGetErrorEstimate(self.nep, i, &rval) )
        return toReal(rval)

    def computeError(self, i: int, etype: ErrorType | None = None) -> float:
        """
        Compute the error associated with the i-th computed eigenpair.

        Collective.

        Compute the error (based on the residual norm) associated with the
        i-th computed eigenpair.

        Parameters
        ----------
        i
            Index of the solution to be considered.
        etype
            The error type to compute.

        Returns
        -------
        float
            The error bound, computed in various ways from the residual norm
            :math:`\|T(\lambda)x\|_2` where :math:`\lambda` is the eigenvalue
            and :math:`x` is the eigenvector.

        Notes
        -----
        The index ``i`` should be a value between ``0`` and ``nconv-1``
        (see `getConverged()`).

        If the computation of left eigenvectors was enabled with `setTwoSided()`,
        then the error will be computed using the maximum of the value above and
        the left residual norm  :math:`\|y^*T(\lambda)\|_2`, where :math:`y`
        is the approximate left eigenvector.

        See Also
        --------
        getErrorEstimate, setTwoSided, slepc.NEPComputeError
        """
        cdef SlepcNEPErrorType et = NEP_ERROR_RELATIVE
        cdef PetscReal rval = 0
        if etype is not None: et = etype
        CHKERR( NEPComputeError(self.nep, i, et, &rval) )
        return toReal(rval)

    def errorView(self, etype: ErrorType | None = None, viewer: petsc4py.PETSc.Viewer | None = None) -> None:
        """
        Display the errors associated with the computed solution.

        Collective.

        Display the errors and the eigenvalues.

        Parameters
        ----------
        etype
            The error type to compute.
        viewer
            Visualization context; if not provided, the standard
            output is used.

        Notes
        -----
        By default, this function checks the error of all eigenpairs and prints
        the eigenvalues if all of them are below the requested tolerance.
        If the viewer has format ``ASCII_INFO_DETAIL`` then a table with
        eigenvalues and corresponding errors is printed.

        See Also
        --------
        solve, valuesView, vectorsView, slepc.NEPErrorView
        """
        cdef SlepcNEPErrorType et = NEP_ERROR_RELATIVE
        if etype is not None: et = etype
        cdef PetscViewer vwr = def_Viewer(viewer)
        CHKERR( NEPErrorView(self.nep, et, vwr) )

    def valuesView(self, viewer: Viewer | None = None) -> None:
        """
        Display the computed eigenvalues in a viewer.

        Collective.

        Parameters
        ----------
        viewer
            Visualization context; if not provided, the standard
            output is used.

        See Also
        --------
        solve, vectorsView, errorView, slepc.NEPValuesView
        """
        cdef PetscViewer vwr = def_Viewer(viewer)
        CHKERR( NEPValuesView(self.nep, vwr) )

    def vectorsView(self, viewer: Viewer | None = None) -> None:
        """
        Output computed eigenvectors to a viewer.

        Collective.

        Parameters
        ----------
        viewer
            Visualization context; if not provided, the standard
            output is used.

        See Also
        --------
        solve, valuesView, errorView, slepc.NEPVectorsView
        """
        cdef PetscViewer vwr = def_Viewer(viewer)
        CHKERR( NEPVectorsView(self.nep, vwr) )

    #

    def setFunction(
        self,
        function: NEPFunction,
        Mat F: petsc4py.PETSc.Mat | None = None,
        Mat P: petsc4py.PETSc.Mat | None = None,
        args: tuple[Any, ...] | None = None,
        kargs: dict[str, Any] | None = None,
    ) -> None:
        """
        Set the function to compute the nonlinear Function :math:`T(\lambda)`.

        Collective.

        Set the function to compute the nonlinear Function :math:`T(\lambda)`
        as well as the location to store the matrix.

        Parameters
        ----------
        function
            Function evaluation routine.
        F
            Function matrix.
        P
            Preconditioner matrix (usually the same as ``F``).

        See Also
        --------
        setJacobian, getFunction, slepc.NEPSetFunction
        """
        cdef PetscMat Fmat = F.mat if F is not None else <PetscMat>NULL
        cdef PetscMat Pmat = P.mat if P is not None else Fmat
        if function is not None:
            if args is None: args = ()
            if kargs is None: kargs = {}
            context = (function, args, kargs)
            self.set_attr('__function__', context)
            CHKERR( NEPSetFunction(self.nep, Fmat, Pmat, NEP_Function, <void*>context) )
        else:
            CHKERR( NEPSetFunction(self.nep, Fmat, Pmat, NULL, NULL) )

    def getFunction(self) -> tuple[petsc4py.PETSc.Mat, petsc4py.PETSc.Mat, NEPFunction]:
        """
        Get the function to compute the nonlinear Function :math:`T(\lambda)`.

        Collective.

        Get the function to compute the nonlinear Function :math:`T(\lambda)`
        and the matrix.

        Returns
        -------
        F: petsc4py.PETSc.Mat
            Function matrix.
        P: petsc4py.PETSc.Mat
            Preconditioner matrix (usually the same as the F).
        function: NEPFunction
            Function evaluation routine.

        See Also
        --------
        setFunction, slepc.NEPGetFunction
        """
        cdef Mat F = Mat()
        cdef Mat P = Mat()
        CHKERR( NEPGetFunction(self.nep, &F.mat, &P.mat, NULL, NULL) )
        CHKERR( PetscINCREF(F.obj) )
        CHKERR( PetscINCREF(P.obj) )
        cdef object function = self.get_attr('__function__')
        return (F, P, function)

    def setJacobian(
        self,
        jacobian: NEPJacobian,
        Mat J: petsc4py.PETSc.Mat | None = None,
        args: tuple[Any, ...] | None = None,
        kargs: dict[str, Any] | None = None,
    ) -> None:
        """
        Set the function to compute the Jacobian :math:`T'(\lambda)`.

        Collective.

        Set the function to compute the Jacobian :math:`T'(\lambda)` as well as
        the location to store the matrix.

        Parameters
        ----------
        jacobian
            Jacobian evaluation routine.
        J
            Jacobian matrix.

        See Also
        --------
        setFunction, getJacobian, slepc.NEPSetJacobian
        """
        cdef PetscMat Jmat = J.mat if J is not None else <PetscMat>NULL
        if jacobian is not None:
            if args is None: args = ()
            if kargs is None: kargs = {}
            context = (jacobian, args, kargs)
            self.set_attr('__jacobian__', context)
            CHKERR( NEPSetJacobian(self.nep, Jmat, NEP_Jacobian, <void*>context) )
        else:
            CHKERR( NEPSetJacobian(self.nep, Jmat, NULL, NULL) )

    def getJacobian(self) -> tuple[petsc4py.PETSc.Mat, NEPJacobian]:
        """
        Get the function to compute the Jacobian :math:`T'(\lambda)` and J.

        Collective.

        Get the function to compute the Jacobian :math:`T'(\lambda)` and the
        matrix.

        Returns
        -------
        J: petsc4py.PETSc.Mat
            Jacobian matrix.
        jacobian: NEPJacobian
            Jacobian evaluation routine.

        See Also
        --------
        setJacobian, slepc.NEPGetJacobian
        """
        cdef Mat J = Mat()
        CHKERR( NEPGetJacobian(self.nep, &J.mat, NULL, NULL) )
        CHKERR( PetscINCREF(J.obj) )
        cdef object jacobian = self.get_attr('__jacobian__')
        return (J, jacobian)

    def setSplitOperator(
        self,
        A: petsc4py.PETSc.Mat | list[petsc4py.PETSc.Mat],
        f: FN | list[FN],
        structure: petsc4py.PETSc.Mat.Structure | None = None,
    ) -> None:
        """
        Set the operator of the nonlinear eigenvalue problem in split form.

        Collective.

        Parameters
        ----------
        A
            Coefficient matrices of the split form.
        f
            Scalar functions of the split form.
        structure
            Structure flag for matrices.

        Notes
        -----
        The nonlinear operator is written as
        :math:`T(\lambda) = \sum_i A_i f_i(\lambda)`, for :math:`i=1,\dots,n`.
        The derivative :math:`T'(\lambda)` can be obtained using the
        derivatives of :math:`f_i`.

        The ``structure`` flag provides information about :math:`A_i`'s
        nonzero pattern.

        This function must be called before `setUp()`. If it is called
        again after `setUp()` then the `NEP` object is reset.

        See Also
        --------
        getSplitOperator, slepc.NEPSetSplitOperator
        """
        if isinstance(A, Mat): A = [A]
        if isinstance(f, FN):  f = [f]
        cdef PetscMat *As = NULL
        cdef SlepcFN  *Fs = NULL
        cdef Py_ssize_t i = 0, n = len(A)
        cdef PetscMatStructure mstr = matstructure(structure)
        assert n == len(f)
        cdef tmp1 = allocate(<size_t>n*sizeof(PetscMat),<void**>&As)
        cdef tmp2 = allocate(<size_t>n*sizeof(SlepcFN),<void**>&Fs)
        for i in range(n):
            As[i] = (<Mat?>A[i]).mat
            Fs[i] = (<FN?>f[i]).fn
        CHKERR( NEPSetSplitOperator(self.nep, <PetscInt>n, As, Fs, mstr) )

    def getSplitOperator(self) -> tuple[list[petsc4py.PETSc.Mat], list[FN], petsc4py.PETSc.Mat.Structure]:
        """
        Get the operator of the nonlinear eigenvalue problem in split form.

        Collective.

        Returns
        -------
        A: list of petsc4py.PETSc.Mat
            Coefficient matrices of the split form.
        f: list of FN
            Scalar functions of the split form.
        structure: petsc4py.PETSc.Mat.Structure
            Structure flag for matrices.

        See Also
        --------
        setSplitOperator, slepc.NEPGetSplitOperatorInfo, slepc.NEPGetSplitOperatorTerm
        """
        cdef Mat A
        cdef FN  f
        cdef PetscMat mat = NULL
        cdef SlepcFN  fn  = NULL
        cdef PetscInt i=0, n=0
        cdef PetscMatStructure mstr
        CHKERR( NEPGetSplitOperatorInfo(self.nep, &n, &mstr) )
        cdef object matrices = []
        cdef object functions = []
        for i in range(n):
            CHKERR( NEPGetSplitOperatorTerm(self.nep, i, &mat, &fn) )
            A = Mat(); A.mat = mat; CHKERR( PetscINCREF(A.obj) )
            f = FN();  f.fn = fn;   CHKERR( PetscINCREF(f.obj) )
            matrices.append(A)
            functions.append(f)
        return (matrices, functions, mstr)

    def setSplitPreconditioner(
        self,
        P: petsc4py.PETSc.Mat | list[petsc4py.PETSc.Mat],
        structure: petsc4py.PETSc.Mat.Structure | None = None,
    ) -> None:
        """
        Set the operator in split form.

        Collective.

        Set the operator in split form from which to build the preconditioner
        to be used when solving the nonlinear eigenvalue problem in split form.

        Parameters
        ----------
        P
            Coefficient matrices of the split preconditioner.
        structure
            Structure flag for matrices.

        See Also
        --------
        getSplitPreconditioner, slepc.NEPSetSplitPreconditioner
        """
        if isinstance(P, Mat): P = [P]
        cdef PetscMat *Ps = NULL
        cdef Py_ssize_t i = 0, n = len(P)
        cdef PetscMatStructure mstr = matstructure(structure)
        cdef tmp1 = allocate(<size_t>n*sizeof(PetscMat),<void**>&Ps)
        for i in range(n):
            Ps[i] = (<Mat?>P[i]).mat
        CHKERR( NEPSetSplitPreconditioner(self.nep, <PetscInt>n, Ps, mstr) )

    def getSplitPreconditioner(self) -> tuple[list[petsc4py.PETSc.Mat], petsc4py.PETSc.Mat.Structure]:
        """
        Get the operator of the split preconditioner.

        Not collective.

        Returns
        -------
        P: list of petsc4py.PETSc.Mat
            Coefficient matrices of the split preconditioner.
        structure: petsc4py.PETSc.Mat.Structure
            Structure flag for matrices.

        See Also
        --------
        setSplitPreconditioner, slepc.NEPGetSplitPreconditionerTerm
        """
        cdef Mat P
        cdef PetscMat mat = NULL
        cdef PetscInt i=0, n=0
        cdef PetscMatStructure mstr
        CHKERR( NEPGetSplitPreconditionerInfo(self.nep, &n, &mstr) )
        cdef object matrices = []
        for i in range(n):
            CHKERR( NEPGetSplitPreconditionerTerm(self.nep, i, &mat) )
            P = Mat(); P.mat = mat; CHKERR( PetscINCREF(P.obj) )
            matrices.append(P)
        return (matrices, mstr)

    def getTwoSided(self) -> bool:
        """
        Get the flag indicating if a two-sided variant is being used.

        Not collective.

        Get the flag indicating whether a two-sided variant of the algorithm
        is being used or not.

        Returns
        -------
        bool
            Whether the two-sided variant is to be used or not.

        See Also
        --------
        setTwoSided, slepc.NEPGetTwoSided
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( NEPGetTwoSided(self.nep, &tval) )
        return toBool(tval)

    def setTwoSided(self, twosided: bool) -> None:
        """
        Set the solver to use a two-sided variant.

        Logically collective.

        Set the solver to use a two-sided variant so that left eigenvectors
        are also computed.

        Parameters
        ----------
        twosided
            Whether the two-sided variant is to be used or not.

        Notes
        -----
        If the user sets ``twosided`` to ``True`` then the solver uses a
        variant of the algorithm that computes both right and left
        eigenvectors. This is usually much more costly. This option is not
        available in all solvers.

        When using two-sided solvers, the problem matrices must have both
        the ``Mat.mult`` and ``Mat.multTranspose`` operations defined.

        See Also
        --------
        getTwoSided, getLeftEigenvector, slepc.NEPSetTwoSided
        """
        cdef PetscBool tval = asBool(twosided)
        CHKERR( NEPSetTwoSided(self.nep, tval) )

    def applyResolvent(
        self,
        omega: Scalar,
        Vec v,
        Vec r,
        RG rg = None,
    ) -> None:
        """
        Apply the resolvent :math:`T^{-1}(z)` to a given vector.

        Collective.

        Parameters
        ----------
        omega
            Value where the resolvent must be evaluated.
        v
            Input vector.
        r
            Placeholder for the result vector.
        rg
            Region.

        Notes
        -----
        The resolvent :math:`T^{-1}(z)=\sum_i(z-\lambda_i)^{-1}x_iy_i^*`
        is evaluated at :math:`z=\omega` and the matrix-vector product
        :math:`r = T^{-1}(\omega) v` is computed. Vectors :math:`x_i,y_i`
        are right and left eigenvectors, respectively, normalized so that
        :math:`y_i^*T'(\lambda_i)x_i=1`. The sum contains only eigenvectors
        that have been previously computed with `solve()`, and if a region
        ``rg`` is given then only those corresponding to eigenvalues inside
        the region are considered.

        See Also
        --------
        solve, getLeftEigenvector, slepc.NEPApplyResolvent
        """
        cdef PetscScalar sval = asScalar(omega)
        cdef SlepcRG region = rg.rg if rg is not None else <SlepcRG>NULL
        CHKERR( NEPApplyResolvent(self.nep, region, sval, v.vec, r.vec) )

    #

    def setRIILagPreconditioner(self, lag: int) -> None:
        """
        Set when the preconditioner is rebuilt in the nonlinear solve.

        Logically collective.

        Parameters
        ----------
        lag
            0 indicates NEVER rebuild, 1 means rebuild every time the Jacobian is
            computed within the nonlinear iteration, 2 means every second time
            the Jacobian is built, etc.

        See Also
        --------
        getRIILagPreconditioner, slepc.NEPRIISetLagPreconditioner
        """
        cdef PetscInt ival = asInt(lag)
        CHKERR( NEPRIISetLagPreconditioner(self.nep, ival) )

    def getRIILagPreconditioner(self) -> int:
        """
        Get how often the preconditioner is rebuilt.

        Not collective.

        Returns
        -------
        int
            The lag parameter.

        See Also
        --------
        setRIILagPreconditioner, slepc.NEPRIIGetLagPreconditioner
        """
        cdef PetscInt ival = 0
        CHKERR( NEPRIIGetLagPreconditioner(self.nep, &ival) )
        return toInt(ival)

    def setRIIConstCorrectionTol(self, cct: bool) -> None:
        """
        Set a flag to keep the tolerance used in the linear solver constant.

        Logically collective.

        Parameters
        ----------
        cct
             If ``True``, the `petsc4py.PETSc.KSP` relative tolerance is constant.

        Notes
        -----
        By default, an exponentially decreasing tolerance is set in the
        ``KSP`` used within the nonlinear iteration, so that each Newton
        iteration requests better accuracy than the previous one. The
        constant correction tolerance flag stops this behavior.

        See Also
        --------
        getRIIConstCorrectionTol, slepc.NEPRIISetConstCorrectionTol
        """
        cdef PetscBool val = asBool(cct)
        CHKERR( NEPRIISetConstCorrectionTol(self.nep, val) )

    def getRIIConstCorrectionTol(self) -> bool:
        """
        Get the constant tolerance flag.

        Not collective.

        Returns
        -------
        bool
            If ``True``, the `petsc4py.PETSc.KSP` relative tolerance is
            constant.

        See Also
        --------
        setRIIConstCorrectionTol, slepc.NEPRIIGetConstCorrectionTol
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( NEPRIIGetConstCorrectionTol(self.nep, &tval) )
        return toBool(tval)

    def setRIIMaximumIterations(self, its: int) -> None:
        """
        Set the max. number of inner iterations to be used in the RII solver.

        Logically collective.

        These are the Newton iterations related to the computation of the
        nonlinear Rayleigh functional.

        Parameters
        ----------
        its
             Maximum inner iterations.

        See Also
        --------
        getRIIMaximumIterations, slepc.NEPRIISetMaximumIterations
        """
        cdef PetscInt ival = asInt(its)
        CHKERR( NEPRIISetMaximumIterations(self.nep, ival) )

    def getRIIMaximumIterations(self) -> int:
        """
        Get the maximum number of inner iterations of RII.

        Not collective.

        Returns
        -------
        int
            Maximum inner iterations.

        See Also
        --------
        setRIIMaximumIterations, slepc.NEPRIIGetMaximumIterations
        """
        cdef PetscInt ival = 0
        CHKERR( NEPRIIGetMaximumIterations(self.nep, &ival) )
        return toInt(ival)

    def setRIIHermitian(self, herm: bool) -> None:
        """
        Set a flag to use the Hermitian version of the solver.

        Logically collective.

        Set a flag to indicate if the Hermitian version of the scalar
        nonlinear equation must be used by the solver.

        Parameters
        ----------
        herm
            If ``True``, the Hermitian version is used.

        Notes
        -----
        By default, the scalar nonlinear equation
        :math:`x^*T(\sigma)^{-1}T(z)x=0` is solved at each step of the
        nonlinear iteration. When this flag is set the simpler form
        :math:`x^*T(z)x=0` is used, which is supposed to be valid only
        for Hermitian problems.

        See Also
        --------
        getRIIHermitian, slepc.NEPRIISetHermitian
        """
        cdef PetscBool val = asBool(herm)
        CHKERR( NEPRIISetHermitian(self.nep, val) )

    def getRIIHermitian(self) -> bool:
        """
        Get if the Hermitian version must be used by the solver.

        Not collective.

        Returns
        -------
        bool
            If ``True``, the Hermitian version is used.

        See Also
        --------
        setRIIHermitian, slepc.NEPRIIGetHermitian
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( NEPRIIGetHermitian(self.nep, &tval) )
        return toBool(tval)

    def setRIIDeflationThreshold(self, deftol: float) -> None:
        """
        Set the threshold used to switch between deflated and non-deflated.

        Logically collective.

        Set the threshold value used to switch between deflated and
        non-deflated iteration.

        Parameters
        ----------
        deftol
            The threshold value.

        Notes
        -----
        Normally, the solver iterates on the extended problem in order
        to deflate previously converged eigenpairs. If this threshold
        is set to a nonzero value, then once the residual error is below
        this threshold the solver will continue the iteration without
        deflation. The intention is to be able to improve the current
        eigenpair further, despite having previous eigenpairs with
        somewhat bad precision.

        See Also
        --------
        getRIIDeflationThreshold, slepc.NEPRIISetDeflationThreshold
        """
        cdef PetscReal val = asReal(deftol)
        CHKERR( NEPRIISetDeflationThreshold(self.nep, val) )

    def getRIIDeflationThreshold(self) -> float:
        """
        Get the threshold value that controls deflation.

        Not collective.

        Returns
        -------
        float
            The threshold value.

        See Also
        --------
        setRIIDeflationThreshold, slepc.NEPRIIGetDeflationThreshold
        """
        cdef PetscReal rval = 0.0
        CHKERR( NEPRIIGetDeflationThreshold(self.nep, &rval) )
        return toReal(rval)

    def setRIIKSP(self, KSP ksp) -> None:
        """
        Set a linear solver object associated to the nonlinear eigensolver.

        Collective.

        Parameters
        ----------
        ksp
            The linear solver object.

        See Also
        --------
        getRIIKSP, slepc.NEPRIISetKSP
        """
        CHKERR( NEPRIISetKSP(self.nep, ksp.ksp) )

    def getRIIKSP(self) -> KSP:
        """
        Get the linear solver object associated with the nonlinear eigensolver.

        Collective.

        Returns
        -------
        petsc4py.PETSc.KSP
            The linear solver object.

        See Also
        --------
        setRIIKSP, slepc.NEPRIIGetKSP
        """
        cdef KSP ksp = KSP()
        CHKERR( NEPRIIGetKSP(self.nep, &ksp.ksp) )
        CHKERR( PetscINCREF(ksp.obj) )
        return ksp

    #

    def setSLPDeflationThreshold(self, deftol: float) -> None:
        """
        Set the threshold used to switch between deflated and non-deflated.

        Logically collective.

        Parameters
        ----------
        deftol
            The threshold value.

        Notes
        -----
        Normally, the solver iterates on the extended problem in order
        to deflate previously converged eigenpairs. If this threshold
        is set to a nonzero value, then once the residual error is below
        this threshold the solver will continue the iteration without
        deflation. The intention is to be able to improve the current
        eigenpair further, despite having previous eigenpairs with
        somewhat bad precision.

        See Also
        --------
        getSLPDeflationThreshold, slepc.NEPSLPSetDeflationThreshold
        """
        cdef PetscReal val = asReal(deftol)
        CHKERR( NEPSLPSetDeflationThreshold(self.nep, val) )

    def getSLPDeflationThreshold(self) -> float:
        """
        Get the threshold value that controls deflation.

        Not collective.

        Returns
        -------
        float
            The threshold value.

        See Also
        --------
        setSLPDeflationThreshold, slepc.NEPSLPGetDeflationThreshold
        """
        cdef PetscReal rval = 0.0
        CHKERR( NEPSLPGetDeflationThreshold(self.nep, &rval) )
        return toReal(rval)

    def setSLPEPS(self, EPS eps) -> None:
        """
        Set a linear eigensolver object associated to the nonlinear eigensolver.

        Collective.

        Parameters
        ----------
        eps
            The linear eigensolver.

        See Also
        --------
        getSLPEPS, slepc.NEPSLPSetEPS
        """
        CHKERR( NEPSLPSetEPS(self.nep, eps.eps) )

    def getSLPEPS(self) -> EPS:
        """
        Get the linear eigensolver object associated with the nonlinear eigensolver.

        Collective.

        Returns
        -------
        EPS
            The linear eigensolver.

        See Also
        --------
        setSLPEPS, slepc.NEPSLPGetEPS
        """
        cdef EPS eps = EPS()
        CHKERR( NEPSLPGetEPS(self.nep, &eps.eps) )
        CHKERR( PetscINCREF(eps.obj) )
        return eps

    def setSLPEPSLeft(self, EPS eps) -> None:
        """
        Set a linear eigensolver object associated to the nonlinear eigensolver.

        Collective.

        Used to compute left eigenvectors in the two-sided variant of SLP.

        Parameters
        ----------
        eps
            The linear eigensolver.

        See Also
        --------
        setTwoSided, setSLPEPS, getSLPEPSLeft, slepc.NEPSLPSetEPSLeft
        """
        CHKERR( NEPSLPSetEPSLeft(self.nep, eps.eps) )

    def getSLPEPSLeft(self) -> EPS:
        """
        Get the left eigensolver.

        Collective.

        Returns
        -------
        EPS
            The linear eigensolver.

        See Also
        --------
        setSLPEPSLeft, slepc.NEPSLPGetEPSLeft
        """
        cdef EPS eps = EPS()
        CHKERR( NEPSLPGetEPSLeft(self.nep, &eps.eps) )
        CHKERR( PetscINCREF(eps.obj) )
        return eps

    def setSLPKSP(self, KSP ksp) -> None:
        """
        Set a linear solver object associated to the nonlinear eigensolver.

        Collective.

        Parameters
        ----------
        ksp
            The linear solver object.

        See Also
        --------
        getSLPKSP, slepc.NEPSLPSetKSP
        """
        CHKERR( NEPSLPSetKSP(self.nep, ksp.ksp) )

    def getSLPKSP(self) -> KSP:
        """
        Get the linear solver object associated with the nonlinear eigensolver.

        Collective.

        Returns
        -------
        petsc4py.PETSc.KSP
            The linear solver object.

        See Also
        --------
        setSLPKSP, slepc.NEPSLPGetKSP
        """
        cdef KSP ksp = KSP()
        CHKERR( NEPSLPGetKSP(self.nep, &ksp.ksp) )
        CHKERR( PetscINCREF(ksp.obj) )
        return ksp

    #

    def setNArnoldiKSP(self, KSP ksp) -> None:
        """
        Set a linear solver object associated to the nonlinear eigensolver.

        Collective.

        Parameters
        ----------
        ksp
            The linear solver object.

        See Also
        --------
        getNArnoldiKSP, slepc.NEPNArnoldiSetKSP
        """
        CHKERR( NEPNArnoldiSetKSP(self.nep, ksp.ksp) )

    def getNArnoldiKSP(self) -> KSP:
        """
        Get the linear solver object associated with the nonlinear eigensolver.

        Collective.

        Returns
        -------
        petsc4py.PETSc.KSP
            The linear solver object.

        See Also
        --------
        setNArnoldiKSP, slepc.NEPNArnoldiGetKSP
        """
        cdef KSP ksp = KSP()
        CHKERR( NEPNArnoldiGetKSP(self.nep, &ksp.ksp) )
        CHKERR( PetscINCREF(ksp.obj) )
        return ksp

    def setNArnoldiLagPreconditioner(self, lag: int) -> None:
        """
        Set when the preconditioner is rebuilt in the nonlinear solve.

        Logically collective.

        Parameters
        ----------
        lag
            0 indicates NEVER rebuild, 1 means rebuild every time the Jacobian is
            computed within the nonlinear iteration, 2 means every second time
            the Jacobian is built, etc.

        Notes
        -----
        The default is 1. The preconditioner is ALWAYS built in the first
        iteration of a nonlinear solve.

        See Also
        --------
        getNArnoldiLagPreconditioner, slepc.NEPNArnoldiSetLagPreconditioner
        """
        cdef PetscInt ival = asInt(lag)
        CHKERR( NEPNArnoldiSetLagPreconditioner(self.nep, ival) )

    def getNArnoldiLagPreconditioner(self) -> int:
        """
        Get how often the preconditioner is rebuilt.

        Not collective.

        Returns
        -------
        int
            The lag parameter.

        See Also
        --------
        setNArnoldiLagPreconditioner, slepc.NEPNArnoldiGetLagPreconditioner
        """
        cdef PetscInt ival = 0
        CHKERR( NEPNArnoldiGetLagPreconditioner(self.nep, &ival) )
        return toInt(ival)

    #

    def setInterpolPEP(self, PEP pep) -> None:
        """
        Set a polynomial eigensolver object associated to the nonlinear eigensolver.

        Collective.

        Parameters
        ----------
        pep
            The polynomial eigensolver.

        See Also
        --------
        getInterpolPEP, slepc.NEPInterpolSetPEP
        """
        CHKERR( NEPInterpolSetPEP(self.nep, pep.pep) )

    def getInterpolPEP(self) -> PEP:
        """
        Get the associated polynomial eigensolver object.

        Collective.

        Returns
        -------
        PEP
            The polynomial eigensolver.

        See Also
        --------
        setInterpolPEP, slepc.NEPInterpolGetPEP
        """
        cdef PEP pep = PEP()
        CHKERR( NEPInterpolGetPEP(self.nep, &pep.pep) )
        CHKERR( PetscINCREF(pep.obj) )
        return pep

    def setInterpolInterpolation(self, tol: float | None = None, deg: int | None = None) -> None:
        """
        Set the tolerance and maximum degree for the interpolation polynomial.

        Collective.

        Parameters
        ----------
        tol
            The tolerance to stop computing polynomial coefficients.
        deg
            The maximum degree of interpolation.

        See Also
        --------
        getInterpolInterpolation, slepc.NEPInterpolSetInterpolation
        """
        cdef PetscReal rval = PETSC_CURRENT
        cdef PetscInt  ival = PETSC_CURRENT
        if tol is not None: rval = asReal(tol)
        if deg is not None: ival = asInt(deg)
        CHKERR( NEPInterpolSetInterpolation(self.nep, rval, ival) )

    def getInterpolInterpolation(self) -> tuple[float, int]:
        """
        Get the tolerance and maximum degree for the interpolation polynomial.

        Not collective.

        Returns
        -------
        tol: float
            The tolerance to stop computing polynomial coefficients.
        deg: int
            The maximum degree of interpolation.

        See Also
        --------
        setInterpolInterpolation, slepc.NEPInterpolGetInterpolation
        """
        cdef PetscReal rval = 0
        cdef PetscInt  ival = 0
        CHKERR( NEPInterpolGetInterpolation(self.nep, &rval, &ival) )
        return (toReal(rval), toInt(ival))

    #

    def setNLEIGSRestart(self, keep: float) -> None:
        """
        Set the restart parameter for the NLEIGS method.

        Logically collective.

        The proportion of basis vectors that must be kept after restart.

        Parameters
        ----------
        keep
            The number of vectors to be kept at restart.

        Notes
        -----
        Allowed values are in the range [0.1,0.9]. The default is 0.5.

        See Also
        --------
        getNLEIGSRestart, slepc.NEPNLEIGSSetRestart
        """
        cdef PetscReal val = asReal(keep)
        CHKERR( NEPNLEIGSSetRestart(self.nep, val) )

    def getNLEIGSRestart(self) -> float:
        """
        Get the restart parameter used in the NLEIGS method.

        Not collective.

        Returns
        -------
        float
            The number of vectors to be kept at restart.

        See Also
        --------
        setNLEIGSRestart, slepc.NEPNLEIGSGetRestart
        """
        cdef PetscReal val = 0
        CHKERR( NEPNLEIGSGetRestart(self.nep, &val) )
        return toReal(val)

    def setNLEIGSLocking(self, lock: bool) -> None:
        """
        Toggle between locking and non-locking variants of the NLEIGS method.

        Logically collective.

        Parameters
        ----------
        lock
            ``True`` if the locking variant must be selected.

        Notes
        -----
        The default is to lock converged eigenpairs when the method restarts.
        This behavior can be changed so that all directions are kept in the
        working subspace even if already converged to working accuracy (the
        non-locking variant).

        See Also
        --------
        getNLEIGSLocking, slepc.NEPNLEIGSSetLocking
        """
        cdef PetscBool val = asBool(lock)
        CHKERR( NEPNLEIGSSetLocking(self.nep, val) )

    def getNLEIGSLocking(self) -> bool:
        """
        Get the locking flag used in the NLEIGS method.

        Not collective.

        Returns
        -------
        bool
            The locking flag.

        See Also
        --------
        setNLEIGSLocking, slepc.NEPNLEIGSGetLocking
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( NEPNLEIGSGetLocking(self.nep, &tval) )
        return toBool(tval)

    def setNLEIGSInterpolation(self, tol: float | None = None, deg: int | None = None) -> None:
        """
        Set the tolerance and maximum degree for the interpolation polynomial.

        Collective.

        Set the tolerance and maximum degree when building the interpolation
        via divided differences.

        Parameters
        ----------
        tol
            The tolerance to stop computing divided differences.
        deg
            The maximum degree of interpolation.

        See Also
        --------
        getNLEIGSInterpolation, slepc.NEPNLEIGSSetInterpolation
        """
        cdef PetscReal rval = PETSC_CURRENT
        cdef PetscInt  ival = PETSC_CURRENT
        if tol is not None: rval = asReal(tol)
        if deg is not None: ival = asInt(deg)
        CHKERR( NEPNLEIGSSetInterpolation(self.nep, rval, ival) )

    def getNLEIGSInterpolation(self) -> tuple[float, int]:
        """
        Get the tolerance and maximum degree for the interpolation polynomial.

        Not collective.

        Get the tolerance and maximum degree when building the interpolation
        via divided differences.

        Returns
        -------
        tol: float
            The tolerance to stop computing divided differences.
        deg: int
            The maximum degree of interpolation.

        See Also
        --------
        setNLEIGSInterpolation, slepc.NEPNLEIGSGetInterpolation
        """
        cdef PetscReal rval = 0
        cdef PetscInt  ival = 0
        CHKERR( NEPNLEIGSGetInterpolation(self.nep, &rval, &ival) )
        return (toReal(rval), toInt(ival))

    def setNLEIGSFullBasis(self, fullbasis: bool = True) -> None:
        """
        Set TOAR-basis (default) or full-basis variants of the NLEIGS method.

        Logically collective.

        Toggle between TOAR-basis (default) and full-basis variants of the
        NLEIGS method.

        Parameters
        ----------
        fullbasis
            ``True`` if the full-basis variant must be selected.

        Notes
        -----
        The default is to use a compact representation of the Krylov basis,
        that is, :math:`V = (I \otimes U) S`, with a `BV` of type `TENSOR`.
        This behavior can be changed so that the full basis :math:`V` is
        explicitly stored and operated with. This variant is more expensive
        in terms of memory and computation, but is necessary in some cases,
        particularly for two-sided computations, see `setTwoSided()`.

        In the full-basis variant, the NLEIGS solver uses an `EPS` object to
        explicitly solve the linearized eigenproblem, see `getNLEIGSEPS()`.

        See Also
        --------
        setTwoSided, getNLEIGSFullBasis, getNLEIGSEPS, slepc.NEPNLEIGSSetFullBasis
        """
        cdef PetscBool val = asBool(fullbasis)
        CHKERR( NEPNLEIGSSetFullBasis(self.nep, val) )

    def getNLEIGSFullBasis(self) -> bool:
        """
        Get the flag that indicates if NLEIGS is using the full-basis variant.

        Not collective.

        Returns
        -------
        bool
            ``True`` if the full-basis variant is selected.

        See Also
        --------
        setNLEIGSFullBasis, slepc.NEPNLEIGSGetFullBasis
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( NEPNLEIGSGetFullBasis(self.nep, &tval) )
        return toBool(tval)

    def setNLEIGSEPS(self, EPS eps) -> None:
        """
        Set a linear eigensolver object associated to the nonlinear eigensolver.

        Collective.

        Parameters
        ----------
        eps
            The linear eigensolver.

        See Also
        --------
        getNLEIGSEPS, slepc.NEPNLEIGSSetEPS
        """
        CHKERR( NEPNLEIGSSetEPS(self.nep, eps.eps) )

    def getNLEIGSEPS(self) -> EPS:
        """
        Get the linear eigensolver object associated with the nonlinear eigensolver.

        Collective.

        Returns
        -------
        EPS
            The linear eigensolver.

        See Also
        --------
        setNLEIGSEPS, slepc.NEPNLEIGSGetEPS
        """
        cdef EPS eps = EPS()
        CHKERR( NEPNLEIGSGetEPS(self.nep, &eps.eps) )
        CHKERR( PetscINCREF(eps.obj) )
        return eps

    def setNLEIGSRKShifts(self, shifts: Sequence[Scalar]) -> None:
        """
        Set a list of shifts to be used in the Rational Krylov method.

        Collective.

        Parameters
        ----------
        shifts
            Values specifying the shifts.

        Notes
        -----
        If only one shift is provided, the built subspace is equivalent
        to shift-and-invert Krylov-Schur (provided that the absolute
        convergence criterion is used). Otherwise, the rational Krylov
        variant is run.

        See Also
        --------
        getNLEIGSRKShifts, getNLEIGSKSPs, slepc.NEPNLEIGSSetRKShifts
        """
        cdef PetscInt na = 0
        cdef PetscScalar *a = NULL
        cdef object tmp1 = iarray_s(shifts, &na, &a)
        CHKERR( NEPNLEIGSSetRKShifts(self.nep, na, a) )

    def getNLEIGSRKShifts(self) -> ArrayScalar:
        """
        Get the list of shifts used in the Rational Krylov method.

        Not collective.

        Returns
        -------
        ArrayScalar
            The shift values.

        See Also
        --------
        setNLEIGSRKShifts, slepc.NEPNLEIGSGetRKShifts
        """
        cdef PetscInt np = 0
        cdef PetscScalar *coeff = NULL
        CHKERR( NEPNLEIGSGetRKShifts(self.nep, &np, &coeff) )
        cdef object ocoeff = None
        try:
            ocoeff = array_s(np, coeff)
        finally:
            CHKERR( PetscFree(coeff) )
        return ocoeff

    def getNLEIGSKSPs(self) -> list[KSP]:
        """
        Get the list of linear solver objects associated with the NLEIGS solver.

        Collective.

        Returns
        -------
        list of `petsc4py.PETSc.KSP`
            The linear solver objects.

        Notes
        -----
        The number of `petsc4py.PETSc.KSP` solvers is equal to the number of
        shifts provided by the user, or 1 if the user did not provide shifts.

        See Also
        --------
        setNLEIGSRKShifts, slepc.NEPNLEIGSGetKSPs
        """
        cdef PetscInt i = 0, n = 0
        cdef PetscKSP *p = NULL
        CHKERR( NEPNLEIGSGetKSPs(self.nep, &n, &p) )
        return [ref_KSP(p[i]) for i from 0 <= i <n]

    #

    def setCISSExtraction(self, extraction: CISSExtraction) -> None:
        """
        Set the extraction technique used in the CISS solver.

        Logically collective.

        Parameters
        ----------
        extraction
            The extraction technique.

        See Also
        --------
        getCISSExtraction, slepc.NEPCISSSetExtraction
        """
        cdef SlepcNEPCISSExtraction val = extraction
        CHKERR( NEPCISSSetExtraction(self.nep, val) )

    def getCISSExtraction(self) -> CISSExtraction:
        """
        Get the extraction technique used in the CISS solver.

        Not collective.

        Returns
        -------
        CISSExtraction
            The extraction technique.

        See Also
        --------
        setCISSExtraction, slepc.NEPCISSGetExtraction
        """
        cdef SlepcNEPCISSExtraction val = NEP_CISS_EXTRACTION_RITZ
        CHKERR( NEPCISSGetExtraction(self.nep, &val) )
        return val

    def setCISSSizes(
        self,
        ip: int | None = None,
        bs: int | None = None,
        ms: int | None = None,
        npart: int | None = None,
        bsmax: int | None = None,
        realmats: bool = False,
    ) -> None:
        """
        Set the values of various size parameters in the CISS solver.

        Logically collective.

        Parameters
        ----------
        ip
            Number of integration points.
        bs
            Block size.
        ms
            Moment size.
        npart
            Number of partitions when splitting the communicator.
        bsmax
            Maximum block size.
        realmats
            ``True`` if A and B are real.

        Notes
        -----
        The default number of partitions is 1. This means the internal
        `petsc4py.PETSc.KSP` object is shared among all processes of the `NEP`
        communicator. Otherwise, the communicator is split into ``npart``
        communicators, so that ``npart`` `petsc4py.PETSc.KSP` solves proceed
        simultaneously.

        See Also
        --------
        getCISSSizes, setCISSThreshold, setCISSRefinement, slepc.NEPCISSSetSizes
        """
        cdef PetscInt  ival1 = PETSC_CURRENT
        cdef PetscInt  ival2 = PETSC_CURRENT
        cdef PetscInt  ival3 = PETSC_CURRENT
        cdef PetscInt  ival4 = PETSC_CURRENT
        cdef PetscInt  ival5 = PETSC_CURRENT
        cdef PetscBool bval  = asBool(realmats)
        if ip    is not None: ival1 = asInt(ip)
        if bs    is not None: ival2 = asInt(bs)
        if ms    is not None: ival3 = asInt(ms)
        if npart is not None: ival4 = asInt(npart)
        if bsmax is not None: ival5 = asInt(bsmax)
        CHKERR( NEPCISSSetSizes(self.nep, ival1, ival2, ival3, ival4, ival5, bval) )

    def getCISSSizes(self) -> tuple[int, int, int, int, int, bool]:
        """
        Get the values of various size parameters in the CISS solver.

        Not collective.

        Returns
        -------
        ip: int
            Number of integration points.
        bs: int
            Block size.
        ms: int
            Moment size.
        npart: int
            Number of partitions when splitting the communicator.
        bsmax: int
            Maximum block size.
        realmats: bool
            ``True`` if A and B are real.

        See Also
        --------
        setCISSSizes, slepc.NEPCISSGetSizes
        """
        cdef PetscInt  ival1 = 0
        cdef PetscInt  ival2 = 0
        cdef PetscInt  ival3 = 0
        cdef PetscInt  ival4 = 0
        cdef PetscInt  ival5 = 0
        cdef PetscBool bval  = PETSC_FALSE
        CHKERR( NEPCISSGetSizes(self.nep, &ival1, &ival2, &ival3, &ival4, &ival5, &bval) )
        return (toInt(ival1), toInt(ival2), toInt(ival3), toInt(ival4), toInt(ival5), toBool(bval))

    def setCISSThreshold(self, delta: float | None = None, spur: float | None = None) -> None:
        """
        Set the values of various threshold parameters in the CISS solver.

        Logically collective.

        Parameters
        ----------
        delta
            Threshold for numerical rank.
        spur
            Spurious threshold (to discard spurious eigenpairs).

        See Also
        --------
        getCISSThreshold, slepc.NEPCISSSetThreshold
        """
        cdef PetscReal rval1 = PETSC_CURRENT
        cdef PetscReal rval2 = PETSC_CURRENT
        if delta is not None: rval1 = asReal(delta)
        if spur  is not None: rval2 = asReal(spur)
        CHKERR( NEPCISSSetThreshold(self.nep, rval1, rval2) )

    def getCISSThreshold(self) -> tuple[float, float]:
        """
        Get the values of various threshold parameters in the CISS solver.

        Not collective.

        Returns
        -------
        delta: float
            Threshold for numerical rank.
        spur: float
            Spurious threshold (to discard spurious eigenpairs.

        See Also
        --------
        setCISSThreshold, slepc.NEPCISSGetThreshold
        """
        cdef PetscReal delta = 0
        cdef PetscReal spur  = 0
        CHKERR( NEPCISSGetThreshold(self.nep, &delta, &spur) )
        return (toReal(delta), toReal(spur))

    def setCISSRefinement(self, inner: int | None = None, blsize: int | None = None) -> None:
        """
        Set the values of various refinement parameters in the CISS solver.

        Logically collective.

        Parameters
        ----------
        inner
            Number of iterative refinement iterations (inner loop).
        blsize
            Number of iterative refinement iterations (blocksize loop).

        See Also
        --------
        getCISSRefinement, slepc.NEPCISSSetRefinement
        """
        cdef PetscInt ival1 = PETSC_CURRENT
        cdef PetscInt ival2 = PETSC_CURRENT
        if inner  is not None: ival1 = asInt(inner)
        if blsize is not None: ival2 = asInt(blsize)
        CHKERR( NEPCISSSetRefinement(self.nep, ival1, ival2) )

    def getCISSRefinement(self) -> tuple[int, int]:
        """
        Get the values of various refinement parameters in the CISS solver.

        Not collective.

        Returns
        -------
        inner: int
            Number of iterative refinement iterations (inner loop).
        blsize: int
            Number of iterative refinement iterations (blocksize loop).

        See Also
        --------
        setCISSRefinement, slepc.NEPCISSGetRefinement
        """
        cdef PetscInt ival1 = 0
        cdef PetscInt ival2 = 0
        CHKERR( NEPCISSGetRefinement(self.nep, &ival1, &ival2) )
        return (toInt(ival1), toInt(ival2))

    def getCISSKSPs(self) -> list[KSP]:
        """
        Get the list of linear solver objects associated with the CISS solver.

        Collective.

        Returns
        -------
        list of `petsc4py.PETSc.KSP`
            The linear solver objects.

        Notes
        -----
        The number of `petsc4py.PETSc.KSP` solvers is equal to the number of
        integration points divided by the number of partitions. This value is
        halved in the case of real matrices with a region centered at the real
        axis.

        See Also
        --------
        setCISSSizes, slepc.NEPCISSGetKSPs
        """
        cdef PetscInt i = 0, n = 0
        cdef PetscKSP *p = NULL
        CHKERR( NEPCISSGetKSPs(self.nep, &n, &p) )
        return [ref_KSP(p[i]) for i from 0 <= i <n]

    property problem_type:
        """The problem type from the NEP object."""
        def __get__(self) -> NEPProblemType:
            return self.getProblemType()
        def __set__(self, value):
            self.setProblemType(value)

    property which:
        """The portion of the spectrum to be sought."""
        def __get__(self) -> NEPWhich:
            return self.getWhichEigenpairs()
        def __set__(self, value):
            self.setWhichEigenpairs(value)

    property target:
        """The value of the target."""
        def __get__(self) -> float:
            return self.getTarget()
        def __set__(self, value):
            self.setTarget(value)

    property tol:
        """The tolerance used by the NEP convergence tests."""
        def __get__(self) -> float:
            return self.getTolerances()[0]
        def __set__(self, value):
            self.setTolerances(tol=value)

    property max_it:
        """The maximum iteration count used by the NEP convergence tests."""
        def __get__(self) -> int:
            return self.getTolerances()[1]
        def __set__(self, value):
            self.setTolerances(max_it=value)

    property track_all:
        """Compute the residual of all approximate eigenpairs."""
        def __get__(self) -> bool:
            return self.getTrackAll()
        def __set__(self, value):
            self.setTrackAll(value)

    property bv:
        """The basis vectors (`BV`) object associated."""
        def __get__(self) -> BV:
            return self.getBV()
        def __set__(self, value):
            self.setBV(value)

    property rg:
        """The region (`RG`) object associated."""
        def __get__(self) -> RG:
            return self.getRG()
        def __set__(self, value):
            self.setRG(value)

    property ds:
        """The direct solver (`DS`) object associated."""
        def __get__(self) -> DS:
            return self.getDS()
        def __set__(self, value):
            self.setDS(value)

# -----------------------------------------------------------------------------

del NEPType
del NEPProblemType
del NEPErrorType
del NEPWhich
del NEPConvergedReason
del NEPRefine
del NEPRefineScheme
del NEPConv
del NEPStop
del NEPCISSExtraction

# -----------------------------------------------------------------------------
