# -----------------------------------------------------------------------------

class EPSType(object):
    """
    EPS type.

    Native sparse eigensolvers.

    - `POWER`:        Power Iteration, Inverse Iteration, RQI.
    - `SUBSPACE`:     Subspace Iteration.
    - `ARNOLDI`:      Arnoldi.
    - `LANCZOS`:      Lanczos.
    - `KRYLOVSCHUR`:  Krylov-Schur (default).
    - `GD`:           Generalized Davidson.
    - `JD`:           Jacobi-Davidson.
    - `RQCG`:         Rayleigh Quotient Conjugate Gradient.
    - `LOBPCG`:       Locally Optimal Block Preconditioned Conjugate Gradient.
    - `CISS`:         Contour Integral Spectrum Slicing.
    - `LYAPII`:       Lyapunov inverse iteration.
    - `LAPACK`:       Wrappers to dense eigensolvers in Lapack.

    Wrappers to external eigensolvers
    (should be enabled during installation of SLEPc)

    - `ARPACK`:
    - `BLOPEX`:
    - `PRIMME`:
    - `FEAST`:
    - `SCALAPACK`:
    - `ELPA`:
    - `ELEMENTAL`:
    - `EVSL`:
    - `CHASE`:
    """
    # provided implementations
    POWER        = S_(EPSPOWER)
    SUBSPACE     = S_(EPSSUBSPACE)
    ARNOLDI      = S_(EPSARNOLDI)
    LANCZOS      = S_(EPSLANCZOS)
    KRYLOVSCHUR  = S_(EPSKRYLOVSCHUR)
    GD           = S_(EPSGD)
    JD           = S_(EPSJD)
    RQCG         = S_(EPSRQCG)
    LOBPCG       = S_(EPSLOBPCG)
    CISS         = S_(EPSCISS)
    LYAPII       = S_(EPSLYAPII)
    LAPACK       = S_(EPSLAPACK)
    # with external libraries
    ARPACK       = S_(EPSARPACK)
    BLOPEX       = S_(EPSBLOPEX)
    PRIMME       = S_(EPSPRIMME)
    FEAST        = S_(EPSFEAST)
    SCALAPACK    = S_(EPSSCALAPACK)
    ELPA         = S_(EPSELPA)
    ELEMENTAL    = S_(EPSELEMENTAL)
    EVSL         = S_(EPSEVSL)
    CHASE        = S_(EPSCHASE)

class EPSProblemType(object):
    """
    EPS problem type.

    - `HEP`:    Hermitian eigenproblem.
    - `NHEP`:   Non-Hermitian eigenproblem.
    - `GHEP`:   Generalized Hermitian eigenproblem.
    - `GNHEP`:  Generalized Non-Hermitian eigenproblem.
    - `PGNHEP`: Generalized Non-Hermitian eigenproblem
                with positive definite :math:`B`.
    - `GHIEP`:  Generalized Hermitian-indefinite eigenproblem.
    - `BSE`:    Structured Bethe-Salpeter eigenproblem.
    - `HAMILT`: Hamiltonian eigenproblem.
    """
    HEP    = EPS_HEP
    NHEP   = EPS_NHEP
    GHEP   = EPS_GHEP
    GNHEP  = EPS_GNHEP
    PGNHEP = EPS_PGNHEP
    GHIEP  = EPS_GHIEP
    BSE    = EPS_BSE
    HAMILT = EPS_HAMILT

class EPSExtraction(object):
    """
    EPS extraction technique.

    - `RITZ`:              Standard Rayleigh-Ritz extraction.
    - `HARMONIC`:          Harmonic extraction.
    - `HARMONIC_RELATIVE`: Harmonic extraction relative to the eigenvalue.
    - `HARMONIC_RIGHT`:    Harmonic extraction for rightmost eigenvalues.
    - `HARMONIC_LARGEST`:  Harmonic extraction for largest magnitude (without
                           target).
    - `REFINED`:           Refined extraction.
    - `REFINED_HARMONIC`:  Refined harmonic extraction.
    """
    RITZ              = EPS_RITZ
    HARMONIC          = EPS_HARMONIC
    HARMONIC_RELATIVE = EPS_HARMONIC_RELATIVE
    HARMONIC_RIGHT    = EPS_HARMONIC_RIGHT
    HARMONIC_LARGEST  = EPS_HARMONIC_LARGEST
    REFINED           = EPS_REFINED
    REFINED_HARMONIC  = EPS_REFINED_HARMONIC

class EPSBalance(object):
    """
    EPS type of balancing used for non-Hermitian problems.

    - `NONE`:     None.
    - `ONESIDE`:  One-sided balancing.
    - `TWOSIDE`:  Two-sided balancing.
    - `USER`:     User-provided balancing matrices.
    """
    NONE    = EPS_BALANCE_NONE
    ONESIDE = EPS_BALANCE_ONESIDE
    TWOSIDE = EPS_BALANCE_TWOSIDE
    USER    = EPS_BALANCE_USER

class EPSErrorType(object):
    """
    EPS error type to assess accuracy of computed solutions.

    - `ABSOLUTE`: Absolute error.
    - `RELATIVE`: Relative error.
    - `BACKWARD`: Backward error.
    """
    ABSOLUTE = EPS_ERROR_ABSOLUTE
    RELATIVE = EPS_ERROR_RELATIVE
    BACKWARD = EPS_ERROR_BACKWARD

class EPSWhich(object):
    """
    EPS desired part of spectrum.

    - `LARGEST_MAGNITUDE`:  Largest magnitude (default).
    - `SMALLEST_MAGNITUDE`: Smallest magnitude.
    - `LARGEST_REAL`:       Largest real parts.
    - `SMALLEST_REAL`:      Smallest real parts.
    - `LARGEST_IMAGINARY`:  Largest imaginary parts in magnitude.
    - `SMALLEST_IMAGINARY`: Smallest imaginary parts in magnitude.
    - `TARGET_MAGNITUDE`:   Closest to target (in magnitude).
    - `TARGET_REAL`:        Real part closest to target.
    - `TARGET_IMAGINARY`:   Imaginary part closest to target.
    - `ALL`:                All eigenvalues in an interval.
    - `USER`:               User defined selection.
    """
    LARGEST_MAGNITUDE  = EPS_LARGEST_MAGNITUDE
    SMALLEST_MAGNITUDE = EPS_SMALLEST_MAGNITUDE
    LARGEST_REAL       = EPS_LARGEST_REAL
    SMALLEST_REAL      = EPS_SMALLEST_REAL
    LARGEST_IMAGINARY  = EPS_LARGEST_IMAGINARY
    SMALLEST_IMAGINARY = EPS_SMALLEST_IMAGINARY
    TARGET_MAGNITUDE   = EPS_TARGET_MAGNITUDE
    TARGET_REAL        = EPS_TARGET_REAL
    TARGET_IMAGINARY   = EPS_TARGET_IMAGINARY
    ALL                = EPS_ALL
    USER               = EPS_WHICH_USER

class EPSConv(object):
    """
    EPS convergence test.

    - `ABS`:  Absolute convergence test.
    - `REL`:  Convergence test relative to the eigenvalue.
    - `NORM`: Convergence test relative to the matrix norms.
    - `USER`: User-defined convergence test.
    """
    ABS  = EPS_CONV_ABS
    REL  = EPS_CONV_REL
    NORM = EPS_CONV_NORM
    USER = EPS_CONV_USER

class EPSStop(object):
    """
    EPS stopping test.

    - `BASIC`:     Default stopping test.
    - `USER`:      User-defined stopping test.
    - `THRESHOLD`: Threshold stopping test.
    """
    BASIC     = EPS_STOP_BASIC
    USER      = EPS_STOP_USER
    THRESHOLD = EPS_STOP_THRESHOLD

class EPSConvergedReason(object):
    """
    EPS convergence reasons.

    - `CONVERGED_TOL`:          All eigenpairs converged to requested tolerance.
    - `CONVERGED_USER`:         User-defined convergence criterion satisfied.
    - `DIVERGED_ITS`:           Maximum number of iterations exceeded.
    - `DIVERGED_BREAKDOWN`:     Solver failed due to breakdown.
    - `DIVERGED_SYMMETRY_LOST`: Lanczos-type method could not preserve symmetry.
    - `CONVERGED_ITERATING`:    Iteration not finished yet.
    """
    CONVERGED_TOL          = EPS_CONVERGED_TOL
    CONVERGED_USER         = EPS_CONVERGED_USER
    DIVERGED_ITS           = EPS_DIVERGED_ITS
    DIVERGED_BREAKDOWN     = EPS_DIVERGED_BREAKDOWN
    DIVERGED_SYMMETRY_LOST = EPS_DIVERGED_SYMMETRY_LOST
    CONVERGED_ITERATING    = EPS_CONVERGED_ITERATING
    ITERATING              = EPS_CONVERGED_ITERATING

class EPSPowerShiftType(object):
    """
    EPS Power shift type.

    - `CONSTANT`:  Constant shift.
    - `RAYLEIGH`:  Rayleigh quotient.
    - `WILKINSON`: Wilkinson shift.
    """
    CONSTANT  = EPS_POWER_SHIFT_CONSTANT
    RAYLEIGH  = EPS_POWER_SHIFT_RAYLEIGH
    WILKINSON = EPS_POWER_SHIFT_WILKINSON

class EPSKrylovSchurBSEType(object):
    """
    EPS Krylov-Schur method for BSE problems.

    - `SHAO`:         Lanczos recurrence for H square.
    - `GRUNING`:      Lanczos recurrence for H.
    - `PROJECTEDBSE`: Lanczos where the projected problem has BSE structure.
    """
    SHAO         = EPS_KRYLOVSCHUR_BSE_SHAO
    GRUNING      = EPS_KRYLOVSCHUR_BSE_GRUNING
    PROJECTEDBSE = EPS_KRYLOVSCHUR_BSE_PROJECTEDBSE

class EPSLanczosReorthogType(object):
    """
    EPS Lanczos reorthogonalization type.

    - `LOCAL`:     Local reorthogonalization only.
    - `FULL`:      Full reorthogonalization.
    - `SELECTIVE`: Selective reorthogonalization.
    - `PERIODIC`:  Periodic reorthogonalization.
    - `PARTIAL`:   Partial reorthogonalization.
    - `DELAYED`:   Delayed reorthogonalization.
    """
    LOCAL     = EPS_LANCZOS_REORTHOG_LOCAL
    FULL      = EPS_LANCZOS_REORTHOG_FULL
    SELECTIVE = EPS_LANCZOS_REORTHOG_SELECTIVE
    PERIODIC  = EPS_LANCZOS_REORTHOG_PERIODIC
    PARTIAL   = EPS_LANCZOS_REORTHOG_PARTIAL
    DELAYED   = EPS_LANCZOS_REORTHOG_DELAYED

class EPSCISSQuadRule(object):
    """
    EPS CISS quadrature rule.

    - `TRAPEZOIDAL`: Trapezoidal rule.
    - `CHEBYSHEV`:   Chebyshev points.
    """
    TRAPEZOIDAL = EPS_CISS_QUADRULE_TRAPEZOIDAL
    CHEBYSHEV   = EPS_CISS_QUADRULE_CHEBYSHEV

class EPSCISSExtraction(object):
    """
    EPS CISS extraction technique.

    - `RITZ`:   Ritz extraction.
    - `HANKEL`: Extraction via Hankel eigenproblem.
    """
    RITZ   = EPS_CISS_EXTRACTION_RITZ
    HANKEL = EPS_CISS_EXTRACTION_HANKEL

# -----------------------------------------------------------------------------

cdef class EPS(Object):

    """EPS."""

    Type            = EPSType
    ProblemType     = EPSProblemType
    Extraction      = EPSExtraction
    Balance         = EPSBalance
    ErrorType       = EPSErrorType
    Which           = EPSWhich
    Conv            = EPSConv
    Stop            = EPSStop
    ConvergedReason = EPSConvergedReason

    PowerShiftType      = EPSPowerShiftType
    KrylovSchurBSEType  = EPSKrylovSchurBSEType
    LanczosReorthogType = EPSLanczosReorthogType
    CISSQuadRule        = EPSCISSQuadRule
    CISSExtraction      = EPSCISSExtraction

    def __cinit__(self):
        self.obj = <PetscObject*> &self.eps
        self.eps = NULL

    def view(self, Viewer viewer=None) -> None:
        """
        Print the EPS data structure.

        Collective.

        Parameters
        ----------
        viewer
            Visualization context; if not provided, the standard
            output is used.
        """
        cdef PetscViewer vwr = def_Viewer(viewer)
        CHKERR( EPSView(self.eps, vwr) )

    def destroy(self) -> Self:
        """
        Destroy the EPS object.

        Collective.
        """
        CHKERR( EPSDestroy(&self.eps) )
        self.eps = NULL
        return self

    def reset(self) -> None:
        """
        Reset the EPS object.

        Collective.
        """
        CHKERR( EPSReset(self.eps) )

    def create(self, comm: Comm | None = None) -> Self:
        """
        Create the EPS object.

        Collective.

        Parameters
        ----------
        comm
            MPI communicator; if not provided, it defaults to all processes.
        """
        cdef MPI_Comm ccomm = def_Comm(comm, SLEPC_COMM_DEFAULT())
        cdef SlepcEPS neweps = NULL
        CHKERR( EPSCreate(ccomm, &neweps) )
        CHKERR( SlepcCLEAR(self.obj) ); self.eps = neweps
        return self

    def setType(self, eps_type: Type | str) -> None:
        """
        Set the particular solver to be used in the EPS object.

        Logically collective.

        Parameters
        ----------
        eps_type
            The solver to be used.

        Notes
        -----
        See `EPS.Type` for available methods. The default is
        `EPS.Type.KRYLOVSCHUR`.  Normally, it is best to use
        `setFromOptions()` and then set the EPS type from the options
        database rather than by using this routine.  Using the options
        database provides the user with maximum flexibility in
        evaluating the different available methods.
        """
        cdef SlepcEPSType cval = NULL
        eps_type = str2bytes(eps_type, &cval)
        CHKERR( EPSSetType(self.eps, cval) )

    def getType(self) -> str:
        """
        Get the EPS type of this object.

        Not collective.

        Returns
        -------
        str
            The solver currently being used.
        """
        cdef SlepcEPSType eps_type = NULL
        CHKERR( EPSGetType(self.eps, &eps_type) )
        return bytes2str(eps_type)

    def getOptionsPrefix(self) -> str:
        """
        Get the prefix used for searching for all EPS options in the database.

        Not collective.

        Returns
        -------
        str
            The prefix string set for this EPS object.
        """
        cdef const char *prefix = NULL
        CHKERR( EPSGetOptionsPrefix(self.eps, &prefix) )
        return bytes2str(prefix)

    def setOptionsPrefix(self, prefix: str | None = None) -> None:
        """
        Set the prefix used for searching for all EPS options in the database.

        Logically collective.

        Parameters
        ----------
        prefix
            The prefix string to prepend to all EPS option requests.

        Notes
        -----
        A hyphen (-) must NOT be given at the beginning of the prefix
        name.  The first character of all runtime options is
        AUTOMATICALLY the hyphen.

        For example, to distinguish between the runtime options for
        two different EPS contexts, one could call::

            E1.setOptionsPrefix("eig1_")
            E2.setOptionsPrefix("eig2_")
        """
        cdef const char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( EPSSetOptionsPrefix(self.eps, cval) )

    def appendOptionsPrefix(self, prefix: str | None = None) -> None:
        """
        Append to the prefix used for searching for all EPS options in the database.

        Logically collective.

        Parameters
        ----------
        prefix
            The prefix string to prepend to all EPS option requests.
        """
        cdef const char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( EPSAppendOptionsPrefix(self.eps, cval) )

    def setFromOptions(self) -> None:
        """
        Set EPS options from the options database.

        Collective.

        This routine must be called before `setUp()` if the user is to be
        allowed to set the solver type.

        Notes
        -----
        To see all options, run your program with the ``-help``
        option.
        """
        CHKERR( EPSSetFromOptions(self.eps) )

    #

    def getProblemType(self) -> ProblemType:
        """
        Get the problem type from the EPS object.

        Not collective.

        Returns
        -------
        ProblemType
            The problem type that was previously set.
        """
        cdef SlepcEPSProblemType val = EPS_NHEP
        CHKERR( EPSGetProblemType(self.eps, &val) )
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
        Allowed values are: Hermitian (HEP), non-Hermitian (NHEP), generalized
        Hermitian (GHEP), generalized non-Hermitian (GNHEP), and generalized
        non-Hermitian with positive semi-definite B (PGNHEP).

        This function must be used to instruct SLEPc to exploit symmetry. If
        no problem type is specified, by default a non-Hermitian problem is
        assumed (either standard or generalized). If the user knows that the
        problem is Hermitian (i.e. :math:`A=A^H`) or generalized Hermitian
        (i.e. :math:`A=A^H`, :math:`B=B^H`, and :math:`B` positive definite)
        then it is recommended to set the problem type so that eigensolver can
        exploit these properties.
        """
        cdef SlepcEPSProblemType val = problem_type
        CHKERR( EPSSetProblemType(self.eps, val) )

    def isGeneralized(self) -> bool:
        """
        Tell if the EPS object corresponds to a generalized eigenproblem.

        Not collective.

        Returns
        -------
        bool
            True if two matrices were set with `setOperators()`.
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( EPSIsGeneralized(self.eps, &tval) )
        return toBool(tval)

    def isHermitian(self) -> bool:
        """
        Tell if the EPS object corresponds to a Hermitian eigenproblem.

        Not collective.

        Returns
        -------
        bool
            True if the problem type set with `setProblemType()` was Hermitian.
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( EPSIsHermitian(self.eps, &tval) )
        return toBool(tval)

    def isPositive(self) -> bool:
        """
        Eigenproblem requiring a positive (semi-) definite matrix :math:`B`.

        Not collective.

        Tell if the EPS corresponds to an eigenproblem requiring a positive
        (semi-) definite matrix :math:`B`.

        Returns
        -------
        bool
            True if the problem type set with `setProblemType()` was positive.
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( EPSIsPositive(self.eps, &tval) )
        return toBool(tval)

    def isStructured(self) -> bool:
        """
        Tell if the EPS object corresponds to a structured eigenvalue problem.

        Not collective.

        Returns
        -------
        bool
            True if the problem type set with `setProblemType()` was structured.
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( EPSIsStructured(self.eps, &tval) )
        return toBool(tval)

    def getBalance(self) -> tuple[Balance, int, float]:
        """
        Get the balancing type used by the EPS, and the associated parameters.

        Not collective.

        Returns
        -------
        balance: Balance
            The balancing method
        iterations: int
            Number of iterations of the balancing algorithm
        cutoff: float
            Cutoff value
        """
        cdef SlepcEPSBalance val = EPS_BALANCE_ONESIDE
        cdef PetscInt ival = 0
        cdef PetscReal rval = 0
        CHKERR( EPSGetBalance(self.eps, &val, &ival, &rval) )
        return (val, toInt(ival), toReal(rval))

    def setBalance(
        self,
        balance: Balance | None = None,
        iterations: int | None = None,
        cutoff: float | None = None,
    ) -> None:
        """
        Set the balancing technique to be used by the eigensolver.

        Logically collective.

        Set the balancing technique to be used by the eigensolver, and some
        parameters associated to it.

        Parameters
        ----------
        balance
            The balancing method
        iterations
            Number of iterations of the balancing algorithm
        cutoff
            Cutoff value
        """
        cdef SlepcEPSBalance val = EPS_BALANCE_NONE
        cdef PetscInt  ival = PETSC_CURRENT
        cdef PetscReal rval = PETSC_CURRENT
        if balance    is not None: val  = balance
        else: CHKERR( EPSGetBalance(self.eps, &val, NULL, NULL) )
        if iterations is not None: ival = asInt(iterations)
        if cutoff     is not None: rval = asReal(cutoff)
        CHKERR( EPSSetBalance(self.eps, val, ival, rval) )

    def getExtraction(self) -> Extraction:
        """
        Get the extraction type used by the EPS object.

        Not collective.

        Returns
        -------
        Extraction
            The method of extraction.
        """
        cdef SlepcEPSExtraction val = EPS_RITZ
        CHKERR( EPSGetExtraction(self.eps, &val) )
        return val

    def setExtraction(self, extraction: Extraction) -> None:
        """
        Set the extraction type used by the EPS object.

        Logically collective.

        Parameters
        ----------
        extraction
            The extraction method to be used by the solver.

        Notes
        -----
        Not all eigensolvers support all types of extraction. See the
        SLEPc documentation for details.

        By default, a standard Rayleigh-Ritz extraction is used. Other
        extractions may be useful when computing interior eigenvalues.

        Harmonic-type extractions are used in combination with a
        *target*. See `setTarget()`.
        """
        cdef SlepcEPSExtraction val = extraction
        CHKERR( EPSSetExtraction(self.eps, val) )

    def getWhichEigenpairs(self) -> Which:
        """
        Get which portion of the spectrum is to be sought.

        Not collective.

        Returns
        -------
        Which
            The portion of the spectrum to be sought by the solver.
        """
        cdef SlepcEPSWhich val = EPS_LARGEST_MAGNITUDE
        CHKERR( EPSGetWhichEigenpairs(self.eps, &val) )
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
        Not all eigensolvers implemented in EPS account for all the
        possible values. Also, some values make sense only for certain
        types of problems. If SLEPc is compiled for real numbers
        `EPS.Which.LARGEST_IMAGINARY` and
        `EPS.Which.SMALLEST_IMAGINARY` use the absolute value of the
        imaginary part for eigenvalue selection.
        """
        cdef SlepcEPSWhich val = which
        CHKERR( EPSSetWhichEigenpairs(self.eps, val) )

    def getThreshold(self) -> tuple[float, bool]:
        """
        Get the threshold used in the threshold stopping test.

        Not collective.

        Returns
        -------
        thres: float
            The threshold.
        rel: bool
            Whether the threshold is relative or not.
        """
        cdef PetscReal rval = 0
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( EPSGetThreshold(self.eps, &rval, &tval) )
        return (toReal(rval), toBool(tval))

    def setThreshold(self, thres: float, rel: bool = False) -> None:
        """
        Set the threshold used in the threshold stopping test.

        Logically collective.

        Parameters
        ----------
        thres
            The threshold.
        rel
            Whether the threshold is relative or not.

        Notes
        -----
        This function internally sets a special stopping test based on
        the threshold, where eigenvalues are computed in sequence
        until one of the computed eigenvalues is below/above the
        threshold (depending on whether largest or smallest eigenvalues
        are computed).
        """
        cdef PetscReal rval = asReal(thres)
        cdef PetscBool tval = asBool(rel)
        CHKERR( EPSSetThreshold(self.eps, rval, tval) )

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
        """
        cdef PetscScalar sval = 0
        CHKERR( EPSGetTarget(self.eps, &sval) )
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
        """
        cdef PetscScalar sval = asScalar(target)
        CHKERR( EPSSetTarget(self.eps, sval) )

    def getInterval(self) -> tuple[float, float]:
        """
        Get the computational interval for spectrum slicing.

        Not collective.

        Returns
        -------
        inta: float
            The left end of the interval.
        intb: float
            The right end of the interval.

        Notes
        -----
        If the interval was not set by the user, then zeros are returned.
        """
        cdef PetscReal inta = 0
        cdef PetscReal intb = 0
        CHKERR( EPSGetInterval(self.eps, &inta, &intb) )
        return (toReal(inta), toReal(intb))

    def setInterval(self, inta: float, intb: float) -> None:
        """
        Set the computational interval for spectrum slicing.

        Logically collective.

        Parameters
        ----------
        inta
            The left end of the interval.
        intb
            The right end of the interval.

        Notes
        -----
        Spectrum slicing is a technique employed for computing all
        eigenvalues of symmetric eigenproblems in a given interval.
        This function provides the interval to be considered. It must
        be used in combination with `EPS.Which.ALL`, see
        `setWhichEigenpairs()`.
        """
        cdef PetscReal rval1 = asReal(inta)
        cdef PetscReal rval2 = asReal(intb)
        CHKERR( EPSSetInterval(self.eps, rval1, rval2) )

    #

    def getTolerances(self) -> tuple[float, int]:
        """
        Get the tolerance and max. iter. count used for convergence tests.

        Not collective.

        Get the tolerance and iteration limit used by the default EPS
        convergence tests.

        Returns
        -------
        tol: float
            The convergence tolerance.
        max_it: int
            The maximum number of iterations.
        """
        cdef PetscReal rval = 0
        cdef PetscInt  ival = 0
        CHKERR( EPSGetTolerances(self.eps, &rval, &ival) )
        return (toReal(rval), toInt(ival))

    def setTolerances(self, tol: float | None = None, max_it: int | None = None) -> None:
        """
        Set the tolerance and max. iter. used by the default EPS convergence tests.

        Logically collective.

        Parameters
        ----------
        tol
            The convergence tolerance.
        max_it
            The maximum number of iterations.

        Notes
        -----
        Use `DECIDE` for maxits to assign a reasonably good value,
        which is dependent on the solution method.
        """
        cdef PetscReal rval = PETSC_CURRENT
        cdef PetscInt  ival = PETSC_CURRENT
        if tol    is not None: rval = asReal(tol)
        if max_it is not None: ival = asInt(max_it)
        CHKERR( EPSSetTolerances(self.eps, rval, ival) )

    def getTwoSided(self) -> bool:
        """
        Get the flag indicating if a two-sided variant of the algorithm is being used.

        Not collective.

        Returns
        -------
        bool
            Whether the two-sided variant is to be used or not.
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( EPSGetTwoSided(self.eps, &tval) )
        return toBool(tval)

    def setTwoSided(self, twosided: bool) -> None:
        """
        Set to use a two-sided variant that also computes left eigenvectors.

        Logically collective.

        Parameters
        ----------
        twosided
            Whether the two-sided variant is to be used or not.
        """
        cdef PetscBool tval = asBool(twosided)
        CHKERR( EPSSetTwoSided(self.eps, tval) )

    def getPurify(self) -> bool:
        """
        Get the flag indicating whether purification is activated or not.

        Not collective.

        Returns
        -------
        bool
            Whether purification is activated or not.
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( EPSGetPurify(self.eps, &tval) )
        return toBool(tval)

    def setPurify(self, purify: bool = True) -> None:
        """
        Set (toggle) eigenvector purification.

        Logically collective.

        Parameters
        ----------
        purify
            True to activate purification (default).
        """
        cdef PetscBool tval = asBool(purify)
        CHKERR( EPSSetPurify(self.eps, tval) )

    def getConvergenceTest(self) -> Conv:
        """
        Get how to compute the error estimate used in the convergence test.

        Not collective.

        Returns
        -------
        Conv
            The method used to compute the error estimate
            used in the convergence test.
        """
        cdef SlepcEPSConv conv = EPS_CONV_REL
        CHKERR( EPSGetConvergenceTest(self.eps, &conv) )
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
        """
        cdef SlepcEPSConv tconv = conv
        CHKERR( EPSSetConvergenceTest(self.eps, tconv) )

    def getTrueResidual(self) -> bool:
        """
        Get the flag indicating if true residual must be computed explicitly.

        Not collective.

        Returns
        -------
        bool
            Whether the solver compute all residuals or not.
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( EPSGetTrueResidual(self.eps, &tval) )
        return toBool(tval)

    def setTrueResidual(self, trueres: bool) -> None:
        """
        Set if the solver must compute the true residual explicitly or not.

        Logically collective.

        Parameters
        ----------
        trueres
            Whether compute the true residual or not.
        """
        cdef PetscBool tval = asBool(trueres)
        CHKERR( EPSSetTrueResidual(self.eps, tval) )

    def getTrackAll(self) -> bool:
        """
        Get the flag indicating if all residual norms must be computed or not.

        Not collective.

        Returns
        -------
        bool
            Whether the solver compute all residuals or not.
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( EPSGetTrackAll(self.eps, &tval) )
        return toBool(tval)

    def setTrackAll(self, trackall: bool) -> None:
        """
        Set if the solver must compute the residual of all approximate eigenpairs.

        Logically collective.

        Parameters
        ----------
        trackall
            Whether compute all residuals or not.
        """
        cdef PetscBool tval = asBool(trackall)
        CHKERR( EPSSetTrackAll(self.eps, tval) )

    def getDimensions(self) -> tuple[int, int, int]:
        """
        Get number of eigenvalues to compute and the dimension of the subspace.

        Not collective.

        Returns
        -------
        nev: int
            Number of eigenvalues to compute.
        ncv: int
            Maximum dimension of the subspace to be used by the solver.
        mpd: int
            Maximum dimension allowed for the projected problem.
        """
        cdef PetscInt ival1 = 0
        cdef PetscInt ival2 = 0
        cdef PetscInt ival3 = 0
        CHKERR( EPSGetDimensions(self.eps, &ival1, &ival2, &ival3) )
        return (toInt(ival1), toInt(ival2), toInt(ival3))

    def setDimensions(
        self,
        nev: int | None = None,
        ncv: int | None = None,
        mpd: int | None = None,
    ) -> None:
        """
        Set number of eigenvalues to compute and the dimension of the subspace.

        Logically collective.

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
        Use `DECIDE` for ``ncv`` and ``mpd`` to assign a reasonably good
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
        """
        cdef PetscInt ival1 = PETSC_CURRENT
        cdef PetscInt ival2 = PETSC_CURRENT
        cdef PetscInt ival3 = PETSC_CURRENT
        if nev is not None: ival1 = asInt(nev)
        if ncv is not None: ival2 = asInt(ncv)
        if mpd is not None: ival3 = asInt(mpd)
        CHKERR( EPSSetDimensions(self.eps, ival1, ival2, ival3) )

    def getST(self) -> ST:
        """
        Get the spectral transformation object associated to the eigensolver.

        Not collective.

        Returns
        -------
        ST
            The spectral transformation.
        """
        cdef ST st = ST()
        CHKERR( EPSGetST(self.eps, &st.st) )
        CHKERR( PetscINCREF(st.obj) )
        return st

    def setST(self, ST st) -> None:
        """
        Set a spectral transformation object associated to the eigensolver.

        Collective.

        Parameters
        ----------
        st
            The spectral transformation.
        """
        CHKERR( EPSSetST(self.eps, st.st) )

    def getBV(self) -> BV:
        """
        Get the basis vector objects associated to the eigensolver.

        Not collective.

        Returns
        -------
        BV
            The basis vectors context.
        """
        cdef BV bv = BV()
        CHKERR( EPSGetBV(self.eps, &bv.bv) )
        CHKERR( PetscINCREF(bv.obj) )
        return bv

    def setBV(self, BV bv) -> None:
        """
        Set a basis vectors object associated to the eigensolver.

        Collective.

        Parameters
        ----------
        bv
            The basis vectors context.
        """
        CHKERR( EPSSetBV(self.eps, bv.bv) )

    def getDS(self) -> DS:
        """
        Get the direct solver associated to the eigensolver.

        Not collective.

        Returns
        -------
        DS
            The direct solver context.
        """
        cdef DS ds = DS()
        CHKERR( EPSGetDS(self.eps, &ds.ds) )
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
        """
        CHKERR( EPSSetDS(self.eps, ds.ds) )

    def getRG(self) -> RG:
        """
        Get the region object associated to the eigensolver.

        Not collective.

        Returns
        -------
        RG
            The region context.
        """
        cdef RG rg = RG()
        CHKERR( EPSGetRG(self.eps, &rg.rg) )
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
        """
        CHKERR( EPSSetRG(self.eps, rg.rg) )

    def getOperators(self) -> tuple[petsc4py.PETSc.Mat, petsc4py.PETSc.Mat] | tuple[petsc4py.PETSc.Mat, None]:
        """
        Get the matrices associated with the eigenvalue problem.

        Collective.

        Returns
        -------
        A: petsc4py.PETSc.Mat
            The matrix associated with the eigensystem.
        B: petsc4py.PETSc.Mat
            The second matrix in the case of generalized eigenproblems.
        """
        cdef Mat A = Mat()
        cdef Mat B = Mat()
        CHKERR( EPSGetOperators(self.eps, &A.mat, &B.mat) )
        CHKERR( PetscINCREF(A.obj) )
        if B.mat:
            CHKERR( PetscINCREF(B.obj) )
            return (A, B)
        else:
            return (A, None)

    def setOperators(self, Mat A, Mat B=None) -> None:
        """
        Set the matrices associated with the eigenvalue problem.

        Collective.

        Parameters
        ----------
        A
            The matrix associated with the eigensystem.
        B
            The second matrix in the case of generalized eigenproblems;
            if not provided, a standard eigenproblem is assumed.
        """
        cdef PetscMat Bmat = B.mat if B is not None else <PetscMat>NULL
        CHKERR( EPSSetOperators(self.eps, A.mat, Bmat) )

    def setDeflationSpace(self, space: Vec | list[Vec]) -> None:
        """
        Add vectors to the basis of the deflation space.

        Collective.

        Parameters
        ----------
        space
            Set of basis vectors to be added to the deflation space.

        Notes
        -----
        When a deflation space is given, the eigensolver seeks the
        eigensolution in the restriction of the problem to the
        orthogonal complement of this space. This can be used for
        instance in the case that an invariant subspace is known
        beforehand (such as the nullspace of the matrix).

        The vectors do not need to be mutually orthonormal, since they
        are explicitly orthonormalized internally.

        These vectors do not persist from one `solve()` call to the other,
        so the deflation space should be set every time.
        """
        if isinstance(space, Vec): space = [space]
        cdef PetscVec* vs = NULL
        cdef Py_ssize_t i = 0, ns = len(space)
        cdef tmp = allocate(<size_t>ns*sizeof(PetscVec),<void**>&vs)
        for i in range(ns): vs[i] = (<Vec?>space[i]).vec
        CHKERR( EPSSetDeflationSpace(self.eps, <PetscInt>ns, vs) )

    #

    def setInitialSpace(self, space: Vec | list[Vec]) -> None:
        """
        Set the initial space from which the eigensolver starts to iterate.

        Collective.

        Parameters
        ----------
        space
            The initial space

        Notes
        -----
        Some solvers start to iterate on a single vector (initial vector).
        In that case, the other vectors are ignored.

        In contrast to `setDeflationSpace()`, these vectors do not persist
        from one `solve()` call to the other, so the initial space should be
        set every time.

        The vectors do not need to be mutually orthonormal, since they are
        explicitly orthonormalized internally.

        Common usage of this function is when the user can provide a rough
        approximation of the wanted eigenspace. Then, convergence may be faster.
        """
        if isinstance(space, Vec): space = [space]
        cdef PetscVec *vs = NULL
        cdef Py_ssize_t i = 0, ns = len(space)
        cdef tmp = allocate(<size_t>ns*sizeof(PetscVec),<void**>&vs)
        for i in range(ns): vs[i] = (<Vec?>space[i]).vec
        CHKERR( EPSSetInitialSpace(self.eps, <PetscInt>ns, vs) )

    def setLeftInitialSpace(self, space: Vec | list[Vec]) -> None:
        """
        Set a left initial space from which the eigensolver starts to iterate.

        Collective.

        Parameters
        ----------
        space
            The left initial space

        Notes
        -----
        Left initial vectors are used to initiate the left search space
        in two-sided eigensolvers. Users should pass here an approximation
        of the left eigenspace, if available.

        The same comments in `setInitialSpace()` are applicable here.
        """
        if isinstance(space, Vec): space = [space]
        cdef PetscVec *vs = NULL
        cdef Py_ssize_t i = 0, ns = len(space)
        cdef tmp = allocate(<size_t>ns*sizeof(PetscVec),<void**>&vs)
        for i in range(ns): vs[i] = (<Vec?>space[i]).vec
        CHKERR( EPSSetLeftInitialSpace(self.eps, <PetscInt>ns, vs) )

    #

    def setStoppingTest(
        self,
        stopping: EPSStoppingFunction | None,
        args: tuple[Any, ...] | None = None,
        kargs: dict[str, Any] | None = None,
    ) -> None:
        """
        Set when to stop the outer iteration of the eigensolver.

        Logically collective.
        """
        if stopping is not None:
            if args is None: args = ()
            if kargs is None: kargs = {}
            self.set_attr('__stopping__', (stopping, args, kargs))
            CHKERR( EPSSetStoppingTestFunction(self.eps, EPS_Stopping, NULL, NULL) )
        else:
            self.set_attr('__stopping__', None)
            CHKERR( EPSSetStoppingTestFunction(self.eps, EPSStoppingBasic, NULL, NULL) )

    def getStoppingTest(self) -> EPSStoppingFunction:
        """
        Get the stopping function.

        Not collective.
        """
        return self.get_attr('__stopping__')

    #

    def setArbitrarySelection(
        self,
        arbitrary: EPSArbitraryFunction | None,
        args: tuple[Any, ...] | None = None,
        kargs: dict[str, Any] | None = None,
    ) -> None:
        """
        Set an arbitrary selection criterion function.

        Logically collective.

        Set a function to look for eigenvalues according to an arbitrary
        selection criterion. This criterion can be based on a computation
        involving the current eigenvector approximation.
        """
        if arbitrary is not None:
            if args is None: args = ()
            if kargs is None: kargs = {}
            self.set_attr('__arbitrary__', (arbitrary, args, kargs))
            ctx = self.get_attr('__arbitrary__')
            CHKERR( EPSSetArbitrarySelection(self.eps, EPS_Arbitrary, <void*>ctx) )
        else:
            self.set_attr('__arbitrary__', None)
            CHKERR( EPSSetArbitrarySelection(self.eps, NULL, NULL) )

    #

    def setEigenvalueComparison(
        self,
        comparison: EPSEigenvalueComparison | None,
        args: tuple[Any, ...] | None = None,
        kargs: dict[str, Any] | None = None,
    ) -> None:
        """
        Set an eigenvalue comparison function.

        Logically collective.

        Specify the eigenvalue comparison function when `setWhichEigenpairs()`
        is set to `EPS.Which.USER`.
        """
        if comparison is not None:
            if args is None: args = ()
            if kargs is None: kargs = {}
            self.set_attr('__comparison__', (comparison, args, kargs))
            ctx = self.get_attr('__comparison__')
            CHKERR( EPSSetEigenvalueComparison(self.eps, EPS_Comparison, <void*>ctx) )
        else:
            self.set_attr('__comparison__', None)
            CHKERR( EPSSetEigenvalueComparison(self.eps, NULL, NULL) )

    def setMonitor(
        self,
        monitor: EPSMonitorFunction | None,
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
            CHKERR( EPSMonitorSet(self.eps, EPS_Monitor, NULL, NULL) )
        if args is None: args = ()
        if kargs is None: kargs = {}
        monitorlist.append((monitor, args, kargs))

    def getMonitor(self) -> EPSMonitorFunction:
        """Get the list of monitor functions."""
        return self.get_attr('__monitor__')

    def cancelMonitor(self) -> None:
        """
        Clear all monitors for an `EPS` object.

        Logically collective.
        """
        CHKERR( EPSMonitorCancel(self.eps) )
        self.set_attr('__monitor__', None)

    #

    def setUp(self) -> None:
        """
        Set up all the internal data structures.

        Collective.

        Set up all the internal data structures necessary for the execution
        of the eigensolver.

        Notes
        -----
        This function need not be called explicitly in most cases,
        since `solve()` calls it. It can be useful when one wants to
        measure the set-up time separately from the solve time.
        """
        CHKERR( EPSSetUp(self.eps) )

    def solve(self) -> None:
        """
        Solve the eigensystem.

        Collective.
        """
        CHKERR( EPSSolve(self.eps) )

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
        CHKERR( EPSGetIterationNumber(self.eps, &ival) )
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
        cdef SlepcEPSConvergedReason val = EPS_CONVERGED_ITERATING
        CHKERR( EPSGetConvergedReason(self.eps, &val) )
        return val

    def getConverged(self) -> int:
        """
        Get the number of converged eigenpairs.

        Not collective.

        Returns
        -------
        int
            Number of converged eigenpairs.

        Notes
        -----
        This function should be called after `solve()` has finished.
        """
        cdef PetscInt ival = 0
        CHKERR( EPSGetConverged(self.eps, &ival) )
        return toInt(ival)

    def getEigenvalue(self, i: int) -> Scalar:
        """
        Get the i-th eigenvalue as computed by `solve()`.

        Not collective.

        Parameters
        ----------
        i
            Index of the solution to be obtained.

        Returns
        -------
        Scalar
            The computed eigenvalue. It will be a real variable in case
            of a Hermitian or generalized Hermitian eigenproblem. Otherwise
            it will be a complex variable (possibly with zero imaginary part).

        Notes
        -----
        The index ``i`` should be a value between ``0`` and ``nconv-1`` (see
        `getConverged()`). Eigenpairs are indexed according to the ordering
        criterion established with `setWhichEigenpairs()`.
        """
        cdef PetscScalar sval1 = 0
        cdef PetscScalar sval2 = 0
        cdef SlepcEPSProblemType ptype
        CHKERR( EPSGetEigenvalue(self.eps, i, &sval1, &sval2) )
        CHKERR( EPSGetProblemType(self.eps, &ptype) )
        if ptype == EPS_HEP or ptype == EPS_GHEP or ptype == EPS_BSE:
            return toReal(PetscRealPart(sval1))
        else:
            return toComplex(sval1, sval2)

    def getEigenvector(self, i: int, Vec Vr = None, Vec Vi = None) -> None:
        """
        Get the i-th eigenvector as computed by `solve()`.

        Collective.

        Parameters
        ----------
        i
            Index of the solution to be obtained.
        Vr
            Placeholder for the returned eigenvector (real part).
        Vi
            Placeholder for the returned eigenvector (imaginary part).

        Notes
        -----
        The index ``i`` should be a value between ``0`` and
        ``nconv-1`` (see `getConverged()`). Eigenpairs are indexed
        according to the ordering criterion established with
        `setWhichEigenpairs()`.
        """
        cdef PetscVec vecr = Vr.vec if Vr is not None else <PetscVec>NULL
        cdef PetscVec veci = Vi.vec if Vi is not None else <PetscVec>NULL
        CHKERR( EPSGetEigenvector(self.eps, i, vecr, veci) )

    def getLeftEigenvector(self, i: int, Vec Wr = None, Vec Wi = None) -> None:
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
        The index ``i`` should be a value between ``0`` and ``nconv-1`` (see
        `getConverged()`). Eigensolutions are indexed according to the
        ordering criterion established with `setWhichEigenpairs()`.

        Left eigenvectors are available only if the twosided flag was set
        with `setTwoSided()`.
        """
        cdef PetscVec vecr = Wr.vec if Wr is not None else <PetscVec>NULL
        cdef PetscVec veci = Wi.vec if Wi is not None else <PetscVec>NULL
        CHKERR( EPSGetLeftEigenvector(self.eps, i, vecr, veci) )

    def getEigenpair(self, i: int, Vec Vr = None, Vec Vi = None) -> Scalar:
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
        e: Scalar
           The computed eigenvalue. It will be a real variable in case
           of a Hermitian or generalized Hermitian eigenproblem. Otherwise
           it will be a complex variable (possibly with zero imaginary part).

        Notes
        -----
        The index ``i`` should be a value between ``0`` and ``nconv-1`` (see
        `getConverged()`). Eigenpairs are indexed according to the ordering
        criterion established with `setWhichEigenpairs()`.
        """
        cdef PetscScalar sval1 = 0
        cdef PetscScalar sval2 = 0
        cdef PetscVec vecr = Vr.vec if Vr is not None else <PetscVec>NULL
        cdef PetscVec veci = Vi.vec if Vi is not None else <PetscVec>NULL
        cdef SlepcEPSProblemType ptype
        CHKERR( EPSGetEigenpair(self.eps, i, &sval1, &sval2, vecr, veci) )
        CHKERR( EPSGetProblemType(self.eps, &ptype) )
        if ptype == EPS_HEP or ptype == EPS_GHEP or ptype == EPS_BSE:
            return toReal(PetscRealPart(sval1))
        else:
            return toComplex(sval1, sval2)

    def getInvariantSubspace(self) -> list[petsc4py.PETSc.Vec]:
        """
        Get an orthonormal basis of the computed invariant subspace.

        Collective.

        Returns
        -------
        list of petsc4py.PETSc.Vec
            Basis of the invariant subspace.

        Notes
        -----
        This function should be called after `solve()` has finished.

        The returned vectors span an invariant subspace associated
        with the computed eigenvalues. An invariant subspace ``X`` of
        ``A` satisfies ``A x`` in ``X`` for all ``x`` in ``X`` (a
        similar definition applies for generalized eigenproblems).
        """
        cdef PetscInt i = 0, ncv = 0
        cdef PetscVec v = NULL, *isp = NULL
        cdef list subspace = []
        CHKERR( EPSGetConverged(self.eps, &ncv) )
        if ncv == 0: return subspace
        cdef PetscMat A = NULL
        CHKERR( EPSGetOperators(self.eps, &A, NULL) )
        CHKERR( MatCreateVecs(A, &v, NULL) )
        cdef Vec V = None
        cdef object tmp = allocate(<size_t>ncv*sizeof(PetscVec),<void**>&isp)
        for i in range(ncv):
            if i == 0: isp[0] = v
            if i >= 1: CHKERR( VecDuplicate(v, &isp[i]) )
            V = Vec(); V.vec = isp[i]; subspace.append(V)
        CHKERR( EPSGetInvariantSubspace(self.eps, isp) )
        return subspace

    #

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
        This is the error estimate used internally by the
        eigensolver. The actual error bound can be computed with
        `computeError()`.
        """
        cdef PetscReal rval = 0
        CHKERR( EPSGetErrorEstimate(self.eps, i, &rval) )
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
            :math:`\|Ax-kBx\|_2` where :math:`k` is the eigenvalue and
            :math:`x` is the eigenvector.

        Notes
        -----
        The index ``i`` should be a value between ``0`` and
        ``nconv-1`` (see `getConverged()`).
        """
        cdef SlepcEPSErrorType et = EPS_ERROR_RELATIVE
        cdef PetscReal rval = 0
        if etype is not None: et = etype
        CHKERR( EPSComputeError(self.eps, i, et, &rval) )
        return toReal(rval)

    def errorView(self, etype: ErrorType | None = None, Viewer viewer or None: petsc4py.PETSc.Viewer | None = None) -> None:
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
        """
        cdef SlepcEPSErrorType et = EPS_ERROR_RELATIVE
        if etype is not None: et = etype
        cdef PetscViewer vwr = def_Viewer(viewer)
        CHKERR( EPSErrorView(self.eps, et, vwr) )

    def valuesView(self, viewer: Viewer | None = None) -> None:
        """
        Display the computed eigenvalues in a viewer.

        Collective.

        Parameters
        ----------
        viewer
            Visualization context; if not provided, the standard
            output is used.
        """
        cdef PetscViewer vwr = def_Viewer(viewer)
        CHKERR( EPSValuesView(self.eps, vwr) )

    def vectorsView(self, viewer: Viewer | None = None) -> None:
        """
        Output computed eigenvectors to a viewer.

        Collective.

        Parameters
        ----------
        viewer
            Visualization context; if not provided, the standard
            output is used.
        """
        cdef PetscViewer vwr = def_Viewer(viewer)
        CHKERR( EPSVectorsView(self.eps, vwr) )

    #

    def setPowerShiftType(self, shift: PowerShiftType) -> None:
        """
        Set the type of shifts used during the power iteration.

        Logically collective.

        This can be used to emulate the Rayleigh Quotient Iteration (RQI)
        method.

        Parameters
        ----------
        shift
            The type of shift.

        Notes
        -----
        This call is only relevant if the type was set to
        `EPS.Type.POWER` with `setType()`.

        By default, shifts are constant
        (`EPS.PowerShiftType.CONSTANT`) and the iteration is the
        simple power method (or inverse iteration if a
        shift-and-invert transformation is being used).

        A variable shift can be specified
        (`EPS.PowerShiftType.RAYLEIGH` or
        `EPS.PowerShiftType.WILKINSON`). In this case, the iteration
        behaves rather like a cubic converging method as RQI.
        """
        cdef SlepcEPSPowerShiftType val = shift
        CHKERR( EPSPowerSetShiftType(self.eps, val) )

    def getPowerShiftType(self) -> PowerShiftType:
        """
        Get the type of shifts used during the power iteration.

        Not collective.

        Returns
        -------
        PowerShiftType
            The type of shift.
        """
        cdef SlepcEPSPowerShiftType val = EPS_POWER_SHIFT_CONSTANT
        CHKERR( EPSPowerGetShiftType(self.eps, &val) )
        return val

    def setArnoldiDelayed(self, delayed: bool) -> None:
        """
        Set (toggle) delayed reorthogonalization in the Arnoldi iteration.

        Logically collective.

        Parameters
        ----------
        delayed
            True if delayed reorthogonalization is to be used.

        Notes
        -----
        This call is only relevant if the type was set to
        `EPS.Type.ARNOLDI` with `setType()`.

        Delayed reorthogonalization is an aggressive optimization for
        the Arnoldi eigensolver than may provide better scalability,
        but sometimes makes the solver converge less than the default
        algorithm.
        """
        cdef PetscBool val = asBool(delayed)
        CHKERR( EPSArnoldiSetDelayed(self.eps, val) )

    def getArnoldiDelayed(self) -> bool:
        """
        Get the type of reorthogonalization used during the Arnoldi iteration.

        Not collective.

        Returns
        -------
        bool
            True if delayed reorthogonalization is to be used.
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( EPSArnoldiGetDelayed(self.eps, &tval) )
        return toBool(tval)

    def setLanczosReorthogType(self, reorthog: LanczosReorthogType) -> None:
        """
        Set the type of reorthogonalization used during the Lanczos iteration.

        Logically collective.

        Parameters
        ----------
        reorthog
            The type of reorthogonalization.

        Notes
        -----
        This call is only relevant if the type was set to
        `EPS.Type.LANCZOS` with `setType()`.
        """
        cdef SlepcEPSLanczosReorthogType val = reorthog
        CHKERR( EPSLanczosSetReorthog(self.eps, val) )

    def getLanczosReorthogType(self) -> LanczosReorthogType:
        """
        Get the type of reorthogonalization used during the Lanczos iteration.

        Not collective.

        Returns
        -------
        LanczosReorthogType
            The type of reorthogonalization.
        """
        cdef SlepcEPSLanczosReorthogType val = \
            EPS_LANCZOS_REORTHOG_LOCAL
        CHKERR( EPSLanczosGetReorthog(self.eps, &val) )
        return val

    #

    def setKrylovSchurBSEType(self, bse: KrylovSchurBSEType) -> None:
        """
        Set the Krylov-Schur variant used for BSE structured eigenproblems.

        Logically collective.

        Parameters
        ----------
        bse
            The BSE method.

        Notes
        -----
        This call is only relevant if the type was set to
        `EPS.Type.KRYLOVSCHUR` with `setType()` and the problem
        type to `EPS.ProblemType.BSE` with `setProblemType()`.
        """
        cdef SlepcEPSKrylovSchurBSEType val = bse
        CHKERR( EPSKrylovSchurSetBSEType(self.eps, val) )

    def getKrylovSchurBSEType(self) -> KrylovSchurBSEType:
        """
        Get the method used for BSE structured eigenproblems (Krylov-Schur).

        Not collective.

        Returns
        -------
        KrylovSchurBSEType
            The BSE method.
        """
        cdef SlepcEPSKrylovSchurBSEType val = EPS_KRYLOVSCHUR_BSE_SHAO
        CHKERR( EPSKrylovSchurGetBSEType(self.eps, &val) )
        return val

    def setKrylovSchurRestart(self, keep: float) -> None:
        """
        Set the restart parameter for the Krylov-Schur method.

        Logically collective.

        It is the proportion of basis vectors that must be kept after restart.

        Parameters
        ----------
        keep
            The number of vectors to be kept at restart.

        Notes
        -----
        Allowed values are in the range [0.1,0.9]. The default is 0.5.
        """
        cdef PetscReal val = asReal(keep)
        CHKERR( EPSKrylovSchurSetRestart(self.eps, val) )

    def getKrylovSchurRestart(self) -> float:
        """
        Get the restart parameter used in the Krylov-Schur method.

        Not collective.

        Returns
        -------
        float
            The number of vectors to be kept at restart.
        """
        cdef PetscReal val = 0
        CHKERR( EPSKrylovSchurGetRestart(self.eps, &val) )
        return toReal(val)

    def setKrylovSchurLocking(self, lock: bool) -> None:
        """
        Set (toggle) locking/non-locking variants of the Krylov-Schur method.

        Logically collective.

        Parameters
        ----------
        lock
            True if the locking variant must be selected.

        Notes
        -----
        The default is to lock converged eigenpairs when the method restarts.
        This behavior can be changed so that all directions are kept in the
        working subspace even if already converged to working accuracy (the
        non-locking variant).
        """
        cdef PetscBool val = asBool(lock)
        CHKERR( EPSKrylovSchurSetLocking(self.eps, val) )

    def getKrylovSchurLocking(self) -> bool:
        """
        Get the locking flag used in the Krylov-Schur method.

        Not collective.

        Returns
        -------
        bool
            The locking flag.
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( EPSKrylovSchurGetLocking(self.eps, &tval) )
        return toBool(tval)

    def setKrylovSchurPartitions(self, npart: int) -> None:
        """
        Set the number of partitions of the communicator (spectrum slicing).

        Logically collective.

        Set the number of partitions for the case of doing spectrum
        slicing for a computational interval with the communicator split
        in several sub-communicators.

        Parameters
        ----------
        npart
            The number of partitions.

        Notes
        -----
        By default, npart=1 so all processes in the communicator participate in
        the processing of the whole interval. If npart>1 then the interval is
        divided into npart subintervals, each of them being processed by a
        subset of processes.
        """
        cdef PetscInt val = asInt(npart)
        CHKERR( EPSKrylovSchurSetPartitions(self.eps, val) )

    def getKrylovSchurPartitions(self) -> int:
        """
        Get the number of partitions of the communicator (spectrum slicing).

        Not collective.

        Returns
        -------
        int
            The number of partitions.
        """
        cdef PetscInt val = 0
        CHKERR( EPSKrylovSchurGetPartitions(self.eps, &val) )
        return toInt(val)

    def setKrylovSchurDetectZeros(self, detect: bool) -> None:
        """
        Set the flag that enforces zero detection in spectrum slicing.

        Logically collective.

        Set a flag to enforce the detection of zeros during the factorizations
        throughout the spectrum slicing computation.

        Parameters
        ----------
        detect
            True if zeros must checked for.

        Notes
        -----
        A zero in the factorization indicates that a shift coincides with
        an eigenvalue.

        This flag is turned off by default, and may be necessary in some cases,
        especially when several partitions are being used. This feature currently
        requires an external package for factorizations with support for zero
        detection, e.g. MUMPS.
        """
        cdef PetscBool val = asBool(detect)
        CHKERR( EPSKrylovSchurSetDetectZeros(self.eps, val) )

    def getKrylovSchurDetectZeros(self) -> bool:
        """
        Get the flag that enforces zero detection in spectrum slicing.

        Not collective.

        Returns
        -------
        bool
            The zero detection flag.
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( EPSKrylovSchurGetDetectZeros(self.eps, &tval) )
        return toBool(tval)

    def setKrylovSchurDimensions(
        self,
        nev: int | None = None,
        ncv: int | None = None,
        mpd: int | None = None,
    ) -> None:
        """
        Set the dimensions used for each subsolve step (spectrum slicing).

        Logically collective.

        Set the dimensions used for each subsolve step in case of doing
        spectrum slicing for a computational interval. The meaning of the
        parameters is the same as in `setDimensions()`.

        Parameters
        ----------
        nev
            Number of eigenvalues to compute.
        ncv
            Maximum dimension of the subspace to be used by the solver.
        mpd
            Maximum dimension allowed for the projected problem.
        """
        cdef PetscInt ival1 = PETSC_CURRENT
        cdef PetscInt ival2 = PETSC_CURRENT
        cdef PetscInt ival3 = PETSC_CURRENT
        if nev is not None: ival1 = asInt(nev)
        if ncv is not None: ival2 = asInt(ncv)
        if mpd is not None: ival3 = asInt(mpd)
        CHKERR( EPSKrylovSchurSetDimensions(self.eps, ival1, ival2, ival3) )

    def getKrylovSchurDimensions(self) -> tuple[int, int, int]:
        """
        Get the dimensions used for each subsolve step (spectrum slicing).

        Not collective.

        Get the dimensions used for each subsolve step in case of doing
        spectrum slicing for a computational interval.

        Returns
        -------
        nev: int
            Number of eigenvalues to compute.
        ncv: int
            Maximum dimension of the subspace to be used by the solver.
        mpd: int
            Maximum dimension allowed for the projected problem.
        """
        cdef PetscInt ival1 = 0
        cdef PetscInt ival2 = 0
        cdef PetscInt ival3 = 0
        CHKERR( EPSKrylovSchurGetDimensions(self.eps, &ival1, &ival2, &ival3) )
        return (toInt(ival1), toInt(ival2), toInt(ival3))

    def getKrylovSchurSubcommInfo(self) -> tuple[int, int, petsc4py.PETSc.Vec]:
        """
        Get information related to the case of doing spectrum slicing.

        Collective on the subcommunicator (if v is given).

        Get information related to the case of doing spectrum slicing
        for a computational interval with multiple communicators.

        Returns
        -------
        k: int
            Number of the subinterval for the calling process.
        n: int
            Number of eigenvalues found in the k-th subinterval.
        v: petsc4py.PETSc.Vec
            A vector owned by processes in the subcommunicator with dimensions
            compatible for locally computed eigenvectors.

        Notes
        -----
        This function is only available for spectrum slicing runs.

        The returned Vec should be destroyed by the user.
        """
        cdef PetscInt ival1 = 0
        cdef PetscInt ival2 = 0
        cdef Vec vec = Vec()
        CHKERR( EPSKrylovSchurGetSubcommInfo(self.eps, &ival1, &ival2, &vec.vec) )
        return (toInt(ival1), toInt(ival2), vec)

    def getKrylovSchurSubcommPairs(self, i: int, Vec V) -> Scalar:
        """
        Get the i-th eigenpair stored in the multi-communicator of the process.

        Collective on the subcommunicator (if v is given).

        Get the i-th eigenpair stored internally in the multi-communicator
        to which the calling process belongs.

        Parameters
        ----------
        i
            Index of the solution to be obtained.
        V
            Placeholder for the returned eigenvector.

        Returns
        -------
        Scalar
            The computed eigenvalue.

        Notes
        -----
        The index ``i`` should be a value between ``0`` and ``n-1``,
        where ``n`` is the number of vectors in the local subinterval,
        see `getKrylovSchurSubcommInfo()`.
        """
        cdef PetscScalar sval = 0
        cdef PetscVec vec = V.vec if V is not None else <PetscVec>NULL
        CHKERR( EPSKrylovSchurGetSubcommPairs(self.eps, i, &sval, vec) )
        return toScalar(sval)

    def getKrylovSchurSubcommMats(self) -> tuple[petsc4py.PETSc.Mat, petsc4py.PETSc.Mat] | tuple[petsc4py.PETSc.Mat, None]:
        """
        Get the eigenproblem matrices stored in the subcommunicator.

        Collective on the subcommunicator.

        Get the eigenproblem matrices stored internally in the subcommunicator
        to which the calling process belongs.

        Returns
        -------
        A: petsc4py.PETSc.Mat
            The matrix associated with the eigensystem.
        B: petsc4py.PETSc.Mat
            The second matrix in the case of generalized eigenproblems.

        Notes
        -----
        This is the analog of `getOperators()`, but returns the matrices distributed
        differently (in the subcommunicator rather than in the parent communicator).

        These matrices should not be modified by the user.
        """
        cdef Mat A = Mat()
        cdef Mat B = Mat()
        CHKERR( EPSKrylovSchurGetSubcommMats(self.eps, &A.mat, &B.mat) )
        CHKERR( PetscINCREF(A.obj) )
        if B.mat:
            CHKERR( PetscINCREF(B.obj) )
            return (A, B)
        else:
            return (A, None)

    def updateKrylovSchurSubcommMats(
        self,
        s: Scalar = 1.0,
        a: Scalar = 1.0,
        Mat Au: petsc4py.PETSc.Mat | None = None,
        t: Scalar = 1.0,
        b: Scalar = 1.0,
        Mat Bu: petsc4py.PETSc.Mat | None = None,
        structure: petsc4py.PETSc.Mat.Structure | None = None,
        globalup: bool = False,
    ) -> None:
        """
        Update the eigenproblem matrices stored internally in the communicator.

        Collective.

        Update the eigenproblem matrices stored internally in the
        subcommunicator to which the calling process belongs.

        Parameters
        ----------
        s
            Scalar that multiplies the existing A matrix.
        a
            Scalar used in the axpy operation on A.
        Au
            The matrix used in the axpy operation on A.
        t
            Scalar that multiplies the existing B matrix.
        b
            Scalar used in the axpy operation on B.
        Bu
            The matrix used in the axpy operation on B.
        structure
            Either same, different, or a subset of the non-zero sparsity pattern.
        globalup
            Whether global matrices must be updated or not.

        Notes
        -----
        This function modifies the eigenproblem matrices at subcommunicator
        level, and optionally updates the global matrices in the parent
        communicator.  The updates are expressed as
        :math:`A \leftarrow s A + a Au`,
        :math:`B \leftarrow t B + b Bu`.

        It is possible to update one of the matrices, or both.

        The matrices ``Au`` and ``Bu`` must be equal in all subcommunicators.

        The ``structure`` flag is passed to the `petsc4py.PETSc.Mat.axpy`
        operations to perform the updates.

        If ``globalup`` is True, communication is carried out to reconstruct
        the updated matrices in the parent communicator.

        """
        cdef PetscMat Amat = Au.mat if Au is not None else <PetscMat>NULL
        cdef PetscMat Bmat = Bu.mat if Bu is not None else <PetscMat>NULL
        cdef PetscMatStructure vstr = matstructure(structure)
        cdef PetscBool tval = globalup
        CHKERR( EPSKrylovSchurUpdateSubcommMats(self.eps, s, a, Amat, t, b, Bmat, vstr, tval) )

    def setKrylovSchurSubintervals(self, subint: Sequence[float]) -> None:
        """
        Set the subinterval boundaries.

        Logically collective.

        Set the subinterval boundaries for spectrum slicing with a
        computational interval.

        Parameters
        ----------
        subint
            Real values specifying subintervals

        Notes
        -----
        This function must be called after setKrylovSchurPartitions().
        For npart partitions, the argument subint must contain npart+1
        real values sorted in ascending order:
        subint_0, subint_1, ..., subint_npart,
        where the first and last values must coincide with the interval
        endpoints set with EPSSetInterval().
        The subintervals are then defined by two consecutive points:
        [subint_0,subint_1], [subint_1,subint_2], and so on.
        """
        cdef PetscBool match = PETSC_FALSE
        CHKERR( PetscObjectTypeCompare(<PetscObject>self.eps, EPSKRYLOVSCHUR, &match) )
        if match == PETSC_FALSE: return
        cdef PetscReal *subintarray = NULL
        cdef Py_ssize_t i = 0, n = len(subint)
        cdef PetscInt nparts = 0
        CHKERR( EPSKrylovSchurGetPartitions(self.eps, &nparts) )
        assert n >= nparts
        cdef tmp = allocate(<size_t>n*sizeof(PetscReal),<void**>&subintarray)
        for i in range(n): subintarray[i] = asReal(subint[i])
        CHKERR( EPSKrylovSchurSetSubintervals(self.eps, subintarray) )

    def getKrylovSchurSubintervals(self) -> ArrayReal:
        """
        Get the points that delimit the subintervals.

        Not collective.

        Get the points that delimit the subintervals used in spectrum slicing
        with several partitions.

        Returns
        -------
        ArrayReal
            Real values specifying subintervals
        """
        cdef PetscReal *subintarray = NULL
        cdef PetscInt nparts = 0
        CHKERR( EPSKrylovSchurGetPartitions(self.eps, &nparts) )
        CHKERR( EPSKrylovSchurGetSubintervals(self.eps, &subintarray) )
        cdef object subint = None
        try:
            subint = array_r(nparts+1, subintarray)
        finally:
            CHKERR( PetscFree(subintarray) )
        return subint

    def getKrylovSchurInertias(self) -> tuple[ArrayReal, ArrayInt]:
        """
        Get the values of the shifts and their corresponding inertias.

        Not collective.

        Get the values of the shifts and their corresponding inertias in case
        of doing spectrum slicing for a computational interval.

        Returns
        -------
        shifts: ArrayReal
            The values of the shifts used internally in the solver.
        inertias: ArrayInt
            The values of the inertia in each shift.
        """
        cdef PetscReal *shiftsarray = NULL
        cdef PetscInt *inertiasarray = NULL
        cdef PetscInt n = 0
        CHKERR(EPSKrylovSchurGetInertias(self.eps, &n, &shiftsarray, &inertiasarray))
        cdef object shifts = None
        cdef object inertias = None
        try:
            shifts = array_r(n, shiftsarray)
            inertias = array_i(n, inertiasarray)
        finally:
            CHKERR( PetscFree(shiftsarray) )
            CHKERR( PetscFree(inertiasarray) )
        return (shifts, inertias)

    def getKrylovSchurKSP(self) -> KSP:
        """
        Get the linear solver object associated with the internal `EPS` object.

        Collective.

        Get the linear solver object associated with the internal `EPS`
        object in case of doing spectrum slicing for a computational interval.

        Returns
        -------
        `petsc4py.PETSc.KSP`
            The linear solver object.
        """
        cdef KSP ksp = KSP()
        CHKERR( EPSKrylovSchurGetKSP(self.eps, &ksp.ksp) )
        CHKERR( PetscINCREF(ksp.obj) )
        return ksp

    #

    def setGDKrylovStart(self, krylovstart: bool = True) -> None:
        """
        Set (toggle) starting the search subspace with a Krylov basis.

        Logically collective.

        Parameters
        ----------
        krylovstart
            True if starting the search subspace with a Krylov basis.
        """
        cdef PetscBool val = asBool(krylovstart)
        CHKERR( EPSGDSetKrylovStart(self.eps, val) )

    def getGDKrylovStart(self) -> bool:
        """
        Get a flag indicating if the search subspace is started with a Krylov basis.

        Not collective.

        Returns
        -------
        bool
            True if starting the search subspace with a Krylov basis.
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( EPSGDGetKrylovStart(self.eps, &tval) )
        return toBool(tval)

    def setGDBlockSize(self, bs: int) -> None:
        """
        Set the number of vectors to be added to the searching space.

        Logically collective.

        Set the number of vectors to be added to the searching space in every
        iteration.

        Parameters
        ----------
        bs
            The number of vectors added to the search space in every iteration.
        """
        cdef PetscInt ival = asInt(bs)
        CHKERR( EPSGDSetBlockSize(self.eps, ival) )

    def getGDBlockSize(self) -> int:
        """
        Get the number of vectors to be added to the searching space.

        Not collective.

        Get the number of vectors to be added to the searching space in every
        iteration.

        Returns
        -------
        int
            The number of vectors added to the search space in every iteration.
        """
        cdef PetscInt ival = 0
        CHKERR( EPSGDGetBlockSize(self.eps, &ival) )
        return toInt(ival)

    def setGDRestart(self, minv: int = None, plusk: int = None) -> None:
        """
        Set the number of vectors of the search space after restart.

        Logically collective.

        Set the number of vectors of the search space after restart and
        the number of vectors saved from the previous iteration.

        Parameters
        ----------
        minv
            The number of vectors of the search subspace after restart.
        plusk
            The number of vectors saved from the previous iteration.
        """
        cdef PetscInt ival1 = PETSC_CURRENT
        cdef PetscInt ival2 = PETSC_CURRENT
        if minv  is not None: ival1 = asInt(minv)
        if plusk is not None: ival2 = asInt(plusk)
        CHKERR( EPSGDSetRestart(self.eps, ival1, ival2) )

    def getGDRestart(self) -> tuple[int, int]:
        """
        Get the number of vectors of the search space after restart.

        Not collective.

        Get the number of vectors of the search space after restart and
        the number of vectors saved from the previous iteration.

        Returns
        -------
        minv: int
            The number of vectors of the search subspace after restart.
        plusk: int
            The number of vectors saved from the previous iteration.
        """
        cdef PetscInt ival1 = 0
        cdef PetscInt ival2 = 0
        CHKERR( EPSGDGetRestart(self.eps, &ival1, &ival2) )
        return (toInt(ival1), toInt(ival2))

    def setGDInitialSize(self, initialsize: int) -> None:
        """
        Set the initial size of the searching space.

        Logically collective.

        Parameters
        ----------
        initialsize
            The number of vectors of the initial searching subspace.
        """
        cdef PetscInt ival = asInt(initialsize)
        CHKERR( EPSGDSetInitialSize(self.eps, ival) )

    def getGDInitialSize(self) -> int:
        """
        Get the initial size of the searching space.

        Not collective.

        Returns
        -------
        int
            The number of vectors of the initial searching subspace.
        """
        cdef PetscInt ival = 0
        CHKERR( EPSGDGetInitialSize(self.eps, &ival) )
        return toInt(ival)

    def setGDBOrth(self, borth: bool) -> int:
        """
        Set the orthogonalization that will be used in the search subspace.

        Logically collective.

        Set the orthogonalization that will be used in the search
        subspace in case of generalized Hermitian problems.

        Parameters
        ----------
        borth
            Whether to B-orthogonalize the search subspace.
        """
        cdef PetscBool tval = asBool(borth)
        CHKERR( EPSGDSetBOrth(self.eps, tval) )

    def getGDBOrth(self) -> bool:
        """
        Get the orthogonalization used in the search subspace.

        Not collective.

        Get the orthogonalization used in the search subspace in
        case of generalized Hermitian problems.

        Returns
        -------
        bool
            Whether to B-orthogonalize the search subspace.
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( EPSGDGetBOrth(self.eps, &tval) )
        return toBool(tval)

    def setGDDoubleExpansion(self, doubleexp: bool) -> None:
        """
        Set that the search subspace is expanded with double expansion.

        Logically collective.

        Set a variant where the search subspace is expanded with
        :math:`K [A x, B x]` (double expansion) instead of the
        classic :math:`K r`, where K is the preconditioner, x the
        selected approximate eigenvector and :math:`r` its associated residual
        vector.

        Parameters
        ----------
        doubleexp
            True if using double expansion.
        """
        cdef PetscBool val = asBool(doubleexp)
        CHKERR( EPSGDSetDoubleExpansion(self.eps, val) )

    def getGDDoubleExpansion(self) -> bool:
        """
        Get a flag indicating whether the double expansion variant is active.

        Not collective.

        Get a flag indicating whether the double expansion variant
        has been activated or not.

        Returns
        -------
        bool
            True if using double expansion.
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( EPSGDGetDoubleExpansion(self.eps, &tval) )
        return toBool(tval)

    #

    def setJDKrylovStart(self, krylovstart: bool = True) -> None:
        """
        Set (toggle) starting the search subspace with a Krylov basis.

        Logically collective.

        Parameters
        ----------
        krylovstart
            True if starting the search subspace with a Krylov basis.
        """
        cdef PetscBool val = asBool(krylovstart)
        CHKERR( EPSJDSetKrylovStart(self.eps, val) )

    def getJDKrylovStart(self) -> bool:
        """
        Get a flag indicating if the search subspace is started with a Krylov basis.

        Not collective.

        Returns
        -------
        bool
            True if starting the search subspace with a Krylov basis.
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( EPSJDGetKrylovStart(self.eps, &tval) )
        return toBool(tval)

    def setJDBlockSize(self, bs: int) -> None:
        """
        Set the number of vectors to be added to the searching space.

        Logically collective.

        Set the number of vectors to be added to the searching space in every
        iteration.

        Parameters
        ----------
        bs
            The number of vectors added to the search space in every iteration.
        """
        cdef PetscInt ival = asInt(bs)
        CHKERR( EPSJDSetBlockSize(self.eps, ival) )

    def getJDBlockSize(self) -> int:
        """
        Get the number of vectors to be added to the searching space.

        Not collective.

        Get the number of vectors to be added to the searching space in every
        iteration.

        Returns
        -------
        int
            The number of vectors added to the search space in every iteration.
        """
        cdef PetscInt ival = 0
        CHKERR( EPSJDGetBlockSize(self.eps, &ival) )
        return toInt(ival)

    def setJDRestart(self, minv: int | None = None, plusk: int | None = None) -> None:
        """
        Set the number of vectors of the search space after restart.

        Logically collective.

        Set the number of vectors of the search space after restart and
        the number of vectors saved from the previous iteration.

        Parameters
        ----------
        minv
            The number of vectors of the search subspace after restart.
        plusk
            The number of vectors saved from the previous iteration.
        """
        cdef PetscInt ival1 = PETSC_CURRENT
        cdef PetscInt ival2 = PETSC_CURRENT
        if minv  is not None: ival1 = asInt(minv)
        if plusk is not None: ival2 = asInt(plusk)
        CHKERR( EPSJDSetRestart(self.eps, ival1, ival2) )

    def getJDRestart(self) -> tuple[int, int]:
        """
        Get the number of vectors of the search space after restart.

        Not collective.

        Get the number of vectors of the search space after restart and
        the number of vectors saved from the previous iteration.

        Returns
        -------
        minv: int
            The number of vectors of the search subspace after restart.
        plusk: int
            The number of vectors saved from the previous iteration.
        """
        cdef PetscInt ival1 = 0
        cdef PetscInt ival2 = 0
        CHKERR( EPSJDGetRestart(self.eps, &ival1, &ival2) )
        return (toInt(ival1), toInt(ival2))

    def setJDInitialSize(self, initialsize: int) -> None:
        """
        Set the initial size of the searching space.

        Logically collective.

        Parameters
        ----------
        initialsize
            The number of vectors of the initial searching subspace.
        """
        cdef PetscInt ival = asInt(initialsize)
        CHKERR( EPSJDSetInitialSize(self.eps, ival) )

    def getJDInitialSize(self) -> int:
        """
        Get the initial size of the searching space.

        Not collective.

        Returns
        -------
        int
            The number of vectors of the initial searching subspace.
        """
        cdef PetscInt ival = 0
        CHKERR( EPSJDGetInitialSize(self.eps, &ival) )
        return toInt(ival)

    def setJDFix(self, fix: float) -> None:
        """
        Set the threshold for changing the target in the correction equation.

        Logically collective.

        Parameters
        ----------
        fix
            The threshold for changing the target.

        Notes
        -----
        The target in the correction equation is fixed at the first iterations.
        When the norm of the residual vector is lower than the fix value,
        the target is set to the corresponding eigenvalue.
        """
        cdef PetscReal val = asReal(fix)
        CHKERR( EPSJDSetFix(self.eps, val) )

    def getJDFix(self) -> float:
        """
        Get the threshold for changing the target in the correction equation.

        Not collective.

        Returns
        -------
        float
            The threshold for changing the target.
        """
        cdef PetscReal val = 0
        CHKERR( EPSJDGetFix(self.eps, &val) )
        return toReal(val)

    def setJDConstCorrectionTol(self, constant: bool) -> None:
        """
        Deactivate the dynamic stopping criterion.

        Logically collective.

        Deactivate the dynamic stopping criterion that sets the
        `petsc4py.PETSc.KSP` relative tolerance to ``0.5**i``, where ``i`` is
        the number of `EPS` iterations from the last converged value.

        Parameters
        ----------
        constant
            If False, the `petsc4py.PETSc.KSP` relative tolerance is set to ``0.5**i``.
        """
        cdef PetscBool tval = asBool(constant)
        CHKERR( EPSJDSetConstCorrectionTol(self.eps, tval) )

    def getJDConstCorrectionTol(self) -> bool:
        """
        Get the flag indicating if the dynamic stopping is being used.

        Not collective.

        Get the flag indicating if the dynamic stopping is being used for
        solving the correction equation.

        Returns
        -------
        bool
            True if the dynamic stopping criterion is not being used.
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( EPSJDGetConstCorrectionTol(self.eps, &tval) )
        return toBool(tval)

    def setJDBOrth(self, borth: bool) -> None:
        """
        Set the orthogonalization that will be used in the search subspace.

        Logically collective.

        Set the orthogonalization that will be used in the search
        subspace in case of generalized Hermitian problems.

        Parameters
        ----------
        borth
            Whether to B-orthogonalize the search subspace.
        """
        cdef PetscBool tval = asBool(borth)
        CHKERR( EPSJDSetBOrth(self.eps, tval) )

    def getJDBOrth(self) -> bool:
        """
        Get the orthogonalization used in the search subspace.

        Not collective.

        Get the orthogonalization used in the search subspace in
        case of generalized Hermitian problems.

        Returns
        -------
        bool
            Whether to B-orthogonalize the search subspace.
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( EPSJDGetBOrth(self.eps, &tval) )
        return toBool(tval)

    #

    def setRQCGReset(self, nrest: int) -> None:
        """
        Set the reset parameter of the RQCG iteration.

        Logically collective.

        Every nrest iterations the solver performs a Rayleigh-Ritz projection
        step.

        Parameters
        ----------
        nrest
            The number of iterations between resets.
        """
        cdef PetscInt val = asInt(nrest)
        CHKERR( EPSRQCGSetReset(self.eps, val) )

    def getRQCGReset(self) -> int:
        """
        Get the reset parameter used in the RQCG method.

        Not collective.

        Returns
        -------
        int
            The number of iterations between resets.
        """
        cdef PetscInt val = 0
        CHKERR( EPSRQCGGetReset(self.eps, &val) )
        return toInt(val)

    def setLOBPCGBlockSize(self, bs: int) -> None:
        """
        Set the block size of the LOBPCG method.

        Logically collective.

        Parameters
        ----------
        bs
            The block size.
        """
        cdef PetscInt ival = asInt(bs)
        CHKERR( EPSLOBPCGSetBlockSize(self.eps, ival) )

    def getLOBPCGBlockSize(self) -> int:
        """
        Get the block size used in the LOBPCG method.

        Not collective.

        Returns
        -------
        int
            The block size.
        """
        cdef PetscInt ival = 0
        CHKERR( EPSLOBPCGGetBlockSize(self.eps, &ival) )
        return toInt(ival)

    def setLOBPCGRestart(self, restart: float) -> None:
        """
        Set the restart parameter for the LOBPCG method.

        Logically collective.

        The meaning of this parameter is the proportion of vectors within the
        current block iterate that must have converged in order to force a
        restart with hard locking.

        Parameters
        ----------
        restart
            The percentage of the block of vectors to force a restart.

        Notes
        -----
        Allowed values are in the range [0.1,1.0]. The default is 0.9.
        """
        cdef PetscReal val = asReal(restart)
        CHKERR( EPSLOBPCGSetRestart(self.eps, val) )

    def getLOBPCGRestart(self) -> float:
        """
        Get the restart parameter used in the LOBPCG method.

        Not collective.

        Returns
        -------
        float
            The restart parameter.
        """
        cdef PetscReal val = 0
        CHKERR( EPSLOBPCGGetRestart(self.eps, &val) )
        return toReal(val)

    def setLOBPCGLocking(self, lock: bool) -> None:
        """
        Toggle between locking and non-locking (LOBPCG method).

        Logically collective.

        Parameters
        ----------
        lock
            True if the locking variant must be selected.

        Notes
        -----
        This flag refers to soft locking (converged vectors within the current
        block iterate), since hard locking is always used (when nev is larger
        than the block size).
        """
        cdef PetscBool val = asBool(lock)
        CHKERR( EPSLOBPCGSetLocking(self.eps, val) )

    def getLOBPCGLocking(self) -> bool:
        """
        Get the locking flag used in the LOBPCG method.

        Not collective.

        Returns
        -------
        bool
            The locking flag.
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( EPSLOBPCGGetLocking(self.eps, &tval) )
        return toBool(tval)

    def setLyapIIRanks(self, rkc: int | None = None, rkl: int | None = None) -> None:
        """
        Set the ranks used in the solution of the Lyapunov equation.

        Logically collective.

        Parameters
        ----------
        rkc
            The compressed rank.
        rkl
            The Lyapunov rank.
        """
        cdef PetscInt ival1 = PETSC_CURRENT
        cdef PetscInt ival2 = PETSC_CURRENT
        if rkc  is not None: ival1 = asInt(rkc)
        if rkl is not None: ival2 = asInt(rkl)
        CHKERR( EPSLyapIISetRanks(self.eps, ival1, ival2) )

    def getLyapIIRanks(self) -> tuple[int, int]:
        """
        Get the rank values used for the Lyapunov step.

        Not collective.

        Returns
        -------
        rkc: int
            The compressed rank.
        rkl: int
            The Lyapunov rank.
        """
        cdef PetscInt ival1 = 0
        cdef PetscInt ival2 = 0
        CHKERR( EPSLyapIIGetRanks(self.eps, &ival1, &ival2) )
        return (toInt(ival1), toInt(ival2))

    #

    def setCISSExtraction(self, extraction: CISSExtraction) -> None:
        """
        Set the extraction technique used in the CISS solver.

        Logically collective.

        Parameters
        ----------
        extraction
            The extraction technique.
        """
        cdef SlepcEPSCISSExtraction val = extraction
        CHKERR( EPSCISSSetExtraction(self.eps, val) )

    def getCISSExtraction(self) -> CISSExtraction:
        """
        Get the extraction technique used in the CISS solver.

        Not collective.

        Returns
        -------
        CISSExtraction
            The extraction technique.
        """
        cdef SlepcEPSCISSExtraction val = EPS_CISS_EXTRACTION_RITZ
        CHKERR( EPSCISSGetExtraction(self.eps, &val) )
        return val

    def setCISSQuadRule(self, quad: CISSQuadRule) -> None:
        """
        Set the quadrature rule used in the CISS solver.

        Logically collective.

        Parameters
        ----------
        quad
            The quadrature rule.
        """
        cdef SlepcEPSCISSQuadRule val = quad
        CHKERR( EPSCISSSetQuadRule(self.eps, val) )

    def getCISSQuadRule(self) -> CISSQuadRule:
        """
        Get the quadrature rule used in the CISS solver.

        Not collective.

        Returns
        -------
        CISSQuadRule
            The quadrature rule.
        """
        cdef SlepcEPSCISSQuadRule val = EPS_CISS_QUADRULE_TRAPEZOIDAL
        CHKERR( EPSCISSGetQuadRule(self.eps, &val) )
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
            True if A and B are real.

        Notes
        -----
        The default number of partitions is 1. This means the internal
        `petsc4py.PETSc.KSP` object is shared among all processes of the
        `EPS` communicator. Otherwise, the communicator is split into npart
        communicators, so that ``npart`` `petsc4py.PETSc.KSP` solves proceed
        simultaneously.
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
        CHKERR( EPSCISSSetSizes(self.eps, ival1, ival2, ival3, ival4, ival5, bval) )

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
            True if A and B are real.
        """
        cdef PetscInt  ival1 = 0
        cdef PetscInt  ival2 = 0
        cdef PetscInt  ival3 = 0
        cdef PetscInt  ival4 = 0
        cdef PetscInt  ival5 = 0
        cdef PetscBool bval  = PETSC_FALSE
        CHKERR( EPSCISSGetSizes(self.eps, &ival1, &ival2, &ival3, &ival4, &ival5, &bval) )
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
        """
        cdef PetscReal rval1 = PETSC_CURRENT
        cdef PetscReal rval2 = PETSC_CURRENT
        if delta is not None: rval1 = asReal(delta)
        if spur  is not None: rval2 = asReal(spur)
        CHKERR( EPSCISSSetThreshold(self.eps, rval1, rval2) )

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
        """
        cdef PetscReal delta = 0
        cdef PetscReal spur  = 0
        CHKERR( EPSCISSGetThreshold(self.eps, &delta, &spur) )
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
        """
        cdef PetscInt ival1 = PETSC_CURRENT
        cdef PetscInt ival2 = PETSC_CURRENT
        if inner  is not None: ival1 = asInt(inner)
        if blsize is not None: ival2 = asInt(blsize)
        CHKERR( EPSCISSSetRefinement(self.eps, ival1, ival2) )

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
        """
        cdef PetscInt ival1 = 0
        cdef PetscInt ival2 = 0
        CHKERR( EPSCISSGetRefinement(self.eps, &ival1, &ival2) )
        return (toInt(ival1), toInt(ival2))

    def setCISSUseST(self, usest: bool) -> None:
        """
        Set a flag indicating that the CISS solver will use the `ST` object.

        Logically collective.

        Set a flag indicating that the CISS solver will use the `ST`
        object for the linear solves.

        Parameters
        ----------
        usest
            Whether to use the `ST` object or not.
        """
        cdef PetscBool tval = asBool(usest)
        CHKERR( EPSCISSSetUseST(self.eps, tval) )

    def getCISSUseST(self) -> bool:
        """
        Get the flag indicating the use of the `ST` object in the CISS solver.

        Not collective.

        Returns
        -------
        bool
            Whether to use the `ST` object or not.
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( EPSCISSGetUseST(self.eps, &tval) )
        return toBool(tval)

    def getCISSKSPs(self) -> list[KSP]:
        """
        Get the array of linear solver objects associated with the CISS solver.

        Not collective.

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
        """
        cdef PetscInt i = 0, n = 0
        cdef PetscKSP *p = NULL
        CHKERR( EPSCISSGetKSPs(self.eps, &n, &p) )
        return [ref_KSP(p[i]) for i from 0 <= i <n]

    #
    property problem_type:
        """The type of the eigenvalue problem."""
        def __get__(self) -> EPSProblemType:
            return self.getProblemType()
        def __set__(self, value):
            self.setProblemType(value)

    property extraction:
        """The type of extraction technique to be employed."""
        def __get__(self) -> EPSExtraction:
            return self.getExtraction()
        def __set__(self, value):
            self.setExtraction(value)

    property which:
        """The portion of the spectrum to be sought."""
        def __get__(self) -> EPSWhich:
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
        """The tolerance."""
        def __get__(self) -> float:
            return self.getTolerances()[0]
        def __set__(self, value):
            self.setTolerances(tol=value)

    property max_it:
        """The maximum iteration count."""
        def __get__(self) -> int:
            return self.getTolerances()[1]
        def __set__(self, value):
            self.setTolerances(max_it=value)

    property two_sided:
        """Two-sided that also computes left eigenvectors."""
        def __get__(self) -> bool:
            return self.getTwoSided()
        def __set__(self, value):
            self.setTwoSided(value)

    property true_residual:
        """Compute the true residual explicitly."""
        def __get__(self) -> bool:
            return self.getTrueResidual()
        def __set__(self, value):
            self.setTrueResidual(value)

    property purify:
        """Eigenvector purification."""
        def __get__(self) -> bool:
            return self.getPurify()
        def __set__(self, value):
            self.setPurify(value)

    property track_all:
        """Compute the residual norm of all approximate eigenpairs."""
        def __get__(self) -> bool:
            return self.getTrackAll()
        def __set__(self, value):
            self.setTrackAll(value)

    property st:
        """The spectral transformation (ST) object associated."""
        def __get__(self) -> ST:
            return self.getST()
        def __set__(self, value):
            self.setST(value)

    property bv:
        """The basis vectors (BV) object associated."""
        def __get__(self) -> BV:
            return self.getBV()
        def __set__(self, value):
            self.setBV(value)

    property rg:
        """The region (RG) object associated."""
        def __get__(self) -> RG:
            return self.getRG()
        def __set__(self, value):
            self.setRG(value)

    property ds:
        """The direct solver (DS) object associated."""
        def __get__(self) -> DS:
            return self.getDS()
        def __set__(self, value):
            self.setDS(value)

# -----------------------------------------------------------------------------

del EPSType
del EPSProblemType
del EPSExtraction
del EPSBalance
del EPSErrorType
del EPSWhich
del EPSConv
del EPSStop
del EPSConvergedReason
del EPSPowerShiftType
del EPSLanczosReorthogType
del EPSCISSQuadRule
del EPSCISSExtraction

# -----------------------------------------------------------------------------
