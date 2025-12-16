# -----------------------------------------------------------------------------

class PEPType(object):
    """
    PEP type.

    - `TOAR`:     Two-level orthogonal Arnoldi.
    - `STOAR`:    Symmetric TOAR.
    - `QARNOLDI`: Q-Arnoldi for quadratic problems.
    - `LINEAR`:   Linearization via EPS.
    - `JD`:       Polynomial Jacobi-Davidson.
    - `CISS`:     Contour integral spectrum slice.

    See Also
    --------
    slepc.PEPType
    """
    TOAR     = S_(PEPTOAR)
    STOAR    = S_(PEPSTOAR)
    QARNOLDI = S_(PEPQARNOLDI)
    LINEAR   = S_(PEPLINEAR)
    JD       = S_(PEPJD)
    CISS     = S_(PEPCISS)

class PEPProblemType(object):
    """
    PEP problem type.

    - `GENERAL`:    No structure.
    - `HERMITIAN`:  Hermitian structure.
    - `HYPERBOLIC`: QEP with Hermitian matrices, :math:`M>0`,
      :math:`(x^TCx)^2 > 4(x^TMx)(x^TKx)`.
    - `GYROSCOPIC`: QEP with :math:`M`, :math:`K`  Hermitian,
      :math:`M>0`, :math:`C` skew-Hermitian.

    See Also
    --------
    slepc.PEPProblemType
    """
    GENERAL    = PEP_GENERAL
    HERMITIAN  = PEP_HERMITIAN
    HYPERBOLIC = PEP_HYPERBOLIC
    GYROSCOPIC = PEP_GYROSCOPIC

class PEPWhich(object):
    """
    PEP desired part of spectrum.

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
    - `USER`:               User-defined criterion.

    See Also
    --------
    slepc.PEPWhich
    """
    LARGEST_MAGNITUDE  = PEP_LARGEST_MAGNITUDE
    SMALLEST_MAGNITUDE = PEP_SMALLEST_MAGNITUDE
    LARGEST_REAL       = PEP_LARGEST_REAL
    SMALLEST_REAL      = PEP_SMALLEST_REAL
    LARGEST_IMAGINARY  = PEP_LARGEST_IMAGINARY
    SMALLEST_IMAGINARY = PEP_SMALLEST_IMAGINARY
    TARGET_MAGNITUDE   = PEP_TARGET_MAGNITUDE
    TARGET_REAL        = PEP_TARGET_REAL
    TARGET_IMAGINARY   = PEP_TARGET_IMAGINARY
    ALL                = PEP_ALL
    USER               = PEP_WHICH_USER

class PEPBasis(object):
    """
    PEP basis type for the representation of the polynomial.

    - `MONOMIAL`:   Monomials (default).
    - `CHEBYSHEV1`: Chebyshev polynomials of the 1st kind.
    - `CHEBYSHEV2`: Chebyshev polynomials of the 2nd kind.
    - `LEGENDRE`:   Legendre polynomials.
    - `LAGUERRE`:   Laguerre polynomials.
    - `HERMITE`:    Hermite polynomials.

    See Also
    --------
    slepc.PEPBasis
    """
    MONOMIAL   = PEP_BASIS_MONOMIAL
    CHEBYSHEV1 = PEP_BASIS_CHEBYSHEV1
    CHEBYSHEV2 = PEP_BASIS_CHEBYSHEV2
    LEGENDRE   = PEP_BASIS_LEGENDRE
    LAGUERRE   = PEP_BASIS_LAGUERRE
    HERMITE    = PEP_BASIS_HERMITE

class PEPScale(object):
    """
    PEP scaling strategy.

    - `NONE`:     No scaling.
    - `SCALAR`:   Parameter scaling.
    - `DIAGONAL`: Diagonal scaling.
    - `BOTH`:     Both parameter and diagonal scaling.

    See Also
    --------
    slepc.PEPScale
    """
    NONE     = PEP_SCALE_NONE
    SCALAR   = PEP_SCALE_SCALAR
    DIAGONAL = PEP_SCALE_DIAGONAL
    BOTH     = PEP_SCALE_BOTH

class PEPRefine(object):
    """
    PEP refinement strategy.

    - `NONE`:     No refinement.
    - `SIMPLE`:   Refine eigenpairs one by one.
    - `MULTIPLE`: Refine all eigenpairs simultaneously (invariant pair).

    See Also
    --------
    slepc.PEPRefine
    """
    NONE     = PEP_REFINE_NONE
    SIMPLE   = PEP_REFINE_SIMPLE
    MULTIPLE = PEP_REFINE_MULTIPLE

class PEPRefineScheme(object):
    """
    PEP scheme for solving linear systems during iterative refinement.

    - `SCHUR`:    Schur complement.
    - `MBE`:      Mixed block elimination.
    - `EXPLICIT`: Build the explicit matrix.

    See Also
    --------
    slepc.PEPRefineScheme
    """
    SCHUR    = PEP_REFINE_SCHEME_SCHUR
    MBE      = PEP_REFINE_SCHEME_MBE
    EXPLICIT = PEP_REFINE_SCHEME_EXPLICIT

class PEPExtract(object):
    """
    PEP extraction strategy used.

    PEP extraction strategy used to obtain eigenvectors of the PEP from the
    eigenvectors of the linearization.

    - `NONE`:       Use the first block.
    - `NORM`:       Use the first or last block depending on norm of H.
    - `RESIDUAL`:   Use the block with smallest residual.
    - `STRUCTURED`: Combine all blocks in a certain way.

    See Also
    --------
    slepc.PEPExtract
    """
    NONE       = PEP_EXTRACT_NONE
    NORM       = PEP_EXTRACT_NORM
    RESIDUAL   = PEP_EXTRACT_RESIDUAL
    STRUCTURED = PEP_EXTRACT_STRUCTURED

class PEPErrorType(object):
    """
    PEP error type to assess accuracy of computed solutions.

    - `ABSOLUTE`: Absolute error.
    - `RELATIVE`: Relative error.
    - `BACKWARD`: Backward error.

    See Also
    --------
    slepc.PEPErrorType
    """
    ABSOLUTE = PEP_ERROR_ABSOLUTE
    RELATIVE = PEP_ERROR_RELATIVE
    BACKWARD = PEP_ERROR_BACKWARD

class PEPConv(object):
    """
    PEP convergence test.

    - `ABS`:  Absolute convergence test.
    - `REL`:  Convergence test relative to the eigenvalue.
    - `NORM`: Convergence test relative to the matrix norms.
    - `USER`: User-defined convergence test.

    See Also
    --------
    slepc.PEPConv
    """
    ABS  = PEP_CONV_ABS
    REL  = PEP_CONV_REL
    NORM = PEP_CONV_NORM
    USER = PEP_CONV_USER

class PEPStop(object):
    """
    PEP stopping test.

    - `BASIC`: Default stopping test.
    - `USER`:  User-defined stopping test.

    See Also
    --------
    slepc.PEPStop
    """
    BASIC = PEP_STOP_BASIC
    USER  = PEP_STOP_USER

class PEPConvergedReason(object):
    """
    PEP convergence reasons.

    - `CONVERGED_TOL`:          All eigenpairs converged to requested tolerance.
    - `CONVERGED_USER`:         User-defined convergence criterion satisfied.
    - `DIVERGED_ITS`:           Maximum number of iterations exceeded.
    - `DIVERGED_BREAKDOWN`:     Solver failed due to breakdown.
    - `DIVERGED_SYMMETRY_LOST`: Lanczos-type method could not preserve symmetry.
    - `CONVERGED_ITERATING`:    Iteration not finished yet.

    See Also
    --------
    slepc.PEPConvergedReason
    """
    CONVERGED_TOL          = PEP_CONVERGED_TOL
    CONVERGED_USER         = PEP_CONVERGED_USER
    DIVERGED_ITS           = PEP_DIVERGED_ITS
    DIVERGED_BREAKDOWN     = PEP_DIVERGED_BREAKDOWN
    DIVERGED_SYMMETRY_LOST = PEP_DIVERGED_SYMMETRY_LOST
    CONVERGED_ITERATING    = PEP_CONVERGED_ITERATING
    ITERATING              = PEP_CONVERGED_ITERATING

class PEPJDProjection(object):
    """
    PEP type of projection to be used in the Jacobi-Davidson solver.

    - `HARMONIC`:   Harmonic projection.
    - `ORTHOGONAL`: Orthogonal projection.

    See Also
    --------
    slepc.PEPJDProjection
    """
    HARMONIC   = PEP_JD_PROJECTION_HARMONIC
    ORTHOGONAL = PEP_JD_PROJECTION_ORTHOGONAL

class PEPCISSExtraction(object):
    """
    PEP CISS extraction technique.

    - `RITZ`:   Ritz extraction.
    - `HANKEL`: Extraction via Hankel eigenproblem.
    - `CAA`:    Communication-avoiding Arnoldi.

    See Also
    --------
    slepc.PEPCISSExtraction
    """
    RITZ   = PEP_CISS_EXTRACTION_RITZ
    HANKEL = PEP_CISS_EXTRACTION_HANKEL
    CAA    = PEP_CISS_EXTRACTION_CAA

# -----------------------------------------------------------------------------

cdef class PEP(Object):

    """
    Polynomial Eigenvalue Problem Solver.

    The Polynomial Eigenvalue Problem (`PEP`) solver is the object provided
    by slepc4py for specifying a polynomial eigenvalue problem. Apart from the
    specific solvers for this type of problems, there is an `EPS`-based solver,
    i.e., it uses a solver from `EPS` to solve a generalized eigenproblem
    obtained after linearization.
    """

    Type            = PEPType
    ProblemType     = PEPProblemType
    Which           = PEPWhich
    Basis           = PEPBasis
    Scale           = PEPScale
    Refine          = PEPRefine
    RefineScheme    = PEPRefineScheme
    Extract         = PEPExtract
    ErrorType       = PEPErrorType
    Conv            = PEPConv
    Stop            = PEPStop
    ConvergedReason = PEPConvergedReason

    JDProjection    = PEPJDProjection
    CISSExtraction  = PEPCISSExtraction

    def __cinit__(self):
        self.obj = <PetscObject*> &self.pep
        self.pep = NULL

    def view(self, Viewer viewer=None) -> None:
        """
        Print the PEP data structure.

        Collective.

        Parameters
        ----------
        viewer
            Visualization context; if not provided, the standard
            output is used.

        See Also
        --------
        slepc.PEPView
        """
        cdef PetscViewer vwr = def_Viewer(viewer)
        CHKERR( PEPView(self.pep, vwr) )

    def destroy(self) -> Self:
        """
        Destroy the PEP object.

        Collective.

        See Also
        --------
        slepc.PEPDestroy
        """
        CHKERR( PEPDestroy(&self.pep) )
        self.pep = NULL
        return self

    def reset(self) -> None:
        """
        Reset the PEP object.

        Collective.

        See Also
        --------
        slepc.PEPReset
        """
        CHKERR( PEPReset(self.pep) )

    def create(self, comm: Comm | None = None) -> Self:
        """
        Create the PEP object.

        Collective.

        Parameters
        ----------
        comm
            MPI communicator. If not provided, it defaults to all processes.

        See Also
        --------
        slepc.PEPCreate
        """
        cdef MPI_Comm ccomm = def_Comm(comm, SLEPC_COMM_DEFAULT())
        cdef SlepcPEP newpep = NULL
        CHKERR( PEPCreate(ccomm, &newpep) )
        CHKERR( SlepcCLEAR(self.obj) ); self.pep = newpep
        return self

    def setType(self, pep_type: Type | str) -> None:
        """
        Set the particular solver to be used in the PEP object.

        Logically collective.

        Parameters
        ----------
        pep_type
            The solver to be used.

        Notes
        -----
        The default is `TOAR`. Normally, it is best to use
        `setFromOptions()` and then set the PEP type from the options
        database rather than by using this routine. Using the options
        database provides the user with maximum flexibility in
        evaluating the different available methods.

        See Also
        --------
        getType, slepc.PEPSetType
        """
        cdef SlepcPEPType cval = NULL
        pep_type = str2bytes(pep_type, &cval)
        CHKERR( PEPSetType(self.pep, cval) )

    def getType(self) -> str:
        """
        Get the PEP type of this object.

        Not collective.

        Returns
        -------
        str
            The solver currently being used.

        See Also
        --------
        setType, slepc.PEPGetType
        """
        cdef SlepcPEPType pep_type = NULL
        CHKERR( PEPGetType(self.pep, &pep_type) )
        return bytes2str(pep_type)

    def getOptionsPrefix(self) -> str:
        """
        Get the prefix used for searching for all PEP options in the database.

        Not collective.

        Returns
        -------
        str
            The prefix string set for this PEP object.

        See Also
        --------
        setOptionsPrefix, appendOptionsPrefix, slepc.PEPGetOptionsPrefix
        """
        cdef const char *prefix = NULL
        CHKERR( PEPGetOptionsPrefix(self.pep, &prefix) )
        return bytes2str(prefix)

    def setOptionsPrefix(self, prefix: str | None = None) -> None:
        """
        Set the prefix used for searching for all PEP options in the database.

        Logically collective.

        Parameters
        ----------
        prefix
            The prefix string to prepend to all PEP option requests.

        Notes
        -----
        A hyphen (-) must NOT be given at the beginning of the prefix
        name.  The first character of all runtime options is
        AUTOMATICALLY the hyphen.

        For example, to distinguish between the runtime options for
        two different PEP contexts, one could call::

            P1.setOptionsPrefix("pep1_")
            P2.setOptionsPrefix("pep2_")

        See Also
        --------
        appendOptionsPrefix, getOptionsPrefix, slepc.PEPGetOptionsPrefix
        """
        cdef const char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( PEPSetOptionsPrefix(self.pep, cval) )

    def appendOptionsPrefix(self, prefix: str | None = None) -> None:
        """
        Append to the prefix used for searching for all PEP options in the database.

        Logically collective.

        Parameters
        ----------
        prefix
            The prefix string to prepend to all PEP option requests.

        See Also
        --------
        setOptionsPrefix, getOptionsPrefix, slepc.PEPAppendOptionsPrefix
        """
        cdef const char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( PEPAppendOptionsPrefix(self.pep, cval) )

    def setFromOptions(self) -> None:
        """
        Set PEP options from the options database.

        Collective.

        Notes
        -----
        To see all options, run your program with the ``-help`` option.

        This routine must be called before `setUp()` if the user is to be
        allowed to set the solver type.

        See Also
        --------
        setOptionsPrefix, slepc.PEPSetFromOptions
        """
        CHKERR( PEPSetFromOptions(self.pep) )

    def getBasis(self) -> Basis:
        """
        Get the type of polynomial basis used.

        Not collective.

        Returns
        -------
        Basis
            The basis that was previously set.

        See Also
        --------
        setBasis, slepc.PEPGetBasis
        """
        cdef SlepcPEPBasis val = PEP_BASIS_MONOMIAL
        CHKERR( PEPGetBasis(self.pep, &val) )
        return val

    def setBasis(self, basis: Basis) -> None:
        """
        Set the type of polynomial basis used.

        Logically collective.

        Set the type of polynomial basis used to describe the polynomial
        eigenvalue problem.

        Parameters
        ----------
        basis
            The basis to be set.

        Notes
        -----
        By default, the coefficient matrices passed via `setOperators()` are
        expressed in the monomial basis, i.e.
        :math:`P(\lambda)=A_0+\lambda A_1+\lambda^2 A_2+\dots+\lambda^d A_d`.
        Other polynomial bases may have better numerical behavior, but the user
        must then pass the coefficient matrices accordingly.

        See Also
        --------
        getBasis, setOperators, slepc.PEPSetBasis
        """
        cdef SlepcPEPBasis val = basis
        CHKERR( PEPSetBasis(self.pep, val) )

    def getProblemType(self) -> ProblemType:
        """
        Get the problem type from the PEP object.

        Not collective.

        Returns
        -------
        ProblemType
            The problem type that was previously set.

        See Also
        --------
        setProblemType, slepc.PEPGetProblemType
        """
        cdef SlepcPEPProblemType val = PEP_GENERAL
        CHKERR( PEPGetProblemType(self.pep, &val) )
        return val

    def setProblemType(self, problem_type: ProblemType) -> None:
        """
        Set the type of the polynomial eigenvalue problem.

        Logically collective.

        Parameters
        ----------
        problem_type
            The problem type to be set.

        Notes
        -----
        This function is used to instruct SLEPc to exploit certain
        structure in the polynomial eigenproblem. By default, no
        particular structure is assumed.

        If the problem matrices are Hermitian (symmetric in the real
        case) or Hermitian/skew-Hermitian then the solver can exploit
        this fact to perform less operations or provide better stability.
        Hyperbolic problems are a particular case of Hermitian problems,
        some solvers may treat them simply as Hermitian.

        See Also
        --------
        setOperators, setType, getProblemType, slepc.PEPSetProblemType
        """
        cdef SlepcPEPProblemType val = problem_type
        CHKERR( PEPSetProblemType(self.pep, val) )

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
        setWhichEigenpairs, slepc.PEPGetWhichEigenpairs
        """
        cdef SlepcPEPWhich val = PEP_LARGEST_MAGNITUDE
        CHKERR( PEPGetWhichEigenpairs(self.pep, &val) )
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
        Not all eigensolvers implemented in PEP account for all the
        possible values. Also, some values make sense only for certain
        types of problems. If SLEPc is compiled for real numbers
        `PEP.Which.LARGEST_IMAGINARY` and
        `PEP.Which.SMALLEST_IMAGINARY` use the absolute value of the
        imaginary part for eigenvalue selection.

        The target is a scalar value provided with `setTarget()`.

        The criterion `PEP.Which.TARGET_IMAGINARY` is available only
        in case PETSc and SLEPc have been built with complex scalars.

        `PEP.Which.ALL` is intended for use in combination with an
        interval (see `setInterval()`), when all eigenvalues within the
        interval are requested, or in the context of the `PEP.Type.CISS`
        solver for computing all eigenvalues in a region.

        See Also
        --------
        getWhichEigenpairs, setTarget, setInterval, slepc.PEPSetWhichEigenpairs
        """
        cdef SlepcPEPWhich val = which
        CHKERR( PEPSetWhichEigenpairs(self.pep, val) )

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
        setTarget, slepc.PEPGetTarget
        """
        cdef PetscScalar sval = 0
        CHKERR( PEPGetTarget(self.pep, &sval) )
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
        getTarget, setWhichEigenpairs, slepc.PEPSetTarget
        """
        cdef PetscScalar sval = asScalar(target)
        CHKERR( PEPSetTarget(self.pep, sval) )

    def getTolerances(self) -> tuple[float, int]:
        """
        Get the tolerance and maximum iteration count.

        Not collective.

        Get the tolerance and maximum iteration count used by the default PEP
        convergence tests.

        Returns
        -------
        tol: float
            The convergence tolerance.
        max_it: int
            The maximum number of iterations.

        See Also
        --------
        setTolerances, slepc.PEPGetTolerances
        """
        cdef PetscReal rval = 0
        cdef PetscInt  ival = 0
        CHKERR( PEPGetTolerances(self.pep, &rval, &ival) )
        return (toReal(rval), toInt(ival))

    def setTolerances(self, tol: float | None = None, max_it: int | None = None) -> None:
        """
        Set the tolerance and maximum iteration count.

        Logically collective.

        Set the tolerance and maximum iteration count used by the default PEP
        convergence tests.

        Parameters
        ----------
        tol
            The convergence tolerance.
        max_it
            The maximum number of iterations

        Notes
        -----
        Use `DETERMINE` for ``max_it`` to assign a reasonably good value,
        which is dependent on the solution method.

        See Also
        --------
        getTolerances, slepc.PEPSetTolerances
        """
        cdef PetscReal rval = PETSC_CURRENT
        cdef PetscInt  ival = PETSC_CURRENT
        if tol    is not None: rval = asReal(tol)
        if max_it is not None: ival = asInt(max_it)
        CHKERR( PEPSetTolerances(self.pep, rval, ival) )

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

        See Also
        --------
        setInterval, slepc.PEPGetInterval
        """
        cdef PetscReal inta = 0
        cdef PetscReal intb = 0
        CHKERR( PEPGetInterval(self.pep, &inta, &intb) )
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
        eigenvalues of symmetric quadratic eigenproblems in a given interval.
        This function provides the interval to be considered. It must
        be used in combination with `PEP.Which.ALL`, see
        `setWhichEigenpairs()`. Note that in polynomial eigenproblems
        spectrum slicing is implemented in `STOAR` only.

        See Also
        --------
        getInterval, slepc.PEPSetInterval
        """
        cdef PetscReal rval1 = asReal(inta)
        cdef PetscReal rval2 = asReal(intb)
        CHKERR( PEPSetInterval(self.pep, rval1, rval2) )

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
        setConvergenceTest, slepc.PEPGetConvergenceTest
        """
        cdef SlepcPEPConv conv = PEP_CONV_REL
        CHKERR( PEPGetConvergenceTest(self.pep, &conv) )
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
        getConvergenceTest, slepc.PEPSetConvergenceTest
        """
        cdef SlepcPEPConv tconv = conv
        CHKERR( PEPSetConvergenceTest(self.pep, tconv) )

    def getRefine(self) -> tuple[Refine, int, float, int, RefineScheme]:
        """
        Get the refinement strategy used by the PEP object.

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
        setRefine, slepc.PEPGetRefine
        """
        cdef SlepcPEPRefine ref = PEP_REFINE_NONE
        cdef PetscInt npart = 1
        cdef PetscReal tol = PETSC_DEFAULT
        cdef PetscInt its = PETSC_DEFAULT
        cdef SlepcPEPRefineScheme scheme = PEP_REFINE_SCHEME_MBE
        CHKERR( PEPGetRefine(self.pep, &ref, &npart, &tol, &its, &scheme) )
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
        Set the refinement strategy used by the PEP object.

        Logically collective.

        Set the refinement strategy used by the PEP object, and the associated
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
        getRefine, slepc.PEPSetRefine
        """
        cdef SlepcPEPRefine tref = ref
        cdef PetscInt tnpart = PETSC_CURRENT
        cdef PetscReal ttol = PETSC_CURRENT
        cdef PetscInt tits = PETSC_CURRENT
        cdef SlepcPEPRefineScheme tscheme = PEP_REFINE_SCHEME_MBE
        if npart is not None: tnpart = asInt(npart)
        if tol is not None: ttol = asReal(tol)
        if its is not None: tits = asInt(its)
        if scheme is not None: tscheme = scheme
        CHKERR( PEPSetRefine(self.pep, tref, tnpart, ttol, tits, tscheme) )

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
        setRefine, slepc.PEPRefineGetKSP
        """
        cdef KSP ksp = KSP()
        CHKERR( PEPRefineGetKSP(self.pep, &ksp.ksp) )
        CHKERR( PetscINCREF(ksp.obj) )
        return ksp

    def setExtract(self, extract: Extract) -> None:
        """
        Set the extraction strategy to be used.

        Logically collective.

        Parameters
        ----------
        extract
            The extraction strategy.

        Notes
        -----
        This is relevant for solvers based on linearization. Once the
        solver has converged, the polynomial eigenvectors can be
        extracted from the eigenvectors of the linearized problem in
        different ways.

        See Also
        --------
        getExtract, slepc.PEPSetExtract
        """
        cdef SlepcPEPExtract val = extract
        CHKERR( PEPSetExtract(self.pep, val) )

    def getExtract(self) -> Extract:
        """
        Get the extraction technique used by the `PEP` object.

        Not collective.

        Returns
        -------
        Extract
            The extraction strategy.

        See Also
        --------
        setExtract, slepc.PEPGetExtract
        """
        cdef SlepcPEPExtract val = PEP_EXTRACT_NONE
        CHKERR( PEPGetExtract(self.pep, &val) )
        return val

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
        setTrackAll, slepc.PEPGetTrackAll
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( PEPGetTrackAll(self.pep, &tval) )
        return toBool(tval)

    def setTrackAll(self, trackall: bool) -> None:
        """
        Set flag to compute the residual of all approximate eigenpairs.

        Logically collective.

        Set if the solver must compute the residual of all approximate
        eigenpairs or not.

        Parameters
        ----------
        trackall
            Whether to compute all residuals or not.

        See Also
        --------
        getTrackAll, slepc.PEPSetTrackAll
        """
        cdef PetscBool tval = trackall
        CHKERR( PEPSetTrackAll(self.pep, tval) )

    def getDimensions(self) -> tuple[int, int, int]:
        """
        Get the number of eigenvalues to compute and the dimension of the subspace.

        Not collective.

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
        setDimensions, slepc.PEPGetDimensions
        """
        cdef PetscInt ival1 = 0
        cdef PetscInt ival2 = 0
        cdef PetscInt ival3 = 0
        CHKERR( PEPGetDimensions(self.pep, &ival1, &ival2, &ival3) )
        return (toInt(ival1), toInt(ival2), toInt(ival3))

    def setDimensions(
        self,
        nev: int | None = None,
        ncv: int | None = None,
        mpd: int | None = None,
    ) -> None:
        """
        Set the number of eigenvalues to compute and the dimension of the subspace.

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
        getDimensions, slepc.PEPSetDimensions
        """
        cdef PetscInt ival1 = PETSC_CURRENT
        cdef PetscInt ival2 = PETSC_CURRENT
        cdef PetscInt ival3 = PETSC_CURRENT
        if nev is not None: ival1 = asInt(nev)
        if ncv is not None: ival2 = asInt(ncv)
        if mpd is not None: ival3 = asInt(mpd)
        CHKERR( PEPSetDimensions(self.pep, ival1, ival2, ival3) )

    def getST(self) -> ST:
        """
        Get the spectral transformation object associated to the eigensolver.

        Not collective.

        Returns
        -------
        ST
            The spectral transformation.

        See Also
        --------
        setST, slepc.PEPGetST
        """
        cdef ST st = ST()
        CHKERR( PEPGetST(self.pep, &st.st) )
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

        See Also
        --------
        getST, slepc.PEPSetST
        """
        CHKERR( PEPSetST(self.pep, st.st) )

    def getScale(
        self,
        Vec Dl: petsc4py.PETSc.Vec | None = None,
        Vec Dr: petsc4py.PETSc.Vec | None = None
    ) -> tuple[Scale, float, int, float]:
        """
        Get the strategy used for scaling the polynomial eigenproblem.

        Not collective.

        Parameters
        ----------
        Dl
            Placeholder for the returned left diagonal matrix.
        Dr
            Placeholder for the returned right diagonal matrix.

        Returns
        -------
        scale: Scale
            The scaling strategy.
        alpha: float
            The scaling factor.
        its: int
            The number of iterations of diagonal scaling.
        lbda: float
            Approximation of the wanted eigenvalues (modulus).

        See Also
        --------
        setScale, slepc.PEPGetScale
        """
        cdef SlepcPEPScale scale = PEP_SCALE_NONE
        cdef PetscReal alpha = 0
        cdef PetscInt its = 0
        cdef PetscReal lbda = 0
        cdef PetscVec vecl = NULL
        cdef PetscVec vecr = NULL
        CHKERR( PEPGetScale(self.pep, &scale, &alpha, &vecl, &vecr, &its, &lbda) )
        if Dl.vec != NULL:
            if vecl != NULL:
                CHKERR( VecCopy(vecl, Dl.vec) )
            else:
                CHKERR( VecSet(Dl.vec, 1.0) )
        if Dr.vec != NULL:
            if vecr != NULL:
                CHKERR( VecCopy(vecr, Dr.vec) )
            else:
                CHKERR( VecSet(Dr.vec, 1.0) )
        CHKERR( VecDestroy(&vecl) )
        CHKERR( VecDestroy(&vecr) )
        return (scale, toReal(alpha), toInt(its), toReal(lbda))

    def setScale(
        self,
        scale: Scale,
        alpha: float | None = None,
        Vec Dl: petsc4py.PETSc.Vec | None = None,
        Vec Dr: petsc4py.PETSc.Vec | None = None,
        its: int | None = None,
        lbda: float | None = None,
    ) -> None:
        """
        Set the scaling strategy to be used.

        Collective.

        Parameters
        ----------
        scale
            The scaling strategy.
        alpha
            The scaling factor.
        Dl
            The left diagonal matrix.
        Dr
            The right diagonal matrix.
        its
            The number of iterations of diagonal scaling.
        lbda
            Approximation of the wanted eigenvalues (modulus).

        See Also
        --------
        getScale, slepc.PEPSetScale
        """
        cdef SlepcPEPScale senum = scale
        cdef PetscReal rval1 = PETSC_CURRENT
        cdef PetscInt ival = PETSC_CURRENT
        cdef PetscReal rval2 = PETSC_CURRENT
        cdef PetscVec vecl = NULL
        cdef PetscVec vecr = NULL
        if alpha is not None: rval1 = asReal(alpha)
        if Dl is not None:    vecl = Dl.vec
        if Dr is not None:    vecr = Dr.vec
        if its is not None:   ival = asInt(its)
        if lbda is not None:  rval2 = asReal(lbda)
        CHKERR( PEPSetScale(self.pep, senum, rval1, vecl, vecr, ival, rval2) )

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
        setBV, slepc.PEPGetBV
        """
        cdef BV bv = BV()
        CHKERR( PEPGetBV(self.pep, &bv.bv) )
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

        See Also
        --------
        getBV, slepc.PEPSetBV
        """
        CHKERR( PEPSetBV(self.pep, bv.bv) )

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
        setRG, slepc.PEPGetRG
        """
        cdef RG rg = RG()
        CHKERR( PEPGetRG(self.pep, &rg.rg) )
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
        getRG, slepc.PEPSetRG
        """
        CHKERR( PEPSetRG(self.pep, rg.rg) )

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
        setDS, slepc.PEPGetDS
        """
        cdef DS ds = DS()
        CHKERR( PEPGetDS(self.pep, &ds.ds) )
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
        getDS, slepc.PEPSetDS
        """
        CHKERR( PEPSetDS(self.pep, ds.ds) )

    def getOperators(self) -> list[Mat]:
        """
        Get the matrices associated with the eigenvalue problem.

        Collective.

        Returns
        -------
        list of petsc4py.PETSc.Mat
            The matrices associated with the eigensystem.

        See Also
        --------
        setOperators, slepc.PEPGetOperators
        """
        cdef Mat A
        cdef PetscMat mat = NULL
        cdef PetscInt k=0, n=0
        CHKERR( PEPGetNumMatrices(self.pep, &n) )
        cdef object operators = []
        for k from 0 <= k < n:
            CHKERR( PEPGetOperators(self.pep, k, &mat) )
            A = Mat(); A.mat = mat; CHKERR( PetscINCREF(A.obj) )
            operators.append(A)
        return tuple(operators)

    def setOperators(self, operators: list[Mat]) -> None:
        """
        Set the matrices associated with the eigenvalue problem.

        Collective.

        Parameters
        ----------
        operators
            The matrices associated with the eigensystem.

        Notes
        -----
        The polynomial eigenproblem is defined as :math:`P(\lambda)x=0`,
        where :math:`\lambda` is the eigenvalue, :math:`x` is the eigenvector,
        and :math:`P` is defined as
        :math:`P(\lambda) = A_0 + \lambda A_1 + \dots + \lambda^d A_d`, with
        :math:`d` = ``nmat``-1 (the degree of :math:`P`). For non-monomial
        bases, this expression is different.

        See Also
        --------
        getOperators, slepc.PEPSetOperators
        """
        operators = tuple(operators)
        cdef PetscMat *mats = NULL
        cdef Py_ssize_t k=0, n = len(operators)
        cdef tmp = allocate(<size_t>n*sizeof(PetscMat),<void**>&mats)
        for k from 0 <= k < n: mats[k] = (<Mat?>operators[k]).mat
        CHKERR( PEPSetOperators(self.pep, <PetscInt>n, mats) )

    #

    def setInitialSpace(self, space: Vec | list[Vec]) -> None:
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
        setUp, slepc.PEPSetInitialSpace
        """
        if isinstance(space, Vec): space = [space]
        cdef PetscVec *vs = NULL
        cdef Py_ssize_t i = 0, ns = len(space)
        cdef tmp = allocate(<size_t>ns*sizeof(PetscVec),<void**>&vs)
        for i in range(ns): vs[i] = (<Vec?>space[i]).vec
        CHKERR( PEPSetInitialSpace(self.pep, <PetscInt>ns, vs) )

    #

    def setStoppingTest(
        self,
        stopping: PEPStoppingFunction | None,
        args: tuple[Any, ...] | None = None,
        kargs: dict[str, Any] | None = None,
    ) -> None:
        """
        Set a function to decide when to stop the outer iteration of the eigensolver.

        Logically collective.

        See Also
        --------
        getStoppingTest, slepc.PEPSetStoppingTestFunction
        """
        if stopping is not None:
            if args is None: args = ()
            if kargs is None: kargs = {}
            self.set_attr('__stopping__', (stopping, args, kargs))
            CHKERR( PEPSetStoppingTestFunction(self.pep, PEP_Stopping, NULL, NULL) )
        else:
            self.set_attr('__stopping__', None)
            CHKERR( PEPSetStoppingTestFunction(self.pep, PEPStoppingBasic, NULL, NULL) )

    def getStoppingTest(self) -> PEPStoppingFunction:
        """
        Get the stopping test function.

        Not collective.

        Returns
        -------
        PEPStoppingFunction
            The stopping test function.

        See Also
        --------
        setStoppingTest
        """
        return self.get_attr('__stopping__')

    def setEigenvalueComparison(
        self,
        comparison: PEPEigenvalueComparison | None,
        args: tuple[Any, ...] | None = None,
        kargs: dict[str, Any] | None = None,
    ) -> None:
        """
        Set an eigenvalue comparison function.

        Logically collective.

        Notes
        -----
        This eigenvalue comparison function is used when `setWhichEigenpairs()`
        is set to `PEP.Which.USER`.

        See Also
        --------
        getEigenvalueComparison, slepc.PEPSetEigenvalueComparison
        """
        if comparison is not None:
            if args is None: args = ()
            if kargs is None: kargs = {}
            self.set_attr('__comparison__', (comparison, args, kargs))
            ctx = self.get_attr('__comparison__')
            CHKERR( PEPSetEigenvalueComparison(self.pep, PEP_Comparison, <void*>ctx) )
        else:
            self.set_attr('__comparison__', None)
            CHKERR( PEPSetEigenvalueComparison(self.pep, NULL, NULL) )

    def getEigenvalueComparison(self) -> PEPEigenvalueComparison:
        """
        Get the eigenvalue comparison function.

        Not collective.

        Returns
        -------
        PEPEigenvalueComparison
            The eigenvalue comparison function.

        See Also
        --------
        setEigenvalueComparison
        """
        return self.get_attr('__comparison__')

    def setMonitor(
        self,
        monitor: PEPMonitorFunction | None,
        args: tuple[Any, ...] | None = None,
        kargs: dict[str, Any] | None = None,
    ) -> None:
        """
        Append a monitor function to the list of monitors.

        Logically collective.

        See Also
        --------
        getMonitor, cancelMonitor, slepc.PEPMonitorSet
        """
        if monitor is None: return
        cdef object monitorlist = self.get_attr('__monitor__')
        if monitorlist is None:
            monitorlist = []
            self.set_attr('__monitor__', monitorlist)
            CHKERR( PEPMonitorSet(self.pep, PEP_Monitor, NULL, NULL) )
        if args is None: args = ()
        if kargs is None: kargs = {}
        monitorlist.append((monitor, args, kargs))

    def getMonitor(self) -> PEPMonitorFunction:
        """
        Get the list of monitor functions.

        Not collective.

        Returns
        -------
        PEPMonitorFunction
            The list of monitor functions.

        See Also
        --------
        setMonitor
        """
        return self.get_attr('__monitor__')

    def cancelMonitor(self) -> None:
        """
        Clear all monitors for a `PEP` object.

        Logically collective.

        See Also
        --------
        slepc.PEPMonitorCancel
        """
        CHKERR( PEPMonitorCancel(self.pep) )
        self.set_attr('__monitor__', None)

    #

    def setUp(self) -> None:
        """
        Set up all the internal data structures.

        Collective.

        Notes
        -----
        Sets up all the internal data structures necessary for the execution
        of the eigensolver. This includes the setup of the internal `ST`
        object.

        This function need not be called explicitly in most cases,
        since `solve()` calls it. It can be useful when one wants to
        measure the set-up time separately from the solve time.

        See Also
        --------
        solve, slepc.PEPSetUp
        """
        CHKERR( PEPSetUp(self.pep) )

    def solve(self) -> None:
        """
        Solve the polynomial eigenproblem.

        Collective.

        Notes
        -----
        The problem matrices are specified with `setOperators()`.

        `solve()` will return without generating an error regardless of
        whether all requested solutions were computed or not. Call
        `getConverged()` to get the actual number of computed solutions,
        and `getConvergedReason()` to determine if the solver converged
        or failed and why.

        See Also
        --------
        setUp, setOperators, getConverged, getConvergedReason, slepc.PEPSolve
        """
        CHKERR( PEPSolve(self.pep) )

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
        getConvergedReason, setTolerances, slepc.PEPGetIterationNumber
        """
        cdef PetscInt ival = 0
        CHKERR( PEPGetIterationNumber(self.pep, &ival) )
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
        setTolerances, solve, slepc.PEPGetConvergedReason
        """
        cdef SlepcPEPConvergedReason val = PEP_CONVERGED_ITERATING
        CHKERR( PEPGetConvergedReason(self.pep, &val) )
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
        setDimensions, solve, getEigenpair, slepc.PEPGetConverged
        """
        cdef PetscInt ival = 0
        CHKERR( PEPGetConverged(self.pep, &ival) )
        return toInt(ival)

    def getEigenpair(self, i: int, Vec Vr = None, Vec Vi = None) -> complex:
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
        solve, getConverged, setWhichEigenpairs, slepc.PEPGetEigenpair
        """
        cdef PetscScalar sval1 = 0
        cdef PetscScalar sval2 = 0
        cdef PetscVec vecr = Vr.vec if Vr is not None else <PetscVec>NULL
        cdef PetscVec veci = Vi.vec if Vi is not None else <PetscVec>NULL
        CHKERR( PEPGetEigenpair(self.pep, i, &sval1, &sval2, vecr, veci) )
        return toComplex(sval1, sval2)

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
        computeError, slepc.PEPGetErrorEstimate
        """
        cdef PetscReal rval = 0
        CHKERR( PEPGetErrorEstimate(self.pep, i, &rval) )
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
            :math:`\|P(\lambda)x\|_2` where :math:`\lambda` is the eigenvalue
            and :math:`x` is the eigenvector.

        Notes
        -----
        The index ``i`` should be a value between ``0`` and ``nconv-1``
        (see `getConverged()`).

        See Also
        --------
        getErrorEstimate, slepc.PEPComputeError
        """
        cdef SlepcPEPErrorType et = PEP_ERROR_BACKWARD
        cdef PetscReal rval = 0
        if etype is not None: et = etype
        CHKERR( PEPComputeError(self.pep, i, et, &rval) )
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
        solve, valuesView, vectorsView, slepc.PEPErrorView
        """
        cdef SlepcPEPErrorType et = PEP_ERROR_RELATIVE
        if etype is not None: et = etype
        cdef PetscViewer vwr = def_Viewer(viewer)
        CHKERR( PEPErrorView(self.pep, et, vwr) )

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
        solve, vectorsView, errorView, slepc.PEPValuesView
        """
        cdef PetscViewer vwr = def_Viewer(viewer)
        CHKERR( PEPValuesView(self.pep, vwr) )

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
        solve, valuesView, errorView, slepc.PEPVectorsView
        """
        cdef PetscViewer vwr = def_Viewer(viewer)
        CHKERR( PEPVectorsView(self.pep, vwr) )

    #

    def setLinearEPS(self, EPS eps) -> None:
        """
        Set an eigensolver object associated to the polynomial eigenvalue solver.

        Collective.

        Parameters
        ----------
        eps
            The linear eigensolver.

        See Also
        --------
        getLinearEPS, slepc.PEPLinearSetEPS
        """
        CHKERR( PEPLinearSetEPS(self.pep, eps.eps) )

    def getLinearEPS(self) -> EPS:
        """
        Get the eigensolver object associated to the polynomial eigenvalue solver.

        Collective.

        Returns
        -------
        EPS
            The linear eigensolver.

        See Also
        --------
        setLinearEPS, slepc.PEPLinearGetEPS
        """
        cdef EPS eps = EPS()
        CHKERR( PEPLinearGetEPS(self.pep, &eps.eps) )
        CHKERR( PetscINCREF(eps.obj) )
        return eps

    def setLinearLinearization(self, alpha: float = 1.0, beta: float = 0.0) -> None:
        """
        Set the coefficients that define the linearization of a quadratic eigenproblem.

        Logically collective.

        Parameters
        ----------
        alpha
            First parameter of the linearization.
        beta
            Second parameter of the linearization.

        See Also
        --------
        getLinearLinearization, slepc.PEPLinearSetLinearization
        """
        cdef PetscReal a = asReal(alpha)
        cdef PetscReal b = asReal(beta)
        CHKERR( PEPLinearSetLinearization(self.pep, a, b) )

    def getLinearLinearization(self) -> tuple[float, float]:
        """
        Get the coeffs. defining the linearization of a quadratic eigenproblem.

        Not collective.

        Returns
        -------
        alpha: float
            First parameter of the linearization.
        beta: float
            Second parameter of the linearization.

        See Also
        --------
        setLinearLinearization, slepc.PEPLinearGetLinearization
        """
        cdef PetscReal a = 0.0
        cdef PetscReal b = 0.0
        CHKERR( PEPLinearGetLinearization(self.pep, &a, &b) )
        return (asReal(a), asReal(b))

    def setLinearExplicitMatrix(self, flag: bool) -> None:
        """
        Set flag to explicitly build the matrices for the linearization.

        Logically collective.

        Parameters
        ----------
        flag
            Boolean flag indicating if the matrices are built explicitly.

        See Also
        --------
        getLinearExplicitMatrix, slepc.PEPLinearSetExplicitMatrix
        """
        cdef PetscBool sval = asBool(flag)
        CHKERR( PEPLinearSetExplicitMatrix(self.pep, sval) )

    def getLinearExplicitMatrix(self) -> bool:
        """
        Get if the matrices for the linearization are built explicitly.

        Not collective.

        Returns
        -------
        bool
            Boolean flag indicating if the matrices are built explicitly.

        See Also
        --------
        getLinearExplicitMatrix, slepc.PEPLinearSetExplicitMatrix
        """
        cdef PetscBool sval = PETSC_FALSE
        CHKERR( PEPLinearGetExplicitMatrix(self.pep, &sval) )
        return toBool(sval)

    #

    def setQArnoldiRestart(self, keep: float) -> None:
        """
        Set the restart parameter for the Q-Arnoldi method.

        Logically collective.

        Set the restart parameter for the Q-Arnoldi method, in
        particular the proportion of basis vectors that must be kept
        after restart.

        Parameters
        ----------
        keep
            The number of vectors to be kept at restart.

        Notes
        -----
        Allowed values are in the range [0.1,0.9]. The default is 0.5.

        See Also
        --------
        getQArnoldiRestart, slepc.PEPQArnoldiSetRestart
        """
        cdef PetscReal val = asReal(keep)
        CHKERR( PEPQArnoldiSetRestart(self.pep, val) )

    def getQArnoldiRestart(self) -> float:
        """
        Get the restart parameter used in the Q-Arnoldi method.

        Not collective.

        Returns
        -------
        float
            The number of vectors to be kept at restart.

        See Also
        --------
        setQArnoldiRestart, slepc.PEPQArnoldiGetRestart
        """
        cdef PetscReal val = 0
        CHKERR( PEPQArnoldiGetRestart(self.pep, &val) )
        return toReal(val)

    def setQArnoldiLocking(self, lock: bool) -> None:
        """
        Toggle between locking and non-locking variants of the Q-Arnoldi method.

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
        getQArnoldiLocking, slepc.PEPQArnoldiSetLocking
        """
        cdef PetscBool val = asBool(lock)
        CHKERR( PEPQArnoldiSetLocking(self.pep, val) )

    def getQArnoldiLocking(self) -> bool:
        """
        Get the locking flag used in the Q-Arnoldi method.

        Not collective.

        Returns
        -------
        bool
            The locking flag.

        See Also
        --------
        setQArnoldiLocking, slepc.PEPQArnoldiGetLocking
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( PEPQArnoldiGetLocking(self.pep, &tval) )
        return toBool(tval)

    #

    def setTOARRestart(self, keep: float) -> None:
        """
        Set the restart parameter for the TOAR method.

        Logically collective.

        Set the restart parameter for the TOAR method, in
        particular the proportion of basis vectors that must be kept
        after restart.

        Parameters
        ----------
        keep
            The number of vectors to be kept at restart.

        Notes
        -----
        Allowed values are in the range [0.1,0.9]. The default is 0.5.

        See Also
        --------
        getTOARRestart, slepc.PEPTOARSetRestart
        """
        cdef PetscReal val = asReal(keep)
        CHKERR( PEPTOARSetRestart(self.pep, val) )

    def getTOARRestart(self) -> float:
        """
        Get the restart parameter used in the TOAR method.

        Not collective.

        Returns
        -------
        float
            The number of vectors to be kept at restart.

        See Also
        --------
        setTOARRestart, slepc.PEPTOARGetRestart
        """
        cdef PetscReal val = 0
        CHKERR( PEPTOARGetRestart(self.pep, &val) )
        return toReal(val)

    def setTOARLocking(self, lock: bool) -> None:
        """
        Toggle between locking and non-locking variants of the TOAR method.

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
        getTOARLocking, slepc.PEPTOARSetLocking
        """
        cdef PetscBool val = asBool(lock)
        CHKERR( PEPTOARSetLocking(self.pep, val) )

    def getTOARLocking(self) -> bool:
        """
        Get the locking flag used in the TOAR method.

        Not collective.

        Returns
        -------
        bool
            The locking flag.

        See Also
        --------
        setTOARLocking, slepc.PEPTOARGetLocking
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( PEPTOARGetLocking(self.pep, &tval) )
        return toBool(tval)

    #

    def setSTOARLinearization(self, alpha: float = 1.0, beta: float = 0.0) -> None:
        """
        Set the coefficients that define the linearization of a quadratic eigenproblem.

        Logically collective.

        Parameters
        ----------
        alpha
            First parameter of the linearization.
        beta
            Second parameter of the linearization.

        See Also
        --------
        getSTOARLinearization, slepc.PEPSTOARSetLinearization
        """
        cdef PetscReal a = asReal(alpha)
        cdef PetscReal b = asReal(beta)
        CHKERR( PEPSTOARSetLinearization(self.pep, a, b) )

    def getSTOARLinearization(self) -> tuple[float, float]:
        """
        Get the coefficients that define the linearization of a quadratic eigenproblem.

        Not collective.

        Returns
        -------
        alpha: float
            First parameter of the linearization.
        beta: float
            Second parameter of the linearization.

        See Also
        --------
        setSTOARLinearization, slepc.PEPSTOARGetLinearization
        """
        cdef PetscReal a = 0.0
        cdef PetscReal b = 0.0
        CHKERR( PEPSTOARGetLinearization(self.pep, &a, &b) )
        return (asReal(a), asReal(b))

    def setSTOARLocking(self, lock: bool) -> None:
        """
        Toggle between locking and non-locking variants of the STOAR method.

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
        getSTOARLocking, slepc.PEPSTOARSetLocking
        """
        cdef PetscBool val = asBool(lock)
        CHKERR( PEPSTOARSetLocking(self.pep, val) )

    def getSTOARLocking(self) -> bool:
        """
        Get the locking flag used in the STOAR method.

        Not collective.

        Returns
        -------
        bool
            The locking flag.

        See Also
        --------
        setSTOARLocking, slepc.PEPSTOARGetLocking
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( PEPSTOARGetLocking(self.pep, &tval) )
        return toBool(tval)

    def setSTOARDetectZeros(self, detect: bool) -> None:
        """
        Set flag to enforce detection of zeros during the factorizations.

        Logically collective.

        Set a flag to enforce detection of zeros during the factorizations
        throughout the spectrum slicing computation.

        Parameters
        ----------
        detect
            ``True`` if zeros must checked for.

        Notes
        -----
        This call makes sense only for spectrum slicing runs, that is, when
        an interval has been given with `setInterval()` and `SINVERT` is set.

        A zero in the factorization indicates that a shift coincides with
        an eigenvalue.

        This flag is turned off by default, and may be necessary in some cases.
        This feature currently requires an external package for factorizations
        with support for zero detection, e.g. MUMPS.

        See Also
        --------
        setInterval, getSTOARDetectZeros, slepc.PEPSTOARSetDetectZeros
        """
        cdef PetscBool val = asBool(detect)
        CHKERR( PEPSTOARSetDetectZeros(self.pep, val) )

    def getSTOARDetectZeros(self) -> bool:
        """
        Get the flag that enforces zero detection in spectrum slicing.

        Not collective.

        Returns
        -------
        bool
            The zero detection flag.

        See Also
        --------
        setSTOARDetectZeros, slepc.PEPSTOARGetDetectZeros
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( PEPSTOARGetDetectZeros(self.pep, &tval) )
        return toBool(tval)

    def setSTOARDimensions(
        self,
        nev: int | None = None,
        ncv: int | None = None,
        mpd: int | None = None,
    ) -> None:
        """
        Set the dimensions used for each subsolve step.

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
        This call makes sense only for spectrum slicing runs, that is, when
        an interval has been given with `setInterval()` and `SINVERT` is set.

        The meaning of the parameters is the same as in `setDimensions()`, but
        the ones here apply to every subsolve done by the child `PEP` object.

        See Also
        --------
        setInterval, setDimensions, getSTOARDimensions, slepc.PEPSTOARSetDimensions
        """
        cdef PetscInt ival1 = PETSC_CURRENT
        cdef PetscInt ival2 = PETSC_CURRENT
        cdef PetscInt ival3 = PETSC_CURRENT
        if nev is not None: ival1 = asInt(nev)
        if ncv is not None: ival2 = asInt(ncv)
        if mpd is not None: ival3 = asInt(mpd)
        CHKERR( PEPSTOARSetDimensions(self.pep, ival1, ival2, ival3) )

    def getSTOARDimensions(self) -> tuple[int, int, int]:
        """
        Get the dimensions used for each subsolve step.

        Not collective.

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
        setSTOARDimensions, slepc.PEPSTOARGetDimensions
        """
        cdef PetscInt ival1 = 0
        cdef PetscInt ival2 = 0
        cdef PetscInt ival3 = 0
        CHKERR( PEPSTOARGetDimensions(self.pep, &ival1, &ival2, &ival3) )
        return (toInt(ival1), toInt(ival2), toInt(ival3))

    def getSTOARInertias(self) -> tuple[ArrayReal, ArrayInt]:
        """
        Get the values of the shifts and their corresponding inertias.

        Not collective.

        Get the values of the shifts and their corresponding inertias
        in case of doing spectrum slicing for a computational interval.

        Returns
        -------
        shifts: ArrayReal
            The values of the shifts used internally in the solver.
        inertias: ArrayInt
            The values of the inertia in each shift.

        Notes
        -----
        This call makes sense only for spectrum slicing runs, that is, when
        an interval has been given with `setInterval()` and `SINVERT` is set.

        If called after `solve()`, all shifts used internally by the solver are
        returned (including both endpoints and any intermediate ones). If called
        before `solve()` and after `setUp()` then only the information of the
        endpoints of subintervals is available.

        See Also
        --------
        setInterval, slepc.PEPSTOARGetInertias
        """
        cdef PetscReal *shiftsarray = NULL
        cdef PetscInt *inertiasarray = NULL
        cdef PetscInt n = 0
        CHKERR(PEPSTOARGetInertias(self.pep, &n, &shiftsarray, &inertiasarray))
        cdef object shifts = None
        cdef object inertias = None
        try:
            shifts = array_r(n, shiftsarray)
            inertias = array_i(n, inertiasarray)
        finally:
            CHKERR( PetscFree(shiftsarray) )
            CHKERR( PetscFree(inertiasarray) )
        return (shifts, inertias)

    def setSTOARCheckEigenvalueType(self, flag: bool) -> None:
        """
        Set flag to check if all eigenvalues have the same definite type.

        Logically collective.

        Set a flag to check that all the eigenvalues obtained throughout the
        spectrum slicing computation have the same definite type.

        Parameters
        ----------
        flag
            Whether the eigenvalue type is checked or not.

        Notes
        -----
        This option is relevant only for spectrum slicing computations, but
        is ignored in `slepc4py.SLEPc.PEP.ProblemType.HYPERBOLIC` problems.

        This flag is turned on by default, to guarantee that the computed
        eigenvalues have the same type (otherwise the computed solution might
        be wrong). But since the check is computationally quite expensive,
        the check may be turned off if the user knows for sure that all
        eigenvalues in the requested interval have the same type.

        See Also
        --------
        getSTOARCheckEigenvalueType, slepc.PEPSTOARSetCheckEigenvalueType
        """
        cdef PetscBool sval = asBool(flag)
        CHKERR( PEPSTOARSetCheckEigenvalueType(self.pep, sval) )

    def getSTOARCheckEigenvalueType(self) -> bool:
        """
        Get the flag for the eigenvalue type check in spectrum slicing.

        Not collective.

        Returns
        -------
        bool
            Whether the eigenvalue type is checked or not.

        See Also
        --------
        setSTOARCheckEigenvalueType, slepc.PEPSTOARGetCheckEigenvalueType
        """
        cdef PetscBool sval = PETSC_FALSE
        CHKERR( PEPSTOARGetCheckEigenvalueType(self.pep, &sval) )
        return toBool(sval)

    #

    def setJDRestart(self, keep: float) -> None:
        """
        Set the restart parameter for the Jacobi-Davidson method.

        Logically collective.

        Set the restart parameter for the Jacobi-Davidson method, in
        particular the proportion of basis vectors that must be kept
        after restart.

        Parameters
        ----------
        keep
            The number of vectors to be kept at restart.

        Notes
        -----
        Allowed values are in the range [0.1,0.9]. The default is 0.5.

        See Also
        --------
        getJDRestart, slepc.PEPJDSetRestart
        """
        cdef PetscReal val = asReal(keep)
        CHKERR( PEPJDSetRestart(self.pep, val) )

    def getJDRestart(self) -> float:
        """
        Get the restart parameter used in the Jacobi-Davidson method.

        Not collective.

        Returns
        -------
        float
            The number of vectors to be kept at restart.

        See Also
        --------
        setJDRestart, slepc.PEPJDGetRestart
        """
        cdef PetscReal val = 0
        CHKERR( PEPJDGetRestart(self.pep, &val) )
        return toReal(val)

    def setJDFix(self, fix: float) -> None:
        """
        Set the threshold for changing the target in the correction equation.

        Logically collective.

        Parameters
        ----------
        fix
            Threshold for changing the target.

        Notes
        -----
        The target in the correction equation is fixed at the first iterations.
        When the norm of the residual vector is lower than the fix value,
        the target is set to the corresponding eigenvalue.

        See Also
        --------
        getJDFix, slepc.PEPJDSetFix
        """
        cdef PetscReal val = asReal(fix)
        CHKERR( PEPJDSetFix(self.pep, val) )

    def getJDFix(self) -> float:
        """
        Get threshold for changing the target in the correction equation.

        Not collective.

        Returns
        -------
        float
            The threshold for changing the target.

        See Also
        --------
        setJDFix, slepc.PEPJDGetFix
        """
        cdef PetscReal val = 0
        CHKERR( PEPJDGetFix(self.pep, &val) )
        return toReal(val)

    def setJDReusePreconditioner(self, flag: bool) -> None:
        """
        Set a flag indicating whether the preconditioner must be reused or not.

        Logically collective.

        Parameters
        ----------
        flag
            The reuse flag.

        Notes
        -----
        The default value is ``False``. If set to ``True``, the
        preconditioner is built only at the beginning, using the
        target value. Otherwise, it may be rebuilt (depending on
        the ``fix`` parameter) at each iteration from the Ritz value.

        See Also
        --------
        setJDFix, getJDReusePreconditioner, slepc.PEPJDSetReusePreconditioner
        """
        cdef PetscBool bval = asBool(flag)
        CHKERR( PEPJDSetReusePreconditioner(self.pep, bval) )

    def getJDReusePreconditioner(self) -> bool:
        """
        Get the flag for reusing the preconditioner.

        Not collective.

        Returns
        -------
        bool
            The reuse flag.

        See Also
        --------
        setJDReusePreconditioner, slepc.PEPJDGetReusePreconditioner
        """
        cdef PetscBool bval = PETSC_FALSE
        CHKERR( PEPJDGetReusePreconditioner(self.pep, &bval) )
        return toBool(bval)

    def setJDMinimalityIndex(self, flag: int) -> None:
        """
        Set the maximum allowed value for the minimality index.

        Logically collective.

        Parameters
        ----------
        flag
            The maximum minimality index.

        Notes
        -----
        The default value is equal to the degree of the polynomial. A
        smaller value can be used if the wanted eigenvectors are known
        to be linearly independent.

        See Also
        --------
        getJDMinimalityIndex, slepc.PEPJDSetMinimalityIndex
        """
        cdef PetscInt ival = asInt(flag)
        CHKERR( PEPJDSetMinimalityIndex(self.pep, ival) )

    def getJDMinimalityIndex(self) -> int:
        """
        Get the maximum allowed value of the minimality index.

        Not collective.

        Returns
        -------
        int
            The maximum minimality index.

        See Also
        --------
        setJDMinimalityIndex, slepc.PEPJDGetMinimalityIndex
        """
        cdef PetscInt ival = 1
        CHKERR( PEPJDGetMinimalityIndex(self.pep, &ival) )
        return toInt(ival)

    def setJDProjection(self, proj: JDProjection) -> None:
        """
        Set the type of projection to be used in the Jacobi-Davidson solver.

        Logically collective.

        Parameters
        ----------
        proj
            The type of projection.

        See Also
        --------
        getJDProjection, slepc.PEPJDSetProjection
        """
        cdef SlepcPEPJDProjection val = proj
        CHKERR( PEPJDSetProjection(self.pep, val) )

    def getJDProjection(self) -> JDProjection:
        """
        Get the type of projection to be used in the Jacobi-Davidson solver.

        Not collective.

        Returns
        -------
        JDProjection
            The type of projection.

        See Also
        --------
        setJDProjection, slepc.PEPJDGetProjection
        """
        cdef SlepcPEPJDProjection val = PEP_JD_PROJECTION_HARMONIC
        CHKERR( PEPJDGetProjection(self.pep, &val) )
        return val

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
        getCISSExtraction, slepc.PEPCISSSetExtraction
        """
        cdef SlepcPEPCISSExtraction val = extraction
        CHKERR( PEPCISSSetExtraction(self.pep, val) )

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
        setCISSExtraction, slepc.PEPCISSGetExtraction
        """
        cdef SlepcPEPCISSExtraction val = PEP_CISS_EXTRACTION_RITZ
        CHKERR( PEPCISSGetExtraction(self.pep, &val) )
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
        `petsc4py.PETSc.KSP` object is shared among all processes of the `PEP`
        communicator. Otherwise, the communicator is split into ``npart``
        communicators, so that ``npart`` `petsc4py.PETSc.KSP` solves proceed
        simultaneously.

        See Also
        --------
        getCISSSizes, setCISSThreshold, setCISSRefinement, slepc.PEPCISSSetSizes
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
        CHKERR( PEPCISSSetSizes(self.pep, ival1, ival2, ival3, ival4, ival5, bval) )

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
        setCISSSizes, slepc.PEPCISSGetSizes
        """
        cdef PetscInt  ival1 = 0
        cdef PetscInt  ival2 = 0
        cdef PetscInt  ival3 = 0
        cdef PetscInt  ival4 = 0
        cdef PetscInt  ival5 = 0
        cdef PetscBool bval  = PETSC_FALSE
        CHKERR( PEPCISSGetSizes(self.pep, &ival1, &ival2, &ival3, &ival4, &ival5, &bval) )
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
        getCISSThreshold, slepc.PEPCISSSetThreshold
        """
        cdef PetscReal rval1 = PETSC_CURRENT
        cdef PetscReal rval2 = PETSC_CURRENT
        if delta is not None: rval1 = asReal(delta)
        if spur  is not None: rval2 = asReal(spur)
        CHKERR( PEPCISSSetThreshold(self.pep, rval1, rval2) )

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
        setCISSThreshold, slepc.PEPCISSGetThreshold
        """
        cdef PetscReal delta = 0
        cdef PetscReal spur  = 0
        CHKERR( PEPCISSGetThreshold(self.pep, &delta, &spur) )
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
        getCISSRefinement, slepc.PEPCISSSetRefinement
        """
        cdef PetscInt ival1 = PETSC_CURRENT
        cdef PetscInt ival2 = PETSC_CURRENT
        if inner  is not None: ival1 = asInt(inner)
        if blsize is not None: ival2 = asInt(blsize)
        CHKERR( PEPCISSSetRefinement(self.pep, ival1, ival2) )

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
        setCISSRefinement, slepc.PEPCISSGetRefinement
        """
        cdef PetscInt ival1 = 0
        cdef PetscInt ival2 = 0
        CHKERR( PEPCISSGetRefinement(self.pep, &ival1, &ival2) )
        return (toInt(ival1), toInt(ival2))

    def getCISSKSPs(self) -> list[KSP]:
        """
        Get the array of linear solver objects associated with the CISS solver.

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
        setCISSSizes, slepc.PEPCISSGetKSPs
        """
        cdef PetscInt i = 0, n = 0
        cdef PetscKSP *p = NULL
        CHKERR( PEPCISSGetKSPs(self.pep, &n, &p) )
        return [ref_KSP(p[i]) for i from 0 <= i <n]

    #

    property problem_type:
        """The type of the eigenvalue problem."""
        def __get__(self) -> PEPProblemType:
            return self.getProblemType()
        def __set__(self, value):
            self.setProblemType(value)

    property which:
        """The portion of the spectrum to be sought."""
        def __get__(self) -> PEPWhich:
            return self.getWhichEigenpairs()
        def __set__(self, value):
            self.setWhichEigenpairs(value)

    property target:
        """The value of the target."""
        def __get__(self) -> float:
            return self.getTarget()
        def __set__(self, value):
            self.setTarget(value)

    property extract:
        """The type of extraction technique to be employed."""
        def __get__(self) -> PEPExtract:
            return self.getExtract()
        def __set__(self, value):
            self.setExtract(value)

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

    property track_all:
        """Compute the residual norm of all approximate eigenpairs."""
        def __get__(self) -> bool:
            return self.getTrackAll()
        def __set__(self, value):
            self.setTrackAll(value)

    property st:
        """The spectral transformation (`ST`) object associated."""
        def __get__(self) -> ST:
            return self.getST()
        def __set__(self, value):
            self.setST(value)

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

del PEPType
del PEPProblemType
del PEPWhich
del PEPBasis
del PEPScale
del PEPRefine
del PEPRefineScheme
del PEPExtract
del PEPErrorType
del PEPConv
del PEPStop
del PEPConvergedReason
del PEPJDProjection
del PEPCISSExtraction

# -----------------------------------------------------------------------------
