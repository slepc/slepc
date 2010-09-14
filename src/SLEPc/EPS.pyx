# -----------------------------------------------------------------------------

class EPSType(object):
    """
    EPS type

    Native sparse eigensolvers.

    - `KRYLOVSCHUR`:  Krylov-Schur (default).
    - `GD`:           Generalized Davidson.
    - `JD`:           Jacobi-Davidson.
    - `LANCZOS`:      Lanczos.
    - `ARNOLDI`:      Arnoldi.
    - `SUBSPACE`:     Subspace Iteration.
    - `POWER`:        Power Iteration, Inverse Iteration, RQI.
    - `LAPACK`:       Wrappers to dense eigensolvers in Lapack.

    Wrappers to sparse eigensolvers
    (should be enabled during installation of SLEPc)

    - `ARPACK`:
    - `BLZPACK`:
    - `TRLAN`:
    - `BLOPEX`:
    - `PRIMME`:
    """
    # provided implementations
    POWER        = S_(EPSPOWER)
    SUBSPACE     = S_(EPSSUBSPACE)
    ARNOLDI      = S_(EPSARNOLDI)
    LANCZOS      = S_(EPSLANCZOS)
    KRYLOVSCHUR  = S_(EPSKRYLOVSCHUR)
    GD           = S_(EPSGD)
    JD           = S_(EPSJD)
    LAPACK       = S_(EPSLAPACK)
    # with external libraries
    ARPACK       = S_(EPSARPACK)
    BLZPACK      = S_(EPSBLZPACK)
    TRLAN        = S_(EPSTRLAN)
    BLOPEX       = S_(EPSBLOPEX)
    PRIMME       = S_(EPSPRIMME)

class EPSProblemType(object):
    """
    EPS problem type

    - `HEP`:    Hermitian eigenproblem.
    - `NHEP`:   Non-Hermitian eigenproblem.
    - `GHEP`:   Generalized Hermitian eigenproblem.
    - `GNHEP`:  Generalized Non-Hermitian eigenproblem.
    - `PGNHEP`: Generalized Non-Hermitian eigenproblem 
                with positive definite ``B``.
    """
    HEP    = EPS_HEP
    NHEP   = EPS_NHEP
    GHEP   = EPS_GHEP
    GNHEP  = EPS_GNHEP
    PGNHEP = EPS_PGNHEP

class EPSExtraction(object):
    """
    EPS extraction technique

    - `RITZ`:              Standard Rayleigh-Ritz extraction.
    - `HARMONIC`:          Harmonic extraction.
    - `HARMONIC_RELATIVE`: Harmonic extraction relative to the eigenvalue.
    - `HARMONIC_RIGHT`:    Harmonic extraction for rightmost eigenvalues.
    - `HARMONIC_LARGEST`:  Harmonic extraction for largest magnitude (without target).
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
    EPS type of balancing used for non-Hermitian problems


    - `NONE`:     None.
    - `ONE_SIDE`: One-sided eigensolver, only right eigenvectors.
    - `TWO_SIDE`: Two-sided eigensolver, right and left eigenvectors.
    - `USER`:     User-defined.
    """
    NONE    = EPS_BALANCE_NONE
    ONESIDE = EPS_BALANCE_ONESIDE
    TWOSIDE = EPS_BALANCE_TWOSIDE
    USER    = EPS_BALANCE_USER

class EPSWhich(object):
    """
    EPS desired piece of spectrum

    - `LARGEST_MAGNITUDE`:  Largest magnitude (default).
    - `LARGEST_REAL`:       Largest real parts.
    - `LARGEST_IMAGINARY`:  Largest imaginary parts in magnitude.
    - `SMALLEST_MAGNITUDE`: Smallest magnitude.
    - `SMALLEST_REAL`:      Smallest real parts.
    - `SMALLEST_IMAGINARY`: Smallest imaginary parts in magnitude.
    - `TARGET_MAGNITUDE`:   Closest to target (in magnitude).
    - `TARGET_REAL`:        Real part closest to target.
    - `TARGET_IMAGINARY`:   Imaginary part closest to target.
    - `USER`:               User defined ordering.
    """
    LARGEST_MAGNITUDE  = EPS_LARGEST_MAGNITUDE
    LARGEST_REAL       = EPS_LARGEST_REAL
    LARGEST_IMAGINARY  = EPS_LARGEST_IMAGINARY
    SMALLEST_MAGNITUDE = EPS_SMALLEST_MAGNITUDE
    SMALLEST_REAL      = EPS_SMALLEST_REAL
    SMALLEST_IMAGINARY = EPS_SMALLEST_IMAGINARY
    TARGET_MAGNITUDE   = EPS_TARGET_MAGNITUDE
    TARGET_REAL        = EPS_TARGET_REAL
    TARGET_IMAGINARY   = EPS_TARGET_IMAGINARY
    USER               = EPS_WHICH_USER

class EPSConvergedReason(object):
    """
    EPS convergence reasons

    - `CONVERGED_TOL`:
    - `DIVERGED_ITS`:
    - `DIVERGED_BREAKDOWN`:
    - `DIVERGED_NONSYMMETRIC`:
    - `CONVERGED_ITERATING`:
    """
    CONVERGED_TOL         = EPS_CONVERGED_TOL
    DIVERGED_ITS          = EPS_DIVERGED_ITS
    DIVERGED_BREAKDOWN    = EPS_DIVERGED_BREAKDOWN
    DIVERGED_NONSYMMETRIC = EPS_DIVERGED_NONSYMMETRIC
    CONVERGED_ITERATING   = EPS_CONVERGED_ITERATING
    ITERATING             = EPS_CONVERGED_ITERATING

class EPSPowerShiftType(object):
    """
    EPS Power shift type.

    - `CONSTANT`:
    - `RAYLEIGH`:
    - `WILKINSON`:
    """
    CONSTANT  = EPS_POWER_SHIFT_CONSTANT
    RAYLEIGH  = EPS_POWER_SHIFT_RAYLEIGH
    WILKINSON = EPS_POWER_SHIFT_WILKINSON

class EPSLanczosReorthogType(object):
    """
    EPS Lanczos reorthogonalization type

    - `LOCAL`:
    - `FULL`:
    - `SELECTIVE`:
    - `PERIODIC`:
    - `PARTIAL`:
    - `DELAYED`:
    """
    LOCAL     =  EPS_LANCZOS_REORTHOG_LOCAL
    FULL      =  EPS_LANCZOS_REORTHOG_FULL
    SELECTIVE =  EPS_LANCZOS_REORTHOG_SELECTIVE
    PERIODIC  =  EPS_LANCZOS_REORTHOG_PERIODIC
    PARTIAL   =  EPS_LANCZOS_REORTHOG_PARTIAL
    DELAYED   =  EPS_LANCZOS_REORTHOG_DELAYED

# --------------------------------------------------------------------

cdef class EPS(Object):

    """
    EPS
    """

    Type            = EPSType
    ProblemType     = EPSProblemType
    Extraction      = EPSExtraction
    Balance         = EPSBalance
    Which           = EPSWhich
    ConvergedReason = EPSConvergedReason

    PowerShiftType      = EPSPowerShiftType
    LanczosReorthogType = EPSLanczosReorthogType

    def __cinit__(self):
        self.obj = <PetscObject*> &self.eps
        self.eps = NULL

    def view(self, Viewer viewer=None):
        """
        Prints the EPS data structure.

        Parameters
        ----------
        viewer: Viewer, optional.
                Visualization context; if not provided, the standard
                output is used.
        """
        cdef PetscViewer vwr = NULL
        if viewer is not None: vwr = viewer.vwr
        CHKERR( EPSView(self.eps, vwr) )

    def destroy(self):
        """
        Destroys the EPS object.
        """
        CHKERR( EPSDestroy(self.eps) )
        self.eps = NULL
        return self

    def create(self, comm=None):
        """
        Creates the EPS object.

        Parameters
        ----------
        comm: MPI_Comm, optional
              MPI communicator; if not provided, it defaults to all
              processes.
        """
        cdef MPI_Comm ccomm = def_Comm(comm, SLEPC_COMM_DEFAULT())
        cdef SlepcEPS neweps = NULL
        CHKERR( EPSCreate(ccomm, &neweps) )
        self.dec_ref(); self.eps = neweps
        return self

    def setType(self, eps_type):
        """
        Selects the particular solver to be used in the EPS object.

        Parameters
        ----------
        eps_type: EPS.Type enumerate
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

    def getType(self):
        """
        Gets the EPS type of this object.

        Returns
        -------
        type: EPS.Type enumerate
              The solver currently being used.
        """
        cdef SlepcEPSType eps_type = NULL
        CHKERR( EPSGetType(self.eps, &eps_type) )
        return bytes2str(eps_type)

    def getOptionsPrefix(self):
        """
        Gets the prefix used for searching for all EPS options in the
        database.

        Returns
        -------
        prefix: string
                The prefix string set for this EPS object.
        """
        cdef const_char *prefix = NULL
        CHKERR( EPSGetOptionsPrefix(self.eps, &prefix) )
        return bytes2str(prefix)

    def setOptionsPrefix(self, prefix):
        """
        Sets the prefix used for searching for all EPS options in the
        database.

        Parameters
        ----------
        prefix: string
                The prefix string to prepend to all EPS option
                requests.

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
        cdef const_char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( EPSSetOptionsPrefix(self.eps, cval) )

    def appendOptionsPrefix(self, prefix):
        """
        Appends to the prefix used for searching for all EPS options
        in the database.

        Parameters
        ----------
        prefix: string
                The prefix string to prepend to all EPS option requests.
        """
        cdef const_char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( EPSAppendOptionsPrefix(self.eps, cval) )

    def setFromOptions(self):
        """
        Sets EPS options from the options database. This routine must
        be called before `setUp()` if the user is to be allowed to set
        the solver type.

        Notes
        -----
        To see all options, run your program with the ``-help``
        option.
        """
        CHKERR( EPSSetFromOptions(self.eps) )

    #

    def getProblemType(self):
        """
        Gets the problem type from the EPS object.

        Returns
        -------
        problem_type: EPS.ProblemType enumerate
                      The problem type that was previously set.
        """
        cdef SlepcEPSProblemType val = EPS_NHEP
        CHKERR( EPSGetProblemType(self.eps, &val) )
        return val

    def setProblemType(self, problem_type):
        """
        Specifies the type of the eigenvalue problem.

        Parameters
        ----------
        problem_type: EPS.ProblemType enumerate
               The problem type to be set.

        Notes
        -----
        Allowed values are: Hermitian (HEP), non-Hermitian (NHEP),
        generalized Hermitian (GHEP), generalized non-Hermitian
        (GNHEP), and generalized non-Hermitian with positive
        semi-definite B (PGNHEP).

        This function must be used to instruct SLEPc to exploit
        symmetry. If no problem type is specified, by default a
        non-Hermitian problem is assumed (either standard or
        generalized). If the user knows that the problem is Hermitian
        (i.e. ``A=A^H``) or generalized Hermitian (i.e. ``A=A^H``,
        ``B=B^H``, and ``B`` positive definite) then it is recommended
        to set the problem type so that eigensolver can exploit these
        properties.
        """
        cdef SlepcEPSProblemType val = problem_type
        CHKERR( EPSSetProblemType(self.eps, val) )

    def isGeneralized(self):
        """
        Tells whether the EPS object corresponds to a generalized
        eigenvalue problem.

        Returns
        -------
        flag: boolean
              True if two matrices were set with `setOperators()`.
        """
        cdef PetscTruth tval = PETSC_FALSE
        CHKERR( EPSIsGeneralized(self.eps, &tval) )
        return <bint> tval

    def isHermitian(self):
        """
        Tells whether the EPS object corresponds to a Hermitian
        eigenvalue problem.

        Returns
        -------
        flag: boolean
              True if the problem type set with `setProblemType()` was
              Hermitian.
        """
        cdef PetscTruth tval = PETSC_FALSE
        CHKERR( EPSIsHermitian(self.eps, &tval) )
        return <bint> tval

    def getBalance(self):
        """
        Gets the balancing type used by the EPS object,
        and the associated parameters.

        Returns
        -------
        balance: EPS.Balance enumerate
                 The balancing method
        iterations: integer
                    Number of iterations of the balancing algorithm
        cutoff: real
                Cutoff value
        """
        cdef SlepcEPSBalance val = EPS_BALANCE_ONESIDE
        cdef PetscInt ival = 0
        cdef PetscReal rval = 0
        CHKERR( EPSGetBalance(self.eps, &val, &ival, &rval) )
        return (val, toInt(ival), toReal(rval))

    def setBalance(self, balance=None, iterations=None, cutoff=None):
        """
        Specifies the balancing technique to be employed by the
        eigensolver, and some parameters associated to it.

        Parameters
        ----------
        balance: EPS.Balance enumerate
                 The balancing method
        iterations: integer
                    Number of iterations of the balancing algorithm
        cutoff: real
                Cutoff value
        """
        cdef SlepcEPSBalance val = <SlepcEPSBalance>PETSC_IGNORE
        cdef PetscInt ival = PETSC_IGNORE
        cdef PetscReal rval = PETSC_IGNORE
        if balance    is not None: val  = balance
        if iterations is not None: ival = asInt(iterations)
        if cutoff     is not None: rval = asReal(cutoff)
        CHKERR( EPSSetBalance(self.eps, val, ival, rval) )

    def getExtraction(self):
        """
        Gets the extraction type used by the EPS object.

        Returns
        -------
        extraction: EPS.Extraction enumerate
                    The method of extraction.
        """
        cdef SlepcEPSExtraction val = EPS_RITZ
        CHKERR( EPSGetExtraction(self.eps, &val) )
        return val

    def setExtraction(self, extraction):
        """
        Sets the extraction type used by the EPS object.

        Parameters
        ----------
        extraction: EPS.Extraction enumerate
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

    def getWhichEigenpairs(self):
        """
        Returns which portion of the spectrum is to be sought.

        Returns
        -------
        which: EPS.Which enumerate
               The portion of the spectrum to be sought by the solver.
        """
        cdef SlepcEPSWhich val = EPS_LARGEST_MAGNITUDE
        CHKERR( EPSGetWhichEigenpairs(self.eps, &val) )
        return val

    def setWhichEigenpairs(self, which):
        """
        Specifies which portion of the spectrum is to be sought.

        Parameters
        ----------
        which: EPS.Which enumerate
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

    def getLeftVectorsWanted(self):
        """
        Returns the flag indicating whether left eigenvectors are
        required or not.

        Returns
        -------
        wanted: boolean
                Whether left eigenvectors are required or not.
        """
        cdef PetscTruth tval = PETSC_FALSE
        CHKERR( EPSGetLeftVectorsWanted(self.eps, &tval) )
        return <bint>tval

    def setLeftVectorsWanted(self, wanted):
        """
        Specifies the flag indicating whether left eigenvectors are
        required or not.

        Parameters
        ----------
        wanted: boolean
                Whether left eigenvectors are required or not.
        """
        cdef PetscTruth tval = wanted
        CHKERR( EPSSetLeftVectorsWanted(self.eps, tval) )

    def getTarget(self):
        """
        Gets the value of the target.

        Returns
        -------
        target: float (real or complex)
                The value of the target.

        Notes
        -----
        If the target was not set by the user, then zero is returned.
        """
        cdef PetscScalar sval = 0
        CHKERR( EPSGetTarget(self.eps, &sval) )
        return toScalar(sval)

    def setTarget(self, target):
        """
        Sets the value of the target.

        Parameters
        ----------
        target: float (real or complex)
                The value of the target.

        Notes
        -----
        The target is a scalar value used to determine the portion of
        the spectrum of interest. It is used in combination with
        `setWhichEigenpairs()`.
        """
        cdef PetscScalar sval = asScalar(target)
        CHKERR( EPSSetTarget(self.eps, sval) )

    #

    def getTolerances(self):
        """
        Gets the tolerance and maximum iteration count used by the
        default EPS convergence tests.

        Returns
        -------
        tol: float
             The convergence tolerance.
        max_it: int
             The maximum number of iterations
        """
        cdef PetscReal rval = 0
        cdef PetscInt  ival = 0
        CHKERR( EPSGetTolerances(self.eps, &rval, &ival) )
        return (toReal(rval), toInt(ival))

    def setTolerances(self, tol=None, max_it=None):
        """
        Sets the tolerance and maximum iteration count used by the
        default EPS convergence tests.

        Parameters
        ----------
        tol: float, optional
             The convergence tolerance.
        max_it: int, optional
             The maximum number of iterations

        Notes
        -----
        Use `DECIDE` for maxits to assign a reasonably good value,
        which is dependent on the solution method.
        """
        cdef PetscReal rval = PETSC_IGNORE
        cdef PetscInt  ival = PETSC_IGNORE
        if tol    is not None: rval = asReal(tol)
        if max_it is not None: ival = asInt(max_it)
        CHKERR( EPSSetTolerances(self.eps, rval, ival) )

    def getTrackAll(self):
        """
        Returns the flag indicating whether all residual norms must be
        computed or not.

        Returns
        -------
        trackall: bool
            Whether the solver compute all residuals or not.
        """
        cdef PetscTruth tval = PETSC_FALSE
        CHKERR( EPSGetTrackAll(self.eps, &tval) )
        return <bint>tval

    def setTrackAll(self, trackall):
        """
        Specifies if the solver must compute the residual of all
        approximate eigenpairs or not.

        Parameters
        ----------
        trackall: bool
            Whether compute all residuals or not.
        """
        cdef PetscTruth tval = trackall
        CHKERR( EPSSetTrackAll(self.eps, tval) )

    def getDimensions(self):
        """
        Gets the number of eigenvalues to compute and the dimension of
        the subspace.

        Returns
        -------
        nev: int
             Number of eigenvalues to compute.
        ncv: int
             Maximum dimension of the subspace to be used by the
             solver.
        mpd: int
             Maximum dimension allowed for the projected problem.
        """
        cdef PetscInt ival1 = 0
        cdef PetscInt ival2 = 0
        cdef PetscInt ival3 = 0
        CHKERR( EPSGetDimensions(self.eps, &ival1, &ival2, &ival3) )
        return (toInt(ival1), toInt(ival2), toInt(ival3))

    def setDimensions(self, nev=None, ncv=None, mpd=None):
        """
        Sets the number of eigenvalues to compute and the dimension of
        the subspace.

        Parameters
        ----------
        nev: int, optional
             Number of eigenvalues to compute.
        ncv: int, optional
             Maximum dimension of the subspace to be used by the
             solver.
        mpd: int, optional
             Maximum dimension allowed for the projected problem.

        Notes
        -----
        Use `DECIDE` for `ncv` and `mpd` to assign a reasonably good
        value, which is dependent on the solution method.

        The parameters `ncv` and `mpd` are intimately related, so that
        the user is advised to set one of them at most. Normal usage
        is the following:

        + In cases where `nev` is small, the user sets `ncv`
          (a reasonable default is 2 * `nev`).

        + In cases where `nev` is large, the user sets `mpd`.

        The value of `ncv` should always be between `nev` and (`nev` +
        `mpd`), typically `ncv` = `nev` + `mpd`. If `nev` is not too
        large, `mpd` = `nev` is a reasonable choice, otherwise a
        smaller value should be used.
        """
        cdef PetscInt ival1 = PETSC_IGNORE
        cdef PetscInt ival2 = PETSC_IGNORE
        cdef PetscInt ival3 = PETSC_IGNORE
        if nev is not None: ival1 = asInt(nev)
        if ncv is not None: ival2 = asInt(ncv)
        if mpd is not None: ival3 = asInt(mpd)
        CHKERR( EPSSetDimensions(self.eps, ival1, ival2, ival3) )

    def getST(self):
        """
        Obtain the spectral transformation (`ST`) object associated to
        the eigensolver object.

        Returns
        -------
        st: ST
            The spectral transformation.
        """
        cdef ST st = ST()
        CHKERR( EPSGetST(self.eps, &st.st) )
        st.inc_ref(); return st

    def setST(self, ST st not None):
        """
        Associates a spectral transformation object to the
        eigensolver.

        Parameters
        ----------
        st: ST
            The spectral transformation.
        """
        CHKERR( EPSSetST(self.eps, st.st) )

    def getIP(self):
        """
        Obtain the inner product associated to the eigensolver.

        Returns
        -------
        ip: IP
            The inner product context.
        """
        cdef IP ip = IP()
        CHKERR( EPSGetIP(self.eps, &ip.ip) )
        ip.inc_ref(); return ip

    def setIP(self, IP ip not None):
        """
        Associates an inner product to the eigensolver.

        Parameters
        ----------
        ip: IP
            The inner product context.
        """
        CHKERR( EPSSetIP(self.eps, ip.ip) )

    def getOperators(self):
        """
        Gets the matrices associated with the eigenvalue problem.

        Returns
        -------
        A: Mat
           The matrix associated with the eigensystem.
        B: Mat
           The second matrix in the case of generalized eigenproblems.
        """
        cdef Mat A = Mat()
        cdef Mat B = Mat()
        CHKERR( EPSGetOperators(self.eps, &A.mat, &B.mat) )
        A.inc_ref(); B.inc_ref(); return (A, B)

    def setOperators(self, Mat A not None, Mat B=None):
        """
        Sets the matrices associated with the eigenvalue problem.

        Parameters
        ----------
        A: Mat
           The matrix associated with the eigensystem.
        B: Mat, optional
           The second matrix in the case of generalized eigenproblems;
           if not provided, a standard eigenproblem is assumed.
        """
        cdef PetscMat Bmat = NULL
        if B is not None: Bmat = B.mat
        CHKERR( EPSSetOperators(self.eps, A.mat, Bmat) )

    def setDeflationSpace(self, space):
        """
        Add vectors to the basis of the deflation space.

        Parameters
        ----------
        space: a Vec or an array of Vec
               Set of basis vectors to be added to the deflation
               space.

        Notes
        -----
        When a deflation space is given, the eigensolver seeks the
        eigensolution in the restriction of the problem to the
        orthogonal complement of this space. This can be used for
        instance in the case that an invariant subspace is known
        beforehand (such as the nullspace of the matrix).

        Basis vectors set by a previous call to `setDeflationSpace()`
        are replaced.

        The vectors do not need to be mutually orthonormal, since they
        are explicitly orthonormalized internally.

        These vectors persist from one `solve()` call to the other,
        use `removeDeflationSpace()` to eliminate them. 
        """
        cdef PetscInt i = 0, nds = 0
        cdef PetscVec* vds = NULL
        if isinstance(space, Vec): space = [space]
        nds = len(space)
        cdef tmp = allocate(nds*sizeof(Vec),<void**>&vds)
        for i in range(nds): vds[i] = (<Vec?>space[<Py_ssize_t>i]).vec
        CHKERR( EPSSetDeflationSpace(self.eps, nds, vds) )

    def removeDeflationSpace(self):
        """
        Removes the deflation space previously set with
        `setDeflationSpace()`.
        """
        CHKERR( EPSRemoveDeflationSpace(self.eps) )

    #

    def setInitialSpace(self, space):
        """
        Sets the initial space from which the eigensolver starts to
        iterate.

        Parameters
        ----------
        space: Vec or sequence of Vec
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
        cdef PetscInt i = 0, ns = 0
        cdef PetscVec *vs = NULL
        if isinstance(space, Vec): space = [space]
        ns = len(space)
        cdef tmp = allocate(ns*sizeof(Vec),<void**>&vs)
        for i in range(ns): vs[i] = (<Vec?>space[<Py_ssize_t>i]).vec
        CHKERR( EPSSetInitialSpace(self.eps, ns, vs) )

    def setInitialSpaceLeft(self, space):
        """
        Sets a basis of vectors that constitute the initial left space, that
        is, the subspace from which the solver starts to iterate for building
        the left subspace (in methods that work with two subspaces).

        Parameters
        ----------
        space: Vec or sequence of Vec
           The initial space
        """
        cdef PetscInt i = 0, ns = 0
        cdef PetscVec *vs = NULL
        if isinstance(space, Vec): space = [space]
        ns = len(space)
        cdef tmp = allocate(ns*sizeof(Vec),<void**>&vs)
        for i in range(ns): vs[i] = (<Vec?>space[<Py_ssize_t>i]).vec
        CHKERR( EPSSetInitialSpaceLeft(self.eps, ns, vs) )

    #

    def cancelMonitor(self):
        """
        Clears all monitors for an EPS object.
        """
        CHKERR( EPSMonitorCancel(self.eps) )

    #

    def setUp(self):
        """
        Sets up all the internal data structures necessary for the
        execution of the eigensolver.

        Notes
        -----
        This function need not be called explicitly in most cases,
        since `solve()` calls it. It can be useful when one wants to
        measure the set-up time separately from the solve time.
        """
        CHKERR( EPSSetUp(self.eps) )

    def solve(self):
        """
        Solves the eigensystem.
        """
        CHKERR( EPSSolve(self.eps) )

    def getIterationNumber(self):
        """
        Gets the current iteration number. If the call to `solve()` is
        complete, then it returns the number of iterations carried out
        by the solution method.

        Returns
        -------
        its: int
             Iteration number.
        """
        cdef PetscInt ival = 0
        CHKERR( EPSGetIterationNumber(self.eps, &ival) )
        return toInt(ival)

    def getConvergedReason(self):
        """
        Gets the reason why the `solve()` iteration was stopped.

        Returns
        -------
        reason: EPS.ConvergedReason enumerate
                Negative value indicates diverged, positive value
                converged.
        """
        cdef SlepcEPSConvergedReason val = EPS_CONVERGED_ITERATING
        CHKERR( EPSGetConvergedReason(self.eps, &val) )
        return val

    def getConverged(self):
        """
        Gets the number of converged eigenpairs.

        Returns
        -------
        nconv: int
               Number of converged eigenpairs.

        Notes
        -----
        This function should be called after `solve()` has finished.
        """
        cdef PetscInt ival = 0
        CHKERR( EPSGetConverged(self.eps, &ival) )
        return toInt(ival)

    def getEigenvalue(self, int i):
        """
        Gets the i-th eigenvalue as computed by `solve()`.

        Parameters
        ----------
        i: int
           Index of the solution to be obtained.

        Returns
        -------
        e: scalar (possibly complex)
           The computed eigenvalue.

        Notes
        -----
        The index ``i`` should be a value between ``0`` and
        ``nconv-1`` (see `getConverged()`. Eigenpairs are indexed
        according to the ordering criterion established with
        `setWhichEigenpairs()`.
        """
        cdef PetscScalar sval1 = 0
        cdef PetscScalar sval2 = 0
        CHKERR( EPSGetEigenvalue(self.eps, i, &sval1, &sval2) )
        return complex(toScalar(sval1), toScalar(sval2))

    def getEigenvector(self, int i, Vec Vr not None, Vec Vi=None):
        """
        Gets the i-th eigenvector as computed by `solve()`.

        Parameters
        ----------
        i: int
           Index of the solution to be obtained.
        Vr: Vec
            Placeholder for the returned eigenvector (real part).
        Vi: Vec, optional
            Placeholder for the returned eigenvector (imaginary part).

        Notes
        -----
        The index ``i`` should be a value between ``0`` and
        ``nconv-1`` (see `getConverged()`. Eigenpairs are indexed
        according to the ordering criterion established with
        `setWhichEigenpairs()`.
        """
        cdef PetscVec vecr = NULL
        cdef PetscVec veci = NULL
        if Vr is not None: vecr = Vr.vec
        if Vi is not None: veci = Vi.vec
        CHKERR( EPSGetEigenvector(self.eps, i, vecr, veci) )

    def getEigenpair(self, int i, Vec Vr=None, Vec Vi=None):
        """
        Gets the i-th solution of the eigenproblem as computed by
        `solve()`.  The solution consists of both the eigenvalue and
        the eigenvector.

        Parameters
        ----------
        i: int
           Index of the solution to be obtained.
        Vr: Vec
            Placeholder for the returned eigenvector (real part).
        Vi: Vec
            Placeholder for the returned eigenvector (imaginary part).

        Returns
        -------
        e: scalar (possibly complex)
           The computed eigenvalue.

        Notes
        -----
        The index ``i`` should be a value between ``0`` and
        ``nconv-1`` (see `getConverged()`. Eigenpairs are indexed
        according to the ordering criterion established with
        `setWhichEigenpairs()`.
        """
        cdef PetscScalar sval1 = 0
        cdef PetscScalar sval2 = 0
        cdef PetscVec vecr = NULL
        cdef PetscVec veci = NULL
        if Vr is not None: vecr = Vr.vec
        if Vi is not None: veci = Vi.vec
        CHKERR( EPSGetEigenpair(self.eps, i, &sval1, &sval2, vecr, veci) )
        return complex(toScalar(sval1), toScalar(sval2))

    ## def getInvariantSubspace(self):
    ##     """
    ##     Gets an orthonormal basis of the computed invariant subspace.
    ##
    ##     Returns
    ##     -------
    ##     subspace: list of Vec
    ##        Basis of the invariant subspace.
    ##
    ##     Notes
    ##     -----
    ##     This function should be called after `solve()` has finished.
    ##
    ##     The returned vectors span an invariant subspace associated
    ##     with the computed eigenvalues. An invariant subspace ``X`` of
    ##     ``A` satisfies ``A x`` in ``X`` for all ``x`` in ``X`` (a
    ##     similar definition applies for generalized eigenproblems).
    ##     """
    ##     cdef PetscInt i = 0, ncv = 0
    ##     cdef PetscVec v = NULL, *isp = NULL
    ##     CHKERR( EPSGetConverged(self.eps, &ncv) )
    ##     CHKERR( EPSGetInitialVector(self.eps, &v) )
    ##     cdef Vec V = None
    ##     cdef list subspace = []
    ##     cdef object tmp = allocate(ncv*sizeof(Vec),<void**>&isp)
    ##     for i in range(ncv):
    ##         V = Vec(); subspace.append(V)
    ##         CHKERR( VecDuplicate(v, &isp[i]) )
    ##         V.vec = isp[i]
    ##     CHKERR( EPSGetInvariantSubspace(self.eps, isp) )
    ##     return subspace

    #

    def getEigenvectorLeft(self, int i, Vec Wr not None, Vec Wi=None):
        """
        Gets the i-th left eigenvector as computed by `solve()`,
        (only available in two-sided eigensolvers).

        Parameters
        ----------
        i: int
           Index of the solution to be obtained.
        Wr: Vec
            Placeholder for the returned left eigenvector (real part).
        Wi: Vec, optional
            Placeholder for the returned left eigenvector (imaginary
            part).

        Notes
        -----
        The index ``i`` should be a value between ``0`` and
        ``nconv-1`` (see `getConverged()`. Eigenpairs are indexed
        according to the ordering criterion established with
        `setWhichEigenpairs()`.
        """
        cdef PetscVec vecr = NULL
        cdef PetscVec veci = NULL
        if Wr is not None: vecr = Wr.vec
        if Wi is not None: veci = Wi.vec
        CHKERR( EPSGetEigenvectorLeft(self.eps, i, vecr, veci) )

    ## def getInvariantSubspaceLeft(self):
    ##     """
    ##     Gets an orthonormal basis of the computed left invariant
    ##     subspace (only available in two-sided eigensolvers).
    ##
    ##     Returns
    ##     -------
    ##     v: array of Vec
    ##        Basis of the left invariant subspace.
    ##
    ##     Notes
    ##     -----
    ##     See `getInvariantSubspace()` for additional information.
    ##     """
    ##     cdef PetscInt i = 0, ncv = 0
    ##     cdef PetscVec w = NULL, *isp = NULL
    ##     CHKERR( EPSGetConverged(self.eps, &ncv) )
    ##     CHKERR( EPSGetLeftInitialVector(self.eps, &w) )
    ##     cdef Vec W = None
    ##     cdef list subspace = []
    ##     cdef object tmp = allocate(ncv*sizeof(Vec),<void**>&isp)
    ##     for i in range(ncv):
    ##         W = Vec(); subspace.append(W)
    ##         CHKERR( VecDuplicate(w, &isp[i]) )
    ##         W.vec = isp[i]
    ##     CHKERR( EPSGetLeftInvariantSubspace(self.eps, isp) )
    ##     return subspace

    #

    def getErrorEstimate(self, int i):
        """
        Returns the error estimate associated to the i-th computed
        eigenpair.

        Parameters
        ----------
        i: int
           Index of the solution to be considered.

        Returns
        -------
        e: real
           Error estimate.

        Notes
        -----
        This is the error estimate used internally by the
        eigensolver. The actual error bound can be computed with
        `computeRelativeError()`.
        """
        cdef PetscReal rval = 0
        CHKERR( EPSGetErrorEstimate(self.eps, i, &rval) )
        return toReal(rval)

    def getErrorEstimateLeft(self, int i):
        """
        Returns the left error estimate associated to the i-th
        computed eigenpair (only available in two-sided eigensolvers).

        Parameters
        ----------
        i: int
           Index of the solution to be considered.

        Returns
        -------
        e: real
           Left error estimate.

        Notes
        -----
        This is the error estimate used internally by the
        eigensolver. The actual error bound can be computed with
        `computeRelativeError()`.
        """
        cdef PetscReal rval = 0
        CHKERR( EPSGetErrorEstimateLeft(self.eps, i, &rval) )
        return toReal(rval)

    def computeRelativeError(self, int i):
        """
        Computes the relative error bound associated with the i-th
        computed eigenpair.

        Parameters
        ----------
        i: int
           Index of the solution to be considered.

        Returns
        -------
        e: real
           The relative error bound, computed as
           ``||Ax-kBx||_2/||kx||_2`` where ``k`` is the eigenvalue and
           ``x`` is the eigenvector.  If ``k=0`` the relative error is
           computed as ``||Ax||_2/||x||_2``.

        Notes
        -----
        The index ``i`` should be a value between ``0`` and
        ``nconv-1`` (see `getConverged()`. Eigenpairs are indexed
        according to the ordering criterion established with
        `setWhichEigenpairs()`.
        """
        cdef PetscReal rval = 0
        CHKERR( EPSComputeRelativeError(self.eps, i, &rval) )
        return toReal(rval)

    def computeRelativeErrorLeft(self, int i):
        """
        Computes the left relative error bound associated with the
        i-th computed eigenpair (only available in two-sided
        eigensolvers).

        Parameters
        ----------
        i: int
           Index of the solution to be considered.

        Returns
        -------
        e: real
           The relative error bound, computed as
           ``||y'A-ky'B||_2/||ky||_2`` where ``k`` is the eigenvalue
           and ``y`` is the left eigenvector.  If ``k=0`` the relative
           error is computed as ``||y'A||_2/||y||_2``.
        """
        cdef PetscReal rval = 0
        CHKERR( EPSComputeRelativeErrorLeft(self.eps, i, &rval) )
        return toReal(rval)

    def computeResidualNorm(self, int i):
        """
        Computes the norm of the residual vector associated with the
        i-th computed eigenpair.

        Parameters
        ----------
        i: int
           Index of the solution to be considered.

        Returns
        -------
        norm: real
              The residual norm, computed as ``||Ax-kBx||_2`` where
              ``k`` is the eigenvalue and ``x`` is the eigenvector.
              If ``k=0`` then the residual norm is computed as
              ``||Ax||_2``.
        """
        cdef PetscReal rval = 0
        CHKERR( EPSComputeResidualNorm(self.eps, i, &rval) )
        return toReal(rval)

    def computeResidualNormLeft(self, int i):
        """
        Computes the norm of the residual vector associated with the
        i-th computed left eigenpair (only available in two-sided
        eigensolvers).

        Parameters
        ----------
        i: int
           Index of the solution to be considered.

        Returns
        -------

        norm: real
              The residual norm, computed as ``||y'A-ky'B||_2`` where
              ``k`` is the eigenvalue and ``y`` is the left
              eigenvector.  If ``k=0`` then the residual norm is
              computed as ``||y'A||_2``.
        """
        cdef PetscReal rval = 0
        CHKERR( EPSComputeResidualNormLeft(self.eps, i, &rval) )
        return toReal(rval)

    def getOperationCounters(self):
        """
        Gets the total number of operator applications, inner product
        operations and linear iterations used by the `ST` object
        during the last `solve()` call.

        Returns
        -------
        ops: int
             number of operator applications.
        dots: int
             number of inner product operations.
        lits: int
             number of linear iterations.

        Notes
        -----
        When the eigensolver algorithm invokes `ST.apply()` then a
        linear system must be solved (except in the case of standard
        eigenproblems and shift transformation). The number of
        iterations required in this solve is accumulated into a
        counter whose value is returned by this function.

        These counters are reset to zero at each successive call to
        `solve()`.
        """
        cdef PetscInt ival1 = 0
        cdef PetscInt ival2 = 0
        cdef PetscInt ival3 = 0
        CHKERR( EPSGetOperationCounters(self.eps, &ival1, &ival2, &ival3) )
        return (toInt(ival1), toInt(ival2), toInt(ival3))

    #

    def setPowerShiftType(self, shift):
        """
        Sets the type of shifts used during the power iteration. This
        can be used to emulate the Rayleigh Quotient Iteration (RQI)
        method.

        Parameters
        ----------
        shift: EPS.PowerShiftType enumerate
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

    def getPowerShiftType(self):
        """
        Gets the type of shifts used during the power iteration.

        Returns
        -------
        shift: EPS.PowerShiftType enumerate
               The type of shift.
        """
        cdef SlepcEPSPowerShiftType val = EPS_POWER_SHIFT_CONSTANT
        CHKERR( EPSPowerGetShiftType(self.eps, &val) )
        return val

    def setArnoldiDelayed(self, delayed):
        """
        Activates or deactivates delayed reorthogonalization in the
        Arnoldi iteration.

        Parameters
        ----------
        delayed: boolean
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
        cdef PetscTruth val = PETSC_FALSE
        if delayed: val = PETSC_TRUE
        CHKERR( EPSArnoldiSetDelayed(self.eps, val) )

    def getArnoldiDelayed(self):
        """
        Gets the type of reorthogonalization used during the Arnoldi
        iteration.

        Returns
        -------
        delayed: boolean
                 True if delayed reorthogonalization is to be used.
        """
        cdef PetscTruth val = PETSC_FALSE
        CHKERR( EPSArnoldiGetDelayed(self.eps, &val) )
        return val

    def setLanczosReorthogType(self, reorthog):
        """
        Sets the type of reorthogonalization used during the Lanczos
        iteration.

        Parameters
        ----------
        reorthog: EPS.LanczosReorthogType enumerate
                  The type of reorthogonalization.

        Notes
        -----
        This call is only relevant if the type was set to
        `EPS.Type.LANCZOS` with `setType()`.
        """
        cdef SlepcEPSLanczosReorthogType val = reorthog
        CHKERR( EPSLanczosSetReorthog(self.eps, val) )

    def getLanczosReorthogType(self):
        """
        Gets the type of reorthogonalization used during the Lanczos
        iteration.

        Returns
        -------
        reorthog: EPS.LanczosReorthogType enumerate
                  The type of reorthogonalization.
        """
        cdef SlepcEPSLanczosReorthogType val = \
            EPS_LANCZOS_REORTHOG_LOCAL
        CHKERR( EPSLanczosGetReorthog(self.eps, &val) )
        return val

    #
    property problem_type:
        def __get__(self):
            return self.getProblemType()
        def __set__(self, value):
            self.setProblemType(value)

    property extraction:
        def __get__(self):
            return self.getExtraction()
        def __set__(self, value):
            self.setExtraction(value)

    property which:
        def __get__(self):
            return self.getWhichEigenpairs()
        def __set__(self, value):
            self.setWhichEigenpairs(value)

    property target:
        def __get__(self):
            return self.getTarget()
        def __set__(self, value):
            self.setTarget(value)

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

    property st:
        def __get__(self):
            return self.getST()
        def __set__(self, value):
            self.setST(value)

    property ip:
        def __get__(self):
            return self.getIP()
        def __set__(self, value):
            self.setIP(value)

# -----------------------------------------------------------------------------

del EPSType
del EPSProblemType
del EPSExtraction
del EPSBalance
del EPSWhich
del EPSConvergedReason
del EPSPowerShiftType
del EPSLanczosReorthogType

# -----------------------------------------------------------------------------
