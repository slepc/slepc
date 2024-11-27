# -----------------------------------------------------------------------------

class SVDType(object):
    """
    SVD types

    - `CROSS`:      Eigenproblem with the cross-product matrix.
    - `CYCLIC`:     Eigenproblem with the cyclic matrix.
    - `LAPACK`:     Wrappers to dense SVD solvers in Lapack.
    - `LANCZOS`:    Lanczos.
    - `TRLANCZOS`:  Thick-restart Lanczos.
    - `RANDOMIZED`: Iterative RSVD for low-rank matrices.

    Wrappers to external SVD solvers
    (should be enabled during installation of SLEPc)

    - `SCALAPACK`:
    - `KSVD`:
    - `ELEMENTAL`:
    - `PRIMME`:
    """
    CROSS      = S_(SVDCROSS)
    CYCLIC     = S_(SVDCYCLIC)
    LAPACK     = S_(SVDLAPACK)
    LANCZOS    = S_(SVDLANCZOS)
    TRLANCZOS  = S_(SVDTRLANCZOS)
    RANDOMIZED = S_(SVDRANDOMIZED)
    SCALAPACK  = S_(SVDSCALAPACK)
    KSVD       = S_(SVDKSVD)
    ELEMENTAL  = S_(SVDELEMENTAL)
    PRIMME     = S_(SVDPRIMME)

class SVDProblemType(object):
    """
    SVD problem type

    - `STANDARD`:    Standard SVD.
    - `GENERALIZED`: Generalized singular value decomposition (GSVD).
    - `HYPERBOLIC` : Hyperbolic singular value decomposition (HSVD).
    """
    STANDARD    = SVD_STANDARD
    GENERALIZED = SVD_GENERALIZED
    HYPERBOLIC  = SVD_HYPERBOLIC

class SVDErrorType(object):
    """
    SVD error type to assess accuracy of computed solutions

    - `ABSOLUTE`: Absolute error.
    - `RELATIVE`: Relative error.
    - `NORM`:     Error relative to the matrix norm.
    """
    ABSOLUTE = SVD_ERROR_ABSOLUTE
    RELATIVE = SVD_ERROR_RELATIVE
    NORM     = SVD_ERROR_NORM

class SVDWhich(object):
    """
    SVD desired part of spectrum

    - `LARGEST`:  Largest singular values.
    - `SMALLEST`: Smallest singular values.
    """
    LARGEST  = SVD_LARGEST
    SMALLEST = SVD_SMALLEST

class SVDConv(object):
    """
    SVD convergence test

    - `ABS`:   Absolute convergence test.
    - `REL`:   Convergence test relative to the singular value.
    - `NORM`:  Convergence test relative to the matrix norms.
    - `MAXIT`: No convergence until maximum number of iterations has been reached.
    - `USER`:  User-defined convergence test.
    """
    ABS   = SVD_CONV_ABS
    REL   = SVD_CONV_REL
    NORM  = SVD_CONV_NORM
    MAXIT = SVD_CONV_MAXIT
    USER  = SVD_CONV_USER

class SVDStop(object):
    """
    SVD stopping test

    - `BASIC`:     Default stopping test.
    - `USER`:      User-defined stopping test.
    - `THRESHOLD`: Threshold stopping test.
    """
    BASIC     = SVD_STOP_BASIC
    USER      = SVD_STOP_USER
    THRESHOLD = SVD_STOP_THRESHOLD

class SVDConvergedReason(object):
    """
    SVD convergence reasons

    - `CONVERGED_TOL`:          All eigenpairs converged to requested tolerance.
    - `CONVERGED_USER`:         User-defined convergence criterion satisfied.
    - `CONVERGED_MAXIT`:        Maximum iterations completed in case MAXIT convergence criterion.
    - `DIVERGED_ITS`:           Maximum number of iterations exceeded.
    - `DIVERGED_BREAKDOWN`:     Solver failed due to breakdown.
    - `DIVERGED_SYMMETRY_LOST`: Underlying indefinite eigensolver was not able to keep symmetry.
    - `CONVERGED_ITERATING`:    Iteration not finished yet.
    """
    CONVERGED_TOL          = SVD_CONVERGED_TOL
    CONVERGED_USER         = SVD_CONVERGED_USER
    CONVERGED_MAXIT        = SVD_CONVERGED_MAXIT
    DIVERGED_ITS           = SVD_DIVERGED_ITS
    DIVERGED_BREAKDOWN     = SVD_DIVERGED_BREAKDOWN
    DIVERGED_SYMMETRY_LOST = SVD_DIVERGED_SYMMETRY_LOST
    CONVERGED_ITERATING    = SVD_CONVERGED_ITERATING
    ITERATING              = SVD_CONVERGED_ITERATING

class SVDTRLanczosGBidiag(object):
    """
    SVD TRLanczos bidiagonalization choices for the GSVD case

    - `SINGLE`: Single bidiagonalization (Qa).
    - `UPPER`:  Joint bidiagonalization, both Qa and Qb in upper bidiagonal form.
    - `LOWER`:  Joint bidiagonalization, Qa lower bidiagonal, Qb upper bidiagonal.
    """
    SINGLE = SVD_TRLANCZOS_GBIDIAG_SINGLE
    UPPER  = SVD_TRLANCZOS_GBIDIAG_UPPER
    LOWER  = SVD_TRLANCZOS_GBIDIAG_LOWER

# -----------------------------------------------------------------------------

cdef class SVD(Object):

    """
    SVD
    """

    Type            = SVDType
    ProblemType     = SVDProblemType
    ErrorType       = SVDErrorType
    Which           = SVDWhich
    Conv            = SVDConv
    Stop            = SVDStop
    ConvergedReason = SVDConvergedReason

    TRLanczosGBidiag = SVDTRLanczosGBidiag

    def __cinit__(self):
        self.obj = <PetscObject*> &self.svd
        self.svd = NULL

    def view(self, Viewer viewer=None):
        """
        Prints the SVD data structure.

        Parameters
        ----------
        viewer: Viewer, optional
                Visualization context; if not provided, the standard
                output is used.
        """
        cdef PetscViewer vwr = def_Viewer(viewer)
        CHKERR( SVDView(self.svd, vwr) )

    def destroy(self):
        """
        Destroys the SVD object.
        """
        CHKERR( SVDDestroy(&self.svd) )
        self.svd = NULL
        return self

    def reset(self):
        """
        Resets the SVD object.
        """
        CHKERR( SVDReset(self.svd) )

    def create(self, comm=None):
        """
        Creates the SVD object.

        Parameters
        ----------
        comm: Comm, optional
              MPI communicator; if not provided, it defaults to all
              processes.
        """
        cdef MPI_Comm ccomm = def_Comm(comm, SLEPC_COMM_DEFAULT())
        cdef SlepcSVD newsvd = NULL
        CHKERR( SVDCreate(ccomm, &newsvd) )
        CHKERR( SlepcCLEAR(self.obj) ); self.svd = newsvd
        return self

    def setType(self, svd_type):
        """
        Selects the particular solver to be used in the SVD object.

        Parameters
        ----------
        svd_type: `SVD.Type` enumerate
                  The solver to be used.

        Notes
        -----
        See `SVD.Type` for available methods. The default is CROSS.
        Normally, it is best to use `setFromOptions()` and then set
        the SVD type from the options database rather than by using
        this routine.  Using the options database provides the user
        with maximum flexibility in evaluating the different available
        methods.
        """
        cdef SlepcSVDType cval = NULL
        svd_type = str2bytes(svd_type, &cval)
        CHKERR( SVDSetType(self.svd, cval) )

    def getType(self):
        """
        Gets the SVD type of this object.

        Returns
        -------
        type: `SVD.Type` enumerate
              The solver currently being used.
        """
        cdef SlepcSVDType svd_type = NULL
        CHKERR( SVDGetType(self.svd, &svd_type) )
        return bytes2str(svd_type)

    def getOptionsPrefix(self):
        """
        Gets the prefix used for searching for all SVD options in the
        database.

        Returns
        -------
        prefix: string
                The prefix string set for this SVD object.
        """
        cdef const char *prefix = NULL
        CHKERR( SVDGetOptionsPrefix(self.svd, &prefix) )
        return bytes2str(prefix)

    def setOptionsPrefix(self, prefix):
        """
        Sets the prefix used for searching for all SVD options in the
        database.

        Parameters
        ----------
        prefix: string
                The prefix string to prepend to all SVD option
                requests.

        Notes
        -----
        A hyphen (-) must NOT be given at the beginning of the prefix
        name.  The first character of all runtime options is
        AUTOMATICALLY the hyphen.

        For example, to distinguish between the runtime options for
        two different SVD contexts, one could call::

            S1.setOptionsPrefix("svd1_")
            S2.setOptionsPrefix("svd2_")
        """
        cdef const char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( SVDSetOptionsPrefix(self.svd, cval) )

    def appendOptionsPrefix(self, prefix):
        """
        Appends to the prefix used for searching for all SVD options
        in the database.

        Parameters
        ----------
        prefix: string
                The prefix string to prepend to all SVD option requests.
        """
        cdef const char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( SVDAppendOptionsPrefix(self.svd, cval) )

    def setFromOptions(self):
        """
        Sets SVD options from the options database. This routine must
        be called before `setUp()` if the user is to be allowed to set
        the solver type.

        Notes
        -----
        To see all options, run your program with the ``-help``
        option.
        """
        CHKERR( SVDSetFromOptions(self.svd) )

    def getProblemType(self):
        """
        Gets the problem type from the SVD object.

        Returns
        -------
        problem_type: `SVD.ProblemType` enumerate
                      The problem type that was previously set.
        """
        cdef SlepcSVDProblemType val = SVD_STANDARD
        CHKERR( SVDGetProblemType(self.svd, &val) )
        return val

    def setProblemType(self, problem_type):
        """
        Specifies the type of the singular value problem.

        Parameters
        ----------
        problem_type: `SVD.ProblemType` enumerate
               The problem type to be set.
        """
        cdef SlepcSVDProblemType val = problem_type
        CHKERR( SVDSetProblemType(self.svd, val) )

    def isGeneralized(self):
        """
        Tells whether the SVD object corresponds to a generalized
        singular value problem.

        Returns
        -------
        flag: bool
              True if two matrices were set with `setOperators()`.
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( SVDIsGeneralized(self.svd, &tval) )
        return toBool(tval)

    def isHyperbolic(self):
        """
        Tells whether the SVD object corresponds to a hyperbolic
        singular value problem.

        Returns
        -------
        flag: bool
              True if the problem was specified as hyperbolic.
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( SVDIsHyperbolic(self.svd, &tval) )
        return toBool(tval)

    #

    def getImplicitTranspose(self):
        """
        Gets the mode used to handle the transpose of the matrix
        associated with the singular value problem.

        Returns
        -------
        impl: bool
              How to handle the transpose (implicitly or not).
        """
        cdef PetscBool val = PETSC_FALSE
        CHKERR( SVDGetImplicitTranspose(self.svd, &val) )
        return toBool(val)

    def setImplicitTranspose(self, mode):
        """
        Indicates how to handle the transpose of the matrix
        associated with the singular value problem.

        Parameters
        ----------
        impl: bool
              How to handle the transpose (implicitly or not).

        Notes
        -----
        By default, the transpose of the matrix is explicitly built
        (if the matrix has defined the MatTranspose operation).

        If this flag is set to true, the solver does not build the
        transpose, but handles it implicitly via MatMultTranspose().
        """
        cdef PetscBool val = asBool(mode)
        CHKERR( SVDSetImplicitTranspose(self.svd, val) )

    def getWhichSingularTriplets(self):
        """
        Returns which singular triplets are to be sought.

        Returns
        -------
        which: `SVD.Which` enumerate
               The singular values to be sought (either largest or
               smallest).
        """
        cdef SlepcSVDWhich val = SVD_LARGEST
        CHKERR( SVDGetWhichSingularTriplets(self.svd, &val) )
        return val

    def setWhichSingularTriplets(self, which):
        """
        Specifies which singular triplets are to be sought.

        Parameters
        ----------
        which: `SVD.Which` enumerate
               The singular values to be sought (either largest or
               smallest).
        """
        cdef SlepcSVDWhich val = which
        CHKERR( SVDSetWhichSingularTriplets(self.svd, val) )

    def getThreshold(self):
        """
        Gets the threshold used in the threshold stopping test.

        Returns
        -------
        thres: float
             The threshold.
        rel: bool
             Whether the threshold is relative or not.
        """
        cdef PetscReal rval = 0
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( SVDGetThreshold(self.svd, &rval, &tval) )
        return (toReal(rval), toBool(tval))

    def setThreshold(self, thres, rel=False):
        """
        Sets the threshold used in the threshold stopping test.

        Parameters
        ----------
        thres: float
             The threshold.
        rel: bool, optional
             Whether the threshold is relative or not.

        Notes
        -----
        This function internally sets a special stopping test based on
        the threshold, where singular values are computed in sequence
        until one of the computed singular values is below/above the
        threshold (depending on whether largest or smallest singular
        values are computed).

        In the case of largest singular values, the threshold can be
        made relative with respect to the largest singular value
        (i.e., the matrix norm).
        """
        cdef PetscReal rval = asReal(thres)
        cdef PetscBool tval = asBool(rel)
        CHKERR( SVDSetThreshold(self.svd, rval, tval) )

    def getTolerances(self):
        """
        Gets the tolerance and maximum iteration count used by the
        default SVD convergence tests.

        Returns
        -------
        tol: float
             The convergence tolerance.
        max_it: int
             The maximum number of iterations
        """
        cdef PetscReal rval = 0
        cdef PetscInt  ival = 0
        CHKERR( SVDGetTolerances(self.svd, &rval, &ival) )
        return (toReal(rval), toInt(ival))

    def setTolerances(self, tol=None, max_it=None):
        """
        Sets the tolerance and maximum iteration count used by the
        default SVD convergence tests.

        Parameters
        ----------
        tol: float, optional
             The convergence tolerance.
        max_it: int, optional
             The maximum number of iterations

        Notes
        -----
        Use `DECIDE` for `max_it` to assign a reasonably good value,
        which is dependent on the solution method.
        """
        cdef PetscReal rval = PETSC_DEFAULT
        cdef PetscInt  ival = PETSC_DEFAULT
        if tol    is not None: rval = asReal(tol)
        if max_it is not None: ival = asInt(max_it)
        CHKERR( SVDSetTolerances(self.svd, rval, ival) )

    def getConvergenceTest(self):
        """
        Return the method used to compute the error estimate
        used in the convergence test.

        Returns
        -------
        conv: SVD.Conv
            The method used to compute the error estimate
            used in the convergence test.
        """
        cdef SlepcSVDConv conv = SVD_CONV_REL
        CHKERR( SVDGetConvergenceTest(self.svd, &conv) )
        return conv

    def setConvergenceTest(self, conv):
        """
        Specifies how to compute the error estimate
        used in the convergence test.

        Parameters
        ----------
        conv: SVD.Conv
            The method used to compute the error estimate
            used in the convergence test.
        """
        cdef SlepcSVDConv tconv = conv
        CHKERR( SVDSetConvergenceTest(self.svd, tconv) )

    def getTrackAll(self):
        """
        Returns the flag indicating whether all residual norms must be
        computed or not.

        Returns
        -------
        trackall: bool
            Whether the solver compute all residuals or not.
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( SVDGetTrackAll(self.svd, &tval) )
        return toBool(tval)

    def setTrackAll(self, trackall):
        """
        Specifies if the solver must compute the residual of all
        approximate singular triplets or not.

        Parameters
        ----------
        trackall: bool
            Whether compute all residuals or not.
        """
        cdef PetscBool tval = asBool(trackall)
        CHKERR( SVDSetTrackAll(self.svd, tval) )

    def getDimensions(self):
        """
        Gets the number of singular values to compute and the
        dimension of the subspace.

        Returns
        -------
        nsv: int
             Number of singular values to compute.
        ncv: int
             Maximum dimension of the subspace to be used by the
             solver.
        mpd: int
             Maximum dimension allowed for the projected problem.
        """
        cdef PetscInt ival1 = 0
        cdef PetscInt ival2 = 0
        cdef PetscInt ival3 = 0
        CHKERR( SVDGetDimensions(self.svd, &ival1, &ival2, &ival3) )
        return (toInt(ival1), toInt(ival2), toInt(ival3))

    def setDimensions(self, nsv=None, ncv=None, mpd=None):
        """
        Sets the number of singular values to compute and the
        dimension of the subspace.

        Parameters
        ----------
        nsv: int, optional
             Number of singular values to compute.
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

         - In cases where `nsv` is small, the user sets `ncv`
           (a reasonable default is 2 * `nsv`).
         - In cases where `nsv` is large, the user sets `mpd`.

        The value of `ncv` should always be between `nsv` and (`nsv` +
        `mpd`), typically `ncv` = `nsv` + `mpd`. If `nsv` is not too
        large, `mpd` = `nsv` is a reasonable choice, otherwise a
        smaller value should be used.
        """
        cdef PetscInt ival1 = PETSC_DEFAULT
        cdef PetscInt ival2 = PETSC_DEFAULT
        cdef PetscInt ival3 = PETSC_DEFAULT
        if nsv is not None: ival1 = asInt(nsv)
        if ncv is not None: ival2 = asInt(ncv)
        if mpd is not None: ival3 = asInt(mpd)
        CHKERR( SVDSetDimensions(self.svd, ival1, ival2, ival3) )

    def getBV(self):
        """
        Obtain the basis vectors objects associated to the SVD object.

        Returns
        -------
        V: BV
            The basis vectors context for right singular vectors.
        U: BV
            The basis vectors context for left singular vectors.
        """
        cdef BV V = BV()
        cdef BV U = BV()
        CHKERR( SVDGetBV(self.svd, &V.bv, &U.bv) )
        CHKERR( PetscINCREF(V.obj) )
        CHKERR( PetscINCREF(U.obj) )
        return (V,U)

    def setBV(self, BV V,BV U=None):
        """
        Associates basis vectors objects to the SVD solver.

        Parameters
        ----------
        V: BV
            The basis vectors context for right singular vectors.
        U: BV
            The basis vectors context for left singular vectors.
        """
        cdef SlepcBV VBV = V.bv
        cdef SlepcBV UBV = U.bv if U is not None else <SlepcBV>NULL
        CHKERR( SVDSetBV(self.svd, VBV, UBV) )

    def getDS(self):
        """
        Obtain the direct solver associated to the singular value solver.

        Returns
        -------
        ds: DS
            The direct solver context.
        """
        cdef DS ds = DS()
        CHKERR( SVDGetDS(self.svd, &ds.ds) )
        CHKERR( PetscINCREF(ds.obj) )
        return ds

    def setDS(self, DS ds):
        """
        Associates a direct solver object to the singular value solver.

        Parameters
        ----------
        ds: DS
            The direct solver context.
        """
        CHKERR( SVDSetDS(self.svd, ds.ds) )

    def getOperators(self):
        """
        Gets the matrices associated with the singular value problem.

        Returns
        -------
        A: Mat
           The matrix associated with the singular value problem.
        B: Mat
           The second matrix in the case of GSVD.
        """
        cdef Mat A = Mat()
        cdef Mat B = Mat()
        CHKERR( SVDGetOperators(self.svd, &A.mat, &B.mat) )
        CHKERR( PetscINCREF(A.obj) )
        CHKERR( PetscINCREF(B.obj) )
        return (A, B)

    def setOperators(self, Mat A, Mat B=None):
        """
        Sets the matrices associated with the singular value problem.

        Parameters
        ----------
        A: Mat
           The matrix associated with the singular value problem.
        B: Mat, optional
           The second matrix in the case of GSVD; if not provided,
           a usual SVD is assumed.
        """
        cdef PetscMat Bmat = B.mat if B is not None else <PetscMat>NULL
        CHKERR( SVDSetOperators(self.svd, A.mat, Bmat) )

    def getSignature(self, Vec omega=None):
        """
        Gets the signature matrix defining a hyperbolic singular value problem.

        Parameters
        ----------
        omega: Vec
           Optional vector to store the diagonal elements of the signature matrix.

        Returns
        -------
        omega: Vec
           A vector containing the diagonal elements of the signature matrix.
        """
        cdef PetscMat A = NULL
        if omega is None:
            omega = Vec()
        if omega.vec == NULL:
            CHKERR( SVDGetOperators(self.svd, &A, <PetscMat*>NULL) )
            CHKERR( MatCreateVecs(A, <PetscVec*>NULL, &omega.vec) )
        CHKERR( SVDGetSignature(self.svd, omega.vec) )
        return omega

    def setSignature(self, Vec omega=None):
        """
        Sets the signature matrix defining a hyperbolic singular value problem.

        Parameters
        ----------
        omega: Vec, optional
           A vector containing the diagonal elements of the signature matrix.
        """
        cdef PetscVec Ovec = omega.vec if omega is not None else <PetscVec>NULL
        CHKERR( SVDSetSignature(self.svd, Ovec) )

    #

    def setInitialSpaces(self, spaceright=None, spaceleft=None):
        """
        Sets the initial spaces from which the SVD solver starts to
        iterate.

        Parameters
        ----------
        spaceright: sequence of Vec
           The right initial space.
        spaceleft: sequence of Vec
           The left initial space.
        """
        cdef Py_ssize_t i = 0
        if spaceright is None: spaceright = []
        elif isinstance(spaceright, Vec): spaceright = [spaceright]
        cdef PetscVec *isr = NULL
        cdef Py_ssize_t nr = len(spaceright)
        cdef tmp1 = allocate(<size_t>nr*sizeof(PetscVec),<void**>&isr)
        for i in range(nr): isr[i] = (<Vec?>spaceright[i]).vec
        if spaceleft is None: spaceright = []
        elif isinstance(spaceleft, Vec): spaceleft = [spaceleft]
        cdef PetscVec *isl = NULL
        cdef Py_ssize_t nl = len(spaceleft)
        cdef tmp2 = allocate(<size_t>nl*sizeof(PetscVec),<void**>&isl)
        for i in range(nl): isl[i] = (<Vec?>spaceleft[i]).vec
        CHKERR( SVDSetInitialSpaces(self.svd, <PetscInt>nr, isr, <PetscInt>nl, isl) )

    #

    def setStoppingTest(self, stopping, args=None, kargs=None):
        """
        Sets a function to decide when to stop the outer iteration of the eigensolver.
        """
        if stopping is not None:
            if args is None: args = ()
            if kargs is None: kargs = {}
            self.set_attr('__stopping__', (stopping, args, kargs))
            CHKERR( SVDSetStoppingTestFunction(self.svd, SVD_Stopping, NULL, NULL) )
        else:
            self.set_attr('__stopping__', None)
            CHKERR( SVDSetStoppingTestFunction(self.svd, SVDStoppingBasic, NULL, NULL) )

    def getStoppingTest(self):
        """
        Gets the stopping function.
        """
        return self.get_attr('__stopping__')

    #

    def setMonitor(self, monitor, args=None, kargs=None):
        """
        Appends a monitor function to the list of monitors.
        """
        if monitor is None: return
        cdef object monitorlist = self.get_attr('__monitor__')
        if monitorlist is None:
            monitorlist = []
            self.set_attr('__monitor__', monitorlist)
            CHKERR( SVDMonitorSet(self.svd, SVD_Monitor, NULL, NULL) )
        if args is None: args = ()
        if kargs is None: kargs = {}
        monitorlist.append((monitor, args, kargs))

    def getMonitor(self):
        """
        Gets the list of monitor functions.
        """
        return self.get_attr('__monitor__')

    def cancelMonitor(self):
        """
        Clears all monitors for an `SVD` object.
        """
        CHKERR( SVDMonitorCancel(self.svd) )
        self.set_attr('__monitor__', None)

    #

    def setUp(self):
        """
        Sets up all the internal data structures necessary for the
        execution of the singular value solver.

        Notes
        -----
        This function need not be called explicitly in most cases,
        since `solve()` calls it. It can be useful when one wants to
        measure the set-up time separately from the solve time.
        """
        CHKERR( SVDSetUp(self.svd) )

    def solve(self):
        """
        Solves the singular value problem.
        """
        CHKERR( SVDSolve(self.svd) )

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
        CHKERR( SVDGetIterationNumber(self.svd, &ival) )
        return toInt(ival)

    def getConvergedReason(self):
        """
        Gets the reason why the `solve()` iteration was stopped.

        Returns
        -------
        reason: `SVD.ConvergedReason` enumerate
                Negative value indicates diverged, positive value
                converged.
        """
        cdef SlepcSVDConvergedReason val = SVD_CONVERGED_ITERATING
        CHKERR( SVDGetConvergedReason(self.svd, &val) )
        return val

    def getConverged(self):
        """
        Gets the number of converged singular triplets.

        Returns
        -------
        nconv: int
               Number of converged singular triplets.

        Notes
        -----
        This function should be called after `solve()` has finished.
        """
        cdef PetscInt ival = 0
        CHKERR( SVDGetConverged(self.svd, &ival) )
        return toInt(ival)

    def getValue(self, int i):
        """
        Gets the i-th singular value as computed by `solve()`.

        Parameters
        ----------
        i: int
           Index of the solution to be obtained.

        Returns
        -------
        s: float
           The computed singular value.

        Notes
        -----
        The index ``i`` should be a value between ``0`` and
        ``nconv-1`` (see `getConverged()`. Singular triplets are
        indexed according to the ordering criterion established with
        `setWhichSingularTriplets()`.
        """
        cdef PetscReal rval = 0
        CHKERR( SVDGetSingularTriplet(self.svd, i, &rval, NULL, NULL) )
        return toReal(rval)

    def getVectors(self, int i, Vec U, Vec V):
        """
        Gets the i-th left and right singular vectors as computed by
        `solve()`.

        Parameters
        ----------
        i: int
           Index of the solution to be obtained.
        U: Vec
           Placeholder for the returned left singular vector.
        V: Vec
           Placeholder for the returned right singular vector.

        Notes
        -----
        The index ``i`` should be a value between ``0`` and
        ``nconv-1`` (see `getConverged()`. Singular triplets are
        indexed according to the ordering criterion established with
        `setWhichSingularTriplets()`.
        """
        cdef PetscReal dummy = 0
        CHKERR( SVDGetSingularTriplet(self.svd, i, &dummy, U.vec, V.vec) )

    def getSingularTriplet(self, int i, Vec U=None, Vec V=None):
        """
        Gets the i-th triplet of the singular value decomposition as
        computed by `solve()`. The solution consists of the singular
        value and its left and right singular vectors.

        Parameters
        ----------
        i: int
           Index of the solution to be obtained.
        U: Vec
           Placeholder for the returned left singular vector.
        V: Vec
           Placeholder for the returned right singular vector.

        Returns
        -------
        s: float
           The computed singular value.

        Notes
        -----
        The index ``i`` should be a value between ``0`` and
        ``nconv-1`` (see `getConverged()`. Singular triplets are
        indexed according to the ordering criterion established with
        `setWhichSingularTriplets()`.
        """
        cdef PetscReal rval = 0
        cdef PetscVec Uvec = U.vec if U is not None else <PetscVec>NULL
        cdef PetscVec Vvec = V.vec if V is not None else <PetscVec>NULL
        CHKERR( SVDGetSingularTriplet(self.svd, i, &rval, Uvec, Vvec) )
        return toReal(rval)

    #

    def computeError(self, int i, etype=None):
        """
        Computes the error (based on the residual norm) associated with the i-th
        singular triplet.

        Parameters
        ----------
        i: int
           Index of the solution to be considered.
        etype: `SVD.ErrorType` enumerate
           The error type to compute.

        Returns
        -------
        e: real
           The relative error bound, computed in various ways from the residual norm
           ``sqrt(n1^2+n2^2)`` where ``n1 = ||A*v-sigma*u||_2``,
           ``n2 = ||A^T*u-sigma*v||_2``, ``sigma`` is the singular
           value, ``u`` and ``v`` are the left and right singular
           vectors.

        Notes
        -----
        The index ``i`` should be a value between ``0`` and
        ``nconv-1`` (see `getConverged()`).
        """
        cdef SlepcSVDErrorType et = SVD_ERROR_RELATIVE
        cdef PetscReal rval = 0
        if etype is not None: et = etype
        CHKERR( SVDComputeError(self.svd, i, et, &rval) )
        return toReal(rval)

    def errorView(self, etype=None, Viewer viewer=None):
        """
        Displays the errors associated with the computed solution
        (as well as the eigenvalues).

        Parameters
        ----------
        etype: `SVD.ErrorType` enumerate, optional
           The error type to compute.
        viewer: Viewer, optional.
                Visualization context; if not provided, the standard
                output is used.

        Notes
        -----
        By default, this function checks the error of all eigenpairs and prints
        the eigenvalues if all of them are below the requested tolerance.
        If the viewer has format ``ASCII_INFO_DETAIL`` then a table with
        eigenvalues and corresponding errors is printed.

        """
        cdef SlepcSVDErrorType et = SVD_ERROR_RELATIVE
        if etype is not None: et = etype
        cdef PetscViewer vwr = def_Viewer(viewer)
        CHKERR( SVDErrorView(self.svd, et, vwr) )

    def valuesView(self, Viewer viewer=None):
        """
        Displays the computed singular values in a viewer.

        Parameters
        ----------
        viewer: Viewer, optional.
                Visualization context; if not provided, the standard
                output is used.
        """
        cdef PetscViewer vwr = def_Viewer(viewer)
        CHKERR( SVDValuesView(self.svd, vwr) )

    def vectorsView(self, Viewer viewer=None):
        """
        Outputs computed singular vectors to a viewer.

        Parameters
        ----------
        viewer: Viewer, optional.
                Visualization context; if not provided, the standard
                output is used.
        """
        cdef PetscViewer vwr = def_Viewer(viewer)
        CHKERR( SVDVectorsView(self.svd, vwr) )

    #

    def setCrossEPS(self, EPS eps):
        """
        Associate an eigensolver object (`EPS`) to the singular value
        solver.

        Parameters
        ----------
        eps: EPS
             The eigensolver object.
        """
        CHKERR( SVDCrossSetEPS(self.svd, eps.eps) )

    def getCrossEPS(self):
        """
        Retrieve the eigensolver object (`EPS`) associated to the
        singular value solver.

        Returns
        -------
        eps: EPS
             The eigensolver object.
        """
        cdef EPS eps = EPS()
        CHKERR( SVDCrossGetEPS(self.svd, &eps.eps) )
        CHKERR( PetscINCREF(eps.obj) )
        return eps

    def setCrossExplicitMatrix(self, flag=True):
        """
        Indicate if the eigensolver operator ``A^T*A`` must be
        computed explicitly.

        Parameters
        ----------
        flag: bool
              True if ``A^T*A`` is built explicitly.
        """
        cdef PetscBool tval = asBool(flag)
        CHKERR( SVDCrossSetExplicitMatrix(self.svd, tval) )

    def getCrossExplicitMatrix(self):
        """
        Returns the flag indicating if ``A^T*A`` is built explicitly.

        Returns
        -------
        flag: bool
              True if ``A^T*A`` is built explicitly.
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( SVDCrossGetExplicitMatrix(self.svd, &tval) )
        return toBool(tval)

    def setCyclicEPS(self, EPS eps):
        """
        Associate an eigensolver object (`EPS`) to the singular value
        solver.

        Parameters
        ----------
        eps: EPS
             The eigensolver object.
        """
        CHKERR( SVDCyclicSetEPS(self.svd, eps.eps) )

    def getCyclicEPS(self):
        """
        Retrieve the eigensolver object (`EPS`) associated to the
        singular value solver.

        Returns
        -------
        eps: EPS
             The eigensolver object.
        """
        cdef EPS eps = EPS()
        CHKERR( SVDCyclicGetEPS(self.svd, &eps.eps) )
        CHKERR( PetscINCREF(eps.obj) )
        return eps

    def setCyclicExplicitMatrix(self, flag=True):
        """
        Indicate if the eigensolver operator ``H(A) = [ 0 A ; A^T 0
        ]`` must be computed explicitly.

        Parameters
        ----------
        flag: bool
              True if ``H(A)`` is built explicitly.
        """
        cdef PetscBool tval = asBool(flag)
        CHKERR( SVDCyclicSetExplicitMatrix(self.svd, tval) )

    def getCyclicExplicitMatrix(self):
        """
        Returns the flag indicating if ``H(A) = [ 0 A ; A^T 0 ]`` is
        built explicitly.

        Returns
        -------
        flag: bool
              True if ``H(A)`` is built explicitly.
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( SVDCyclicGetExplicitMatrix(self.svd, &tval) )
        return toBool(tval)

    def setLanczosOneSide(self, flag=True):
        """
        Indicate if the variant of the Lanczos method to be used is
        one-sided or two-sided.

        Parameters
        ----------
        flag: bool
              True if the method is one-sided.

        Notes
        -----
        By default, a two-sided variant is selected, which is
        sometimes slightly more robust. However, the one-sided variant
        is faster because it avoids the orthogonalization associated
        to left singular vectors. It also saves the memory required
        for storing such vectors.
        """
        cdef PetscBool tval = asBool(flag)
        CHKERR( SVDLanczosSetOneSide(self.svd, tval) )

    def getLanczosOneSide(self):
        """
        Gets if the variant of the Lanczos method to be used is
        one-sided or two-sided.

        Returns
        -------
        delayed: bool
                 True if the method is one-sided.
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( SVDLanczosGetOneSide(self.svd, &tval) )
        return toBool(tval)

    def setTRLanczosOneSide(self, flag=True):
        """
        Indicate if the variant of the thick-restart Lanczos method to
        be used is one-sided or two-sided.

        Parameters
        ----------
        flag: bool
              True if the method is one-sided.

        Notes
        -----
        By default, a two-sided variant is selected, which is
        sometimes slightly more robust. However, the one-sided variant
        is faster because it avoids the orthogonalization associated
        to left singular vectors.
        """
        cdef PetscBool tval = asBool(flag)
        CHKERR( SVDLanczosSetOneSide(self.svd, tval) )

    def getTRLanczosOneSide(self):
        """
        Gets if the variant of the thick-restart Lanczos method to be
        used is one-sided or two-sided.

        Returns
        -------
        delayed: bool
                 True if the method is one-sided.
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( SVDTRLanczosGetOneSide(self.svd, &tval) )
        return toBool(tval)

    def setTRLanczosGBidiag(self, bidiag):
        """
        Sets the bidiagonalization choice to use in the GSVD
        TRLanczos solver.

        Parameters
        ----------
        bidiag: `SVD.TRLanczosGBidiag` enumerate
               The bidiagonalization choice.
        """
        cdef SlepcSVDTRLanczosGBidiag val = bidiag
        CHKERR( SVDTRLanczosSetGBidiag(self.svd, val) )

    def getTRLanczosGBidiag(self):
        """
        Returns bidiagonalization choice used in the GSVD
        TRLanczos solver.

        Returns
        -------
        bidiag: `SVD.TRLanczosGBidiag` enumerate
               The bidiagonalization choice.
        """
        cdef SlepcSVDTRLanczosGBidiag val = SVD_TRLANCZOS_GBIDIAG_LOWER
        CHKERR( SVDTRLanczosGetGBidiag(self.svd, &val) )
        return val

    def setTRLanczosRestart(self, keep):
        """
        Sets the restart parameter for the thick-restart Lanczos method, in
        particular the proportion of basis vectors that must be kept
        after restart.

        Parameters
        ----------
        keep: float
              The number of vectors to be kept at restart.

        Notes
        -----
        Allowed values are in the range [0.1,0.9]. The default is 0.5.
        """
        cdef PetscReal val = asReal(keep)
        CHKERR( SVDTRLanczosSetRestart(self.svd, val) )

    def getTRLanczosRestart(self):
        """
        Gets the restart parameter used in the thick-restart Lanczos method.

        Returns
        -------
        keep: float
              The number of vectors to be kept at restart.
        """
        cdef PetscReal val = 0
        CHKERR( SVDTRLanczosGetRestart(self.svd, &val) )
        return toReal(val)

    def setTRLanczosLocking(self, lock):
        """
        Choose between locking and non-locking variants of the
        thick-restart Lanczos method.

        Parameters
        ----------
        lock: bool
              True if the locking variant must be selected.

        Notes
        -----
        The default is to lock converged singular triplets when the method restarts.
        This behaviour can be changed so that all directions are kept in the
        working subspace even if already converged to working accuracy (the
        non-locking variant).
        """
        cdef PetscBool val = asBool(lock)
        CHKERR( SVDTRLanczosSetLocking(self.svd, val) )

    def getTRLanczosLocking(self):
        """
        Gets the locking flag used in the thick-restart Lanczos method.

        Returns
        -------
        lock: bool
              The locking flag.
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( SVDTRLanczosGetLocking(self.svd, &tval) )
        return toBool(tval)

    def setTRLanczosKSP(self, KSP ksp):
        """
        Associate a linear solver object to the SVD solver.

        Parameters
        ----------
        ksp: `KSP`
             The linear solver object.
        """
        CHKERR( SVDTRLanczosSetKSP(self.svd, ksp.ksp) )

    def getTRLanczosKSP(self):
        """
        Retrieve the linear solver object associated with the SVD solver.

        Returns
        -------
        ksp: `KSP`
             The linear solver object.
        """
        cdef KSP ksp = KSP()
        CHKERR( SVDTRLanczosGetKSP(self.svd, &ksp.ksp) )
        CHKERR( PetscINCREF(ksp.obj) )
        return ksp

    def setTRLanczosExplicitMatrix(self, flag=True):
        """
        Indicate if the matrix ``Z=[A;B]`` must be built explicitly.

        Parameters
        ----------
        flag: bool
              True if ``Z=[A;B]`` is built explicitly.
        """
        cdef PetscBool tval = asBool(flag)
        CHKERR( SVDTRLanczosSetExplicitMatrix(self.svd, tval) )

    def getTRLanczosExplicitMatrix(self):
        """
        Returns the flag indicating if ``Z=[A;B]`` is built explicitly.

        Returns
        -------
        flag: bool
              True if ``Z=[A;B]`` is built explicitly.
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( SVDTRLanczosGetExplicitMatrix(self.svd, &tval) )
        return toBool(tval)

    setOperator = setOperators  # backward compatibility

    #

    property problem_type:
        def __get__(self):
            return self.getProblemType()
        def __set__(self, value):
            self.setProblemType(value)

    property transpose_mode:
        def __get__(self):
            return self.getTransposeMode()
        def __set__(self, value):
            self.setTransposeMode(value)

    property which:
        def __get__(self):
            return self.getWhichSingularTriplets()
        def __set__(self, value):
            self.setWhichSingularTriplets(value)

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

    property track_all:
        def __get__(self):
            return self.getTrackAll()
        def __set__(self, value):
            self.setTrackAll(value)

    property ds:
        def __get__(self):
            return self.getDS()
        def __set__(self, value):
            self.setDS(value)

# -----------------------------------------------------------------------------

del SVDType
del SVDProblemType
del SVDErrorType
del SVDWhich
del SVDConv
del SVDStop
del SVDConvergedReason
del SVDTRLanczosGBidiag

# -----------------------------------------------------------------------------
