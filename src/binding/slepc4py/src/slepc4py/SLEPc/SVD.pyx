# -----------------------------------------------------------------------------

class SVDType(object):
    """
    SVD type.

    Native singular value solvers.

    - `CROSS`:      Eigenproblem with the cross-product matrix.
    - `CYCLIC`:     Eigenproblem with the cyclic matrix.
    - `LANCZOS`:    Explicitly restarted Lanczos.
    - `TRLANCZOS`:  Thick-restart Lanczos.
    - `RANDOMIZED`: Iterative RSVD for low-rank matrices.

    Wrappers to external SVD solvers
    (should be enabled during installation of SLEPc).

    - `LAPACK`:     Sequential dense SVD solver.
    - `SCALAPACK`:  Parallel dense SVD solver.
    - `KSVD`:       Parallel dense SVD solver.
    - `ELEMENTAL`:  Parallel dense SVD solver.
    - `PRIMME`:     Iterative SVD solvers of Davidson type.

    See Also
    --------
    slepc.SVDType
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
    SVD problem type.

    - `STANDARD`:    Standard SVD.
    - `GENERALIZED`: Generalized singular value decomposition (GSVD).
    - `HYPERBOLIC` : Hyperbolic singular value decomposition (HSVD).

    See Also
    --------
    slepc.SVDProblemType
    """
    STANDARD    = SVD_STANDARD
    GENERALIZED = SVD_GENERALIZED
    HYPERBOLIC  = SVD_HYPERBOLIC

class SVDErrorType(object):
    """
    SVD error type to assess accuracy of computed solutions.

    - `ABSOLUTE`: Absolute error.
    - `RELATIVE`: Relative error.
    - `NORM`:     Error relative to the matrix norm.

    See Also
    --------
    slepc.SVDErrorType
    """
    ABSOLUTE = SVD_ERROR_ABSOLUTE
    RELATIVE = SVD_ERROR_RELATIVE
    NORM     = SVD_ERROR_NORM

class SVDWhich(object):
    """
    SVD desired part of spectrum.

    - `LARGEST`:  Largest singular values.
    - `SMALLEST`: Smallest singular values.

    See Also
    --------
    slepc.SVDWhich
    """
    LARGEST  = SVD_LARGEST
    SMALLEST = SVD_SMALLEST

class SVDConv(object):
    """
    SVD convergence test.

    - `ABS`:   Absolute convergence test.
    - `REL`:   Convergence test relative to the singular value.
    - `NORM`:  Convergence test relative to the matrix norms.
    - `MAXIT`: No convergence until maximum number of iterations has been reached.
    - `USER`:  User-defined convergence test.

    See Also
    --------
    slepc.SVDConv
    """
    ABS   = SVD_CONV_ABS
    REL   = SVD_CONV_REL
    NORM  = SVD_CONV_NORM
    MAXIT = SVD_CONV_MAXIT
    USER  = SVD_CONV_USER

class SVDStop(object):
    """
    SVD stopping test.

    - `BASIC`:     Default stopping test.
    - `USER`:      User-defined stopping test.
    - `THRESHOLD`: Threshold stopping test.

    See Also
    --------
    slepc.SVDStop
    """
    BASIC     = SVD_STOP_BASIC
    USER      = SVD_STOP_USER
    THRESHOLD = SVD_STOP_THRESHOLD

class SVDConvergedReason(object):
    """
    SVD convergence reasons.

    - `CONVERGED_TOL`: All eigenpairs converged to requested tolerance.
    - `CONVERGED_USER`: User-defined convergence criterion satisfied.
    - `CONVERGED_MAXIT`: Maximum iterations completed in case MAXIT
      convergence criterion.
    - `DIVERGED_ITS`: Maximum number of iterations exceeded.
    - `DIVERGED_BREAKDOWN`: Solver failed due to breakdown.
    - `DIVERGED_SYMMETRY_LOST`: Underlying indefinite eigensolver was not able
      to keep symmetry.
    - `CONVERGED_ITERATING`: Iteration not finished yet.

    See Also
    --------
    slepc.SVDConvergedReason
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
    SVD TRLanczos bidiagonalization choices for the GSVD case.

    - `SINGLE`: Single bidiagonalization (:math:`Q_A`).
    - `UPPER`: Joint bidiagonalization, both :math:`Q_A` and :math:`Q_B`
      in upper bidiagonal form.
    - `LOWER`: Joint bidiagonalization, :math:`Q_A` lower bidiagonal,
      :math:`Q_B` upper bidiagonal.

    See Also
    --------
    slepc.SVDTRLanczosGBidiag
    """
    SINGLE = SVD_TRLANCZOS_GBIDIAG_SINGLE
    UPPER  = SVD_TRLANCZOS_GBIDIAG_UPPER
    LOWER  = SVD_TRLANCZOS_GBIDIAG_LOWER

# -----------------------------------------------------------------------------

cdef class SVD(Object):

    """
    Singular Value Decomposition Solver.

    The Singular Value Decomposition Solver (`SVD`) is very similar to the
    `EPS` object, but intended for the computation of the partial SVD of a
    rectangular matrix. With this type of object, the user can specify an
    SVD problem and solve it with any of the different solvers encapsulated
    by the package. Some of these solvers are actually implemented through
    calls to `EPS` eigensolvers.
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

    def view(self, Viewer viewer=None) -> None:
        """
        Print the SVD data structure.

        Collective.

        Parameters
        ----------
        viewer
            Visualization context; if not provided, the standard
            output is used.

        See Also
        --------
        slepc.SVDView
        """
        cdef PetscViewer vwr = def_Viewer(viewer)
        CHKERR( SVDView(self.svd, vwr) )

    def destroy(self) -> Self:
        """
        Destroy the SVD object.

        Collective.

        See Also
        --------
        slepc.SVDDestroy
        """
        CHKERR( SVDDestroy(&self.svd) )
        self.svd = NULL
        return self

    def reset(self) -> None:
        """
        Reset the SVD object.

        Collective.

        See Also
        --------
        slepc.SVDReset
        """
        CHKERR( SVDReset(self.svd) )

    def create(self, comm: Comm | None = None) -> Self:
        """
        Create the SVD object.

        Collective.

        Parameters
        ----------
        comm
            MPI communicator; if not provided, it defaults to all processes.

        See Also
        --------
        slepc.SVDCreate
        """
        cdef MPI_Comm ccomm = def_Comm(comm, SLEPC_COMM_DEFAULT())
        cdef SlepcSVD newsvd = NULL
        CHKERR( SVDCreate(ccomm, &newsvd) )
        CHKERR( SlepcCLEAR(self.obj) ); self.svd = newsvd
        return self

    def setType(self, svd_type: Type | str) -> None:
        """
        Set the particular solver to be used in the SVD object.

        Logically collective.

        Parameters
        ----------
        svd_type
            The solver to be used.

        Notes
        -----
        The default is `CROSS`. Normally, it is best to use
        `setFromOptions()` and then set the SVD type from the options
        database rather than by using this routine. Using the options
        database provides the user with maximum flexibility in
        evaluating the different available methods.

        See Also
        --------
        getType, slepc.SVDSetType
        """
        cdef SlepcSVDType cval = NULL
        svd_type = str2bytes(svd_type, &cval)
        CHKERR( SVDSetType(self.svd, cval) )

    def getType(self) -> str:
        """
        Get the SVD type of this object.

        Not collective.

        Returns
        -------
        str
            The solver currently being used.

        See Also
        --------
        setType, slepc.SVDGetType
        """
        cdef SlepcSVDType svd_type = NULL
        CHKERR( SVDGetType(self.svd, &svd_type) )
        return bytes2str(svd_type)

    def getOptionsPrefix(self) -> str:
        """
        Get the prefix used for searching for all SVD options in the database.

        Not collective.

        Returns
        -------
        str
            The prefix string set for this SVD object.

        See Also
        --------
        setOptionsPrefix, appendOptionsPrefix, slepc.SVDGetOptionsPrefix
        """
        cdef const char *prefix = NULL
        CHKERR( SVDGetOptionsPrefix(self.svd, &prefix) )
        return bytes2str(prefix)

    def setOptionsPrefix(self, prefix: str | None = None) -> None:
        """
        Set the prefix used for searching for all SVD options in the database.

        Logically collective.

        Parameters
        ----------
        prefix
            The prefix string to prepend to all SVD option requests.

        Notes
        -----
        A hyphen (-) must NOT be given at the beginning of the prefix
        name.  The first character of all runtime options is
        AUTOMATICALLY the hyphen.

        For example, to distinguish between the runtime options for
        two different SVD contexts, one could call::

            S1.setOptionsPrefix("svd1_")
            S2.setOptionsPrefix("svd2_")

        See Also
        --------
        appendOptionsPrefix, getOptionsPrefix, slepc.SVDGetOptionsPrefix
        """
        cdef const char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( SVDSetOptionsPrefix(self.svd, cval) )

    def appendOptionsPrefix(self, prefix: str | None = None) -> None:
        """
        Append to the prefix used for searching for all SVD options in the database.

        Logically collective.

        Parameters
        ----------
        prefix
            The prefix string to prepend to all SVD option requests.

        See Also
        --------
        setOptionsPrefix, getOptionsPrefix, slepc.SVDAppendOptionsPrefix
        """
        cdef const char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( SVDAppendOptionsPrefix(self.svd, cval) )

    def setFromOptions(self) -> None:
        """
        Set SVD options from the options database.

        Collective.

        Notes
        -----
        To see all options, run your program with the ``-help`` option.

        This routine must be called before `setUp()` if the user is to be
        allowed to set the solver type.

        See Also
        --------
        setOptionsPrefix, slepc.SVDSetFromOptions
        """
        CHKERR( SVDSetFromOptions(self.svd) )

    def getProblemType(self) -> ProblemType:
        """
        Get the problem type from the SVD object.

        Not collective.

        Returns
        -------
        ProblemType
            The problem type that was previously set.

        See Also
        --------
        setProblemType, slepc.SVDGetProblemType
        """
        cdef SlepcSVDProblemType val = SVD_STANDARD
        CHKERR( SVDGetProblemType(self.svd, &val) )
        return val

    def setProblemType(self, problem_type: ProblemType) -> None:
        """
        Set the type of the singular value problem.

        Logically collective.

        Parameters
        ----------
        problem_type
            The problem type to be set.

        Notes
        -----
        The GSVD requires that two matrices have been passed via
        `setOperators()`. The HSVD requires that a signature matrix
        has been passed via `setSignature()`.

        See Also
        --------
        setOperators, setSignature, getProblemType, slepc.SVDSetProblemType
        """
        cdef SlepcSVDProblemType val = problem_type
        CHKERR( SVDSetProblemType(self.svd, val) )

    def isGeneralized(self) -> bool:
        """
        Tell if the SVD corresponds to a generalized singular value problem.

        Not collective.

        Returns
        -------
        bool
            ``True`` if two matrices were set with `setOperators()`.

        See Also
        --------
        setProblemType, isHyperbolic, slepc.SVDIsGeneralized
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( SVDIsGeneralized(self.svd, &tval) )
        return toBool(tval)

    def isHyperbolic(self) -> bool:
        """
        Tell whether the SVD object corresponds to a hyperbolic singular value problem.

        Not collective.

        Returns
        -------
        bool
            ``True`` if the problem was specified as hyperbolic.

        See Also
        --------
        setProblemType, isGeneralized, slepc.SVDIsHyperbolic
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( SVDIsHyperbolic(self.svd, &tval) )
        return toBool(tval)

    #

    def getImplicitTranspose(self) -> bool:
        """
        Get the mode used to handle the transpose of the associated matrix.

        Not collective.

        Returns
        -------
        bool
            How to handle the transpose (implicitly or not).

        See Also
        --------
        setImplicitTranspose, slepc.SVDGetImplicitTranspose
        """
        cdef PetscBool val = PETSC_FALSE
        CHKERR( SVDGetImplicitTranspose(self.svd, &val) )
        return toBool(val)

    def setImplicitTranspose(self, mode: bool) -> None:
        """
        Set how to handle the transpose of the associated matrix.

        Logically collective.

        Parameters
        ----------
        impl
            How to handle the transpose (implicitly or not).

        Notes
        -----
        By default, the transpose of the matrix is explicitly built
        (if the matrix has defined the ``Mat.transpose()`` operation).

        If this flag is set to ``True``, the solver does not build the
        transpose, but handles it implicitly via ``Mat.multTranspose()``
        (or ``Mat.multHermitianTranspose()`` in the complex case).

        See Also
        --------
        getImplicitTranspose, slepc.SVDSetImplicitTranspose
        """
        cdef PetscBool val = asBool(mode)
        CHKERR( SVDSetImplicitTranspose(self.svd, val) )

    def getWhichSingularTriplets(self) -> Which:
        """
        Get which singular triplets are to be sought.

        Not collective.

        Returns
        -------
        Which
            The singular values to be sought (either largest or smallest).

        See Also
        --------
        setWhichSingularTriplets, slepc.SVDGetWhichSingularTriplets
        """
        cdef SlepcSVDWhich val = SVD_LARGEST
        CHKERR( SVDGetWhichSingularTriplets(self.svd, &val) )
        return val

    def setWhichSingularTriplets(self, which: Which) -> None:
        """
        Set which singular triplets are to be sought.

        Logically collective.

        Parameters
        ----------
        which
            The singular values to be sought (either largest or smallest).

        See Also
        --------
        getWhichSingularTriplets, slepc.SVDSetWhichSingularTriplets
        """
        cdef SlepcSVDWhich val = which
        CHKERR( SVDSetWhichSingularTriplets(self.svd, val) )

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

        See Also
        --------
        setThreshold, slepc.SVDGetThreshold
        """
        cdef PetscReal rval = 0
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( SVDGetThreshold(self.svd, &rval, &tval) )
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
        the threshold, where singular values are computed in sequence
        until one of the computed singular values is below/above the
        threshold (depending on whether largest or smallest singular
        values are computed).

        In the case of largest singular values, the threshold can be
        made relative with respect to the largest singular value
        (i.e., the matrix norm).

        The details are given in `slepc.SVDSetThreshold`.

        See Also
        --------
        setStoppingTest, getThreshold, slepc.SVDSetThreshold
        """
        cdef PetscReal rval = asReal(thres)
        cdef PetscBool tval = asBool(rel)
        CHKERR( SVDSetThreshold(self.svd, rval, tval) )

    def getTolerances(self) -> tuple[float, int]:
        """
        Get the tolerance and maximum iteration count.

        Not collective.

        Get the tolerance and maximum iteration count used by the default SVD
        convergence tests.

        Returns
        -------
        tol: float
            The convergence tolerance.
        max_it: int
            The maximum number of iterations.

        See Also
        --------
        setTolerances, slepc.SVDGetTolerances
        """
        cdef PetscReal rval = 0
        cdef PetscInt  ival = 0
        CHKERR( SVDGetTolerances(self.svd, &rval, &ival) )
        return (toReal(rval), toInt(ival))

    def setTolerances(self, tol: float | None = None, max_it: int | None = None) -> None:
        """
        Set the tolerance and maximum iteration count used.

        Logically collective.

        Set the tolerance and maximum iteration count used by the default SVD
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
        getTolerances, slepc.SVDSetTolerances
        """
        cdef PetscReal rval = PETSC_CURRENT
        cdef PetscInt  ival = PETSC_CURRENT
        if tol    is not None: rval = asReal(tol)
        if max_it is not None: ival = asInt(max_it)
        CHKERR( SVDSetTolerances(self.svd, rval, ival) )

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
        setConvergenceTest, slepc.SVDGetConvergenceTest
        """
        cdef SlepcSVDConv conv = SVD_CONV_REL
        CHKERR( SVDGetConvergenceTest(self.svd, &conv) )
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
        getConvergenceTest, slepc.SVDSetConvergenceTest
        """
        cdef SlepcSVDConv tconv = conv
        CHKERR( SVDSetConvergenceTest(self.svd, tconv) )

    def getTrackAll(self) -> bool:
        """
        Get the flag indicating if all residual norms must be computed or not.

        Not collective.

        Returns
        -------
        bool
            Whether the solver computes all residuals or not.

        See Also
        --------
        setTrackAll, slepc.SVDGetTrackAll
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( SVDGetTrackAll(self.svd, &tval) )
        return toBool(tval)

    def setTrackAll(self, trackall: bool) -> None:
        """
        Set flag to compute the residual of all singular triplets.

        Logically collective.

        Set if the solver must compute the residual of all approximate
        singular triplets or not.

        Parameters
        ----------
        trackall
            Whether to compute all residuals or not.

        See Also
        --------
        getTrackAll, slepc.SVDSetTrackAll
        """
        cdef PetscBool tval = asBool(trackall)
        CHKERR( SVDSetTrackAll(self.svd, tval) )

    def getDimensions(self) -> tuple[int, int, int]:
        """
        Get the number of singular values to compute and the dimension of the subspace.

        Not collective.

        Returns
        -------
        nsv: int
            Number of singular values to compute.
        ncv: int
            Maximum dimension of the subspace to be used by the solver.
        mpd: int
            Maximum dimension allowed for the projected problem.

        See Also
        --------
        setDimensions, slepc.SVDGetDimensions
        """
        cdef PetscInt ival1 = 0
        cdef PetscInt ival2 = 0
        cdef PetscInt ival3 = 0
        CHKERR( SVDGetDimensions(self.svd, &ival1, &ival2, &ival3) )
        return (toInt(ival1), toInt(ival2), toInt(ival3))

    def setDimensions(
        self,
        nsv: int | None = None,
        ncv: int | None = None,
        mpd: int | None = None,
    ) -> None:
        """
        Set the number of singular values to compute and the dimension of the subspace.

        Logically collective.

        Parameters
        ----------
        nsv
            Number of singular values to compute.
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

        + In cases where ``nsv`` is small, the user sets ``ncv``
          (a reasonable default is 2 * ``nsv``).

        + In cases where ``nsv`` is large, the user sets ``mpd``.

        The value of ``ncv`` should always be between ``nsv`` and (``nsv`` +
        ``mpd``), typically ``ncv`` = ``nsv`` + ``mpd``. If ``nsv`` is not too
        large, ``mpd`` = ``nsv`` is a reasonable choice, otherwise a
        smaller value should be used.

        See Also
        --------
        getDimensions, slepc.SVDSetDimensions
        """
        cdef PetscInt ival1 = PETSC_CURRENT
        cdef PetscInt ival2 = PETSC_CURRENT
        cdef PetscInt ival3 = PETSC_CURRENT
        if nsv is not None: ival1 = asInt(nsv)
        if ncv is not None: ival2 = asInt(ncv)
        if mpd is not None: ival3 = asInt(mpd)
        CHKERR( SVDSetDimensions(self.svd, ival1, ival2, ival3) )

    def getBV(self) -> tuple[BV, BV]:
        """
        Get the basis vectors objects associated to the SVD object.

        Not collective.

        Returns
        -------
        V: BV
            The basis vectors context for right singular vectors.
        U: BV
            The basis vectors context for left singular vectors.

        See Also
        --------
        setBV, slepc.SVDGetBV
        """
        cdef BV V = BV()
        cdef BV U = BV()
        CHKERR( SVDGetBV(self.svd, &V.bv, &U.bv) )
        CHKERR( PetscINCREF(V.obj) )
        CHKERR( PetscINCREF(U.obj) )
        return (V,U)

    def setBV(self, BV V,BV U=None) -> None:
        """
        Set basis vectors objects associated to the SVD solver.

        Collective.

        Parameters
        ----------
        V
            The basis vectors context for right singular vectors.
        U
            The basis vectors context for left singular vectors.

        See Also
        --------
        getBV, slepc.SVDSetBV
        """
        cdef SlepcBV VBV = V.bv
        cdef SlepcBV UBV = U.bv if U is not None else <SlepcBV>NULL
        CHKERR( SVDSetBV(self.svd, VBV, UBV) )

    def getDS(self) -> DS:
        """
        Get the direct solver associated to the singular value solver.

        Not collective.

        Returns
        -------
        DS
            The direct solver context.

        See Also
        --------
        setDS, slepc.SVDGetDS
        """
        cdef DS ds = DS()
        CHKERR( SVDGetDS(self.svd, &ds.ds) )
        CHKERR( PetscINCREF(ds.obj) )
        return ds

    def setDS(self, DS ds) -> None:
        """
        Set a direct solver object associated to the singular value solver.

        Collective.

        Parameters
        ----------
        ds
            The direct solver context.

        See Also
        --------
        getDS, slepc.SVDSetDS
        """
        CHKERR( SVDSetDS(self.svd, ds.ds) )

    def getOperators(self) -> tuple[Mat, Mat] | tuple[Mat, None]:
        """
        Get the matrices associated with the singular value problem.

        Collective.

        Returns
        -------
        A: petsc4py.PETSc.Mat
            The matrix associated with the singular value problem.
        B: petsc4py.PETSc.Mat
            The second matrix in the case of GSVD.

        See Also
        --------
        setOperators, slepc.SVDGetOperators
        """
        cdef Mat A = Mat()
        cdef Mat B = Mat()
        CHKERR( SVDGetOperators(self.svd, &A.mat, &B.mat) )
        CHKERR( PetscINCREF(A.obj) )
        if B.mat:
            CHKERR( PetscINCREF(B.obj) )
            return (A, B)
        else:
            return (A, None)

    def setOperators(self, Mat A, Mat B=None) -> None:
        """
        Set the matrices associated with the singular value problem.

        Collective.

        Parameters
        ----------
        A
            The matrix associated with the singular value problem.
        B
            The second matrix in the case of GSVD.

        See Also
        --------
        getOperators, slepc.SVDSetOperators
        """
        cdef PetscMat Bmat = B.mat if B is not None else <PetscMat>NULL
        CHKERR( SVDSetOperators(self.svd, A.mat, Bmat) )

    def getSignature(self, Vec omega = None) -> Vec:
        """
        Get the signature matrix defining a hyperbolic singular value problem.

        Collective.

        Parameters
        ----------
        omega
            Optional vector to store the diagonal elements of the signature matrix.

        Returns
        -------
        petsc4py.PETSc.Vec
            A vector containing the diagonal elements of the signature matrix.

        See Also
        --------
        setSignature, slepc.SVDGetSignature
        """
        cdef PetscMat A = NULL
        if omega is None:
            omega = Vec()
        if omega.vec == NULL:
            CHKERR( SVDGetOperators(self.svd, &A, <PetscMat*>NULL) )
            CHKERR( MatCreateVecs(A, <PetscVec*>NULL, &omega.vec) )
        CHKERR( SVDGetSignature(self.svd, omega.vec) )
        return omega

    def setSignature(self, Vec omega=None) -> None:
        """
        Set the signature matrix defining a hyperbolic singular value problem.

        Collective.

        Parameters
        ----------
        omega
            A vector containing the diagonal elements of the signature matrix.

        See Also
        --------
        getSignature, slepc.SVDSetSignature
        """
        cdef PetscVec Ovec = omega.vec if omega is not None else <PetscVec>NULL
        CHKERR( SVDSetSignature(self.svd, Ovec) )

    #

    def setInitialSpace(
        self,
        spaceright: list[Vec] | None = None,
        spaceleft: list[Vec] | None = None,
    ) -> None:
        """
        Set the initial spaces from which the SVD solver starts to iterate.

        Collective.

        Parameters
        ----------
        spaceright
            The right initial space.
        spaceleft
            The left initial space.

        Notes
        -----
        The initial right and left spaces are rough approximations to the
        right and/or left singular subspaces from which the solver starts
        to iterate. It is not necessary to provide both sets of vectors.

        Some solvers start to iterate on a single vector (initial vector).
        In that case, the other vectors are ignored.

        These vectors do not persist from one `solve()` call to the other,
        so the initial spaces should be set every time.

        The vectors do not need to be mutually orthonormal, since they are
        explicitly orthonormalized internally.

        Common usage of this function is when the user can provide a rough
        approximation of the wanted singular spaces. Then, convergence may
        be faster.

        See Also
        --------
        slepc.SVDSetInitialSpaces
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

    def setStoppingTest(
        self,
        stopping: SVDStoppingFunction | None,
        args: tuple[Any, ...] | None = None,
        kargs: dict[str, Any] | None = None,
    ) -> None:
        """
        Set a function to decide when to stop the outer iteration of the eigensolver.

        Logically collective.

        See Also
        --------
        getStoppingTest, slepc.SVDSetStoppingTestFunction
        """
        if stopping is not None:
            if args is None: args = ()
            if kargs is None: kargs = {}
            self.set_attr('__stopping__', (stopping, args, kargs))
            CHKERR( SVDSetStoppingTestFunction(self.svd, SVD_Stopping, NULL, NULL) )
        else:
            self.set_attr('__stopping__', None)
            CHKERR( SVDSetStoppingTestFunction(self.svd, SVDStoppingBasic, NULL, NULL) )

    def getStoppingTest(self) -> SVDStoppingFunction:
        """
        Get the stopping test function.

        Not collective.

        Returns
        -------
        SVDStoppingFunction
            The stopping test function.

        See Also
        --------
        setStoppingTest
        """
        return self.get_attr('__stopping__')

    #

    def setMonitor(
        self,
        monitor: SVDMonitorFunction | None,
        args: tuple[Any, ...] | None = None,
        kargs: dict[str, Any] | None = None,
    ) -> None:
        """
        Append a monitor function to the list of monitors.

        Logically collective.

        See Also
        --------
        getMonitor, cancelMonitor, slepc.SVDMonitorSet
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

    def getMonitor(self) -> SVDMonitorFunction:
        """
        Get the list of monitor functions.

        Not collective.

        Returns
        -------
        SVDMonitorFunction
            The list of monitor functions.

        See Also
        --------
        setMonitor
        """
        return self.get_attr('__monitor__')

    def cancelMonitor(self) -> None:
        """
        Clear all monitors for an `SVD` object.

        Logically collective.

        See Also
        --------
        slepc.SVDMonitorCancel
        """
        CHKERR( SVDMonitorCancel(self.svd) )
        self.set_attr('__monitor__', None)

    #

    def setUp(self) -> None:
        """
        Set up all the internal data structures.

        Collective.

        Notes
        -----
        Sets up all the internal data structures necessary for the execution
        of the singular value solver.

        This function need not be called explicitly in most cases,
        since `solve()` calls it. It can be useful when one wants to
        measure the set-up time separately from the solve time.

        See Also
        --------
        solve, slepc.SVDSetUp
        """
        CHKERR( SVDSetUp(self.svd) )

    def solve(self) -> None:
        """
        Solve the singular value problem.

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
        setUp, setOperators, getConverged, getConvergedReason, slepc.SVDSolve
        """
        CHKERR( SVDSolve(self.svd) )

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
        getConvergedReason, setTolerances, slepc.SVDGetIterationNumber
        """
        cdef PetscInt ival = 0
        CHKERR( SVDGetIterationNumber(self.svd, &ival) )
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
        setTolerances, solve, slepc.SVDGetConvergedReason
        """
        cdef SlepcSVDConvergedReason val = SVD_CONVERGED_ITERATING
        CHKERR( SVDGetConvergedReason(self.svd, &val) )
        return val

    def getConverged(self) -> int:
        """
        Get the number of converged singular triplets.

        Not collective.

        Returns
        -------
        nconv: int
            Number of converged singular triplets.

        Notes
        -----
        This function should be called after `solve()` has finished.

        The value ``nconv`` may be different from the number of requested
        solutions ``nsv``, but not larger than ``ncv``, see `setDimensions()`.

        See Also
        --------
        setDimensions, solve, getValue, slepc.SVDGetConverged
        """
        cdef PetscInt ival = 0
        CHKERR( SVDGetConverged(self.svd, &ival) )
        return toInt(ival)

    def getValue(self, i: int) -> float:
        """
        Get the i-th singular value as computed by `solve()`.

        Collective.

        Parameters
        ----------
        i
            Index of the solution to be obtained.

        Returns
        -------
        float
            The computed singular value.

        Notes
        -----
        The index ``i`` should be a value between ``0`` and
        ``nconv-1`` (see `getConverged()`. Singular triplets are
        indexed according to the ordering criterion established with
        `setWhichSingularTriplets()`.

        See Also
        --------
        getConverged, setWhichSingularTriplets, slepc.SVDGetSingularTriplet
        """
        cdef PetscReal rval = 0
        CHKERR( SVDGetSingularTriplet(self.svd, i, &rval, NULL, NULL) )
        return toReal(rval)

    def getVectors(self, i: int, Vec U, Vec V) -> None:
        """
        Get the i-th left and right singular vectors as computed by `solve()`.

        Collective.

        Parameters
        ----------
        i
            Index of the solution to be obtained.
        U
            Placeholder for the returned left singular vector.
        V
            Placeholder for the returned right singular vector.

        Notes
        -----
        The index ``i`` should be a value between ``0`` and
        ``nconv-1`` (see `getConverged()`. Singular triplets are
        indexed according to the ordering criterion established with
        `setWhichSingularTriplets()`.

        See Also
        --------
        getConverged, setWhichSingularTriplets, slepc.SVDGetSingularTriplet
        """
        cdef PetscReal dummy = 0
        CHKERR( SVDGetSingularTriplet(self.svd, i, &dummy, U.vec, V.vec) )

    def getSingularTriplet(self, i: int, Vec U=None, Vec V=None) -> float:
        """
        Get the i-th triplet of the singular value decomposition.

        Collective.

        Get the i-th triplet of the singular value decomposition as computed
        by `solve()`. The solution consists of the singular value and its left
        and right singular vectors.

        Parameters
        ----------
        i
            Index of the solution to be obtained.
        U
            Placeholder for the returned left singular vector.
        V
            Placeholder for the returned right singular vector.

        Returns
        -------
        float
            The computed singular value.

        Notes
        -----
        The index ``i`` should be a value between ``0`` and
        ``nconv-1`` (see `getConverged()`. Singular triplets are
        indexed according to the ordering criterion established with
        `setWhichSingularTriplets()`.

        See Also
        --------
        getConverged, setWhichSingularTriplets, slepc.SVDGetSingularTriplet
        """
        cdef PetscReal rval = 0
        cdef PetscVec Uvec = U.vec if U is not None else <PetscVec>NULL
        cdef PetscVec Vvec = V.vec if V is not None else <PetscVec>NULL
        CHKERR( SVDGetSingularTriplet(self.svd, i, &rval, Uvec, Vvec) )
        return toReal(rval)

    #

    def computeError(self, i: int, etype: ErrorType | None = None) -> float:
        """
        Compute the error associated with the i-th singular triplet.

        Collective.

        Compute the error (based on the residual norm) associated with the
        i-th singular triplet.

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
            :math:`\sqrt{\eta_1^2+\eta_2^2}` where
            :math:`\eta_1 = \|A v - \sigma u\|_2`,
            :math:`\eta_2 = \|A^* u - \sigma v\|_2`, :math:`\sigma` is the
            approximate singular value, :math:`u` and :math:`v` are the left
            and right singular vectors.

        Notes
        -----
        The index ``i`` should be a value between ``0`` and ``nconv-1``
        (see `getConverged()`).

        In the case of the GSVD, the two components of the residual norm are
        :math:`\eta_1 = \|s^2 A^*u-cB^*Bx\|_2` and
        :math:`\eta_2 = ||c^2 B^*v-sA^*Ax||_2`, where :math:`(\sigma,u,v,x)`
        is the approximate generalized singular quadruple, with
        :math:`\sigma=c/s`.

        See Also
        --------
        solve, slepc.SVDComputeError
        """
        cdef SlepcSVDErrorType et = SVD_ERROR_RELATIVE
        cdef PetscReal rval = 0
        if etype is not None: et = etype
        CHKERR( SVDComputeError(self.svd, i, et, &rval) )
        return toReal(rval)

    def errorView(self, etype: ErrorType | None = None, viewer: petsc4py.PETSc.Viewer | None = None) -> None:
        """
        Display the errors associated with the computed solution.

        Collective.

        Display the errors and the singular values.

        Parameters
        ----------
        etype
            The error type to compute.
        viewer
            Visualization context; if not provided, the standard
            output is used.

        Notes
        -----
        By default, this function checks the error of all singular triplets and
        prints the singular values if all of them are below the requested
        tolerance. If the viewer has format ``ASCII_INFO_DETAIL`` then a table
        with singular values and corresponding errors is printed.

        See Also
        --------
        solve, valuesView, vectorsView, slepc.SVDErrorView
        """
        cdef SlepcSVDErrorType et = SVD_ERROR_RELATIVE
        if etype is not None: et = etype
        cdef PetscViewer vwr = def_Viewer(viewer)
        CHKERR( SVDErrorView(self.svd, et, vwr) )

    def valuesView(self, viewer: Viewer | None = None) -> None:
        """
        Display the computed singular values in a viewer.

        Collective.

        Parameters
        ----------
        viewer
            Visualization context; if not provided, the standard
            output is used.

        See Also
        --------
        solve, vectorsView, errorView, slepc.SVDValuesView
        """
        cdef PetscViewer vwr = def_Viewer(viewer)
        CHKERR( SVDValuesView(self.svd, vwr) )

    def vectorsView(self, viewer: Viewer | None = None) -> None:
        """
        Output computed singular vectors to a viewer.

        Collective.

        Parameters
        ----------
        viewer
            Visualization context; if not provided, the standard
            output is used.

        See Also
        --------
        solve, valuesView, errorView, slepc.SVDVectorsView
        """
        cdef PetscViewer vwr = def_Viewer(viewer)
        CHKERR( SVDVectorsView(self.svd, vwr) )

    #

    def setCrossEPS(self, EPS eps) -> None:
        """
        Set an eigensolver object associated to the singular value solver.

        Collective.

        Parameters
        ----------
        eps
            The eigensolver object.

        See Also
        --------
        getCrossEPS, slepc.SVDCrossSetEPS
        """
        CHKERR( SVDCrossSetEPS(self.svd, eps.eps) )

    def getCrossEPS(self) -> EPS:
        """
        Get the eigensolver object associated to the singular value solver.

        Collective.

        Returns
        -------
        EPS
            The eigensolver object.

        See Also
        --------
        setCrossEPS, slepc.SVDCrossGetEPS
        """
        cdef EPS eps = EPS()
        CHKERR( SVDCrossGetEPS(self.svd, &eps.eps) )
        CHKERR( PetscINCREF(eps.obj) )
        return eps

    def setCrossExplicitMatrix(self, flag: bool = True) -> None:
        """
        Set if the eigensolver operator :math:`A^*A` must be computed.

        Logically collective.

        Parameters
        ----------
        flag
            ``True`` to build :math:`A^*A` explicitly.

        Notes
        -----
        In GSVD there are two cross product matrices, :math:`A^*A` and
        :math:`B^*B`. In HSVD the expression for the cross product matrix
        is different, :math:`A^*\Omega A`.

        By default the matrices are not built explicitly, but handled as
        shell matrices

        See Also
        --------
        getCrossExplicitMatrix, slepc.SVDCrossSetExplicitMatrix
        """
        cdef PetscBool tval = asBool(flag)
        CHKERR( SVDCrossSetExplicitMatrix(self.svd, tval) )

    def getCrossExplicitMatrix(self) -> bool:
        """
        Get the flag indicating if :math:`A^*A` is built explicitly.

        Not collective.

        Returns
        -------
        bool
            ``True`` if :math:`A^*A` is built explicitly.

        See Also
        --------
        setCrossExplicitMatrix, slepc.SVDCrossGetExplicitMatrix
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( SVDCrossGetExplicitMatrix(self.svd, &tval) )
        return toBool(tval)

    def setCyclicEPS(self, EPS eps) -> None:
        """
        Set an eigensolver object associated to the singular value solver.

        Collective.

        Parameters
        ----------
        eps
            The eigensolver object.

        See Also
        --------
        getCyclicEPS, slepc.SVDCyclicSetEPS
        """
        CHKERR( SVDCyclicSetEPS(self.svd, eps.eps) )

    def getCyclicEPS(self) -> EPS:
        """
        Get the eigensolver object associated to the singular value solver.

        Collective.

        Returns
        -------
        EPS
            The eigensolver object.

        See Also
        --------
        setCyclicEPS, slepc.SVDCyclicGetEPS
        """
        cdef EPS eps = EPS()
        CHKERR( SVDCyclicGetEPS(self.svd, &eps.eps) )
        CHKERR( PetscINCREF(eps.obj) )
        return eps

    def setCyclicExplicitMatrix(self, flag: bool = True) -> None:
        """
        Set if the eigensolver operator :math:`H(A)` must be computed explicitly.

        Logically collective.

        Set if the eigensolver operator :math:`H(A) = [ 0\; A ; A^T\; 0 ]`
        must be computed explicitly.

        Parameters
        ----------
        flag
            ``True`` if :math:`H(A)` must be built explicitly.

        Notes
        -----
        In GSVD and HSVD the equivalent eigenvalue problem has
        generalized form, and hence two matrices are built.

        By default the matrices are not built explicitly, but handled as
        shell matrices.

        See Also
        --------
        getCyclicExplicitMatrix, slepc.SVDCyclicSetExplicitMatrix
        """
        cdef PetscBool tval = asBool(flag)
        CHKERR( SVDCyclicSetExplicitMatrix(self.svd, tval) )

    def getCyclicExplicitMatrix(self) -> bool:
        """
        Get the flag indicating if :math:`H(A)` is built explicitly.

        Not collective.

        Get the flag indicating if :math:`H(A) = [ 0\; A ; A^T\; 0 ]`
        is built explicitly.

        Returns
        -------
        bool
            ``True`` if :math:`H(A)` is built explicitly.

        See Also
        --------
        setCyclicExplicitMatrix, slepc.SVDCyclicGetExplicitMatrix
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( SVDCyclicGetExplicitMatrix(self.svd, &tval) )
        return toBool(tval)

    def setLanczosOneSide(self, flag: bool = True) -> None:
        """
        Set if the variant of the Lanczos method to be used is one-sided or two-sided.

        Logically collective.

        Parameters
        ----------
        flag
            ``True`` if the method is one-sided.

        Notes
        -----
        By default, a two-sided variant is selected, which is
        sometimes slightly more robust. However, the one-sided variant
        is faster because it avoids the orthogonalization associated
        to left singular vectors. It also saves the memory required
        for storing such vectors.

        See Also
        --------
        getLanczosOneSide, slepc.SVDLanczosSetOneSide
        """
        cdef PetscBool tval = asBool(flag)
        CHKERR( SVDLanczosSetOneSide(self.svd, tval) )

    def getLanczosOneSide(self) -> bool:
        """
        Get if the variant of the Lanczos method to be used is one-sided or two-sided.

        Not collective.

        Returns
        -------
        bool
            ``True`` if the method is one-sided.

        See Also
        --------
        setLanczosOneSide, slepc.SVDLanczosGetOneSide
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( SVDLanczosGetOneSide(self.svd, &tval) )
        return toBool(tval)

    def setTRLanczosOneSide(self, flag: bool = True) -> None:
        """
        Set if the variant of the method to be used is one-sided or two-sided.

        Logically collective.

        Set if the variant of the thick-restart Lanczos method to be used is
        one-sided or two-sided.

        Parameters
        ----------
        flag
            ``True`` if the method is one-sided.

        Notes
        -----
        By default, a two-sided variant is selected, which is
        sometimes slightly more robust. However, the one-sided variant
        is faster because it avoids the orthogonalization associated
        to left singular vectors.

        See Also
        --------
        getTRLanczosOneSide, slepc.SVDLanczosSetOneSide
        """
        cdef PetscBool tval = asBool(flag)
        CHKERR( SVDLanczosSetOneSide(self.svd, tval) )

    def getTRLanczosOneSide(self) -> bool:
        """
        Get if the variant of the method to be used is one-sided or two-sided.

        Not collective.

        Get if the variant of the thick-restart Lanczos method to be used is
        one-sided or two-sided.

        Returns
        -------
        bool
            ``True`` if the method is one-sided.

        See Also
        --------
        setTRLanczosOneSide, slepc.SVDLanczosGetOneSide
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( SVDTRLanczosGetOneSide(self.svd, &tval) )
        return toBool(tval)

    def setTRLanczosGBidiag(self, bidiag: TRLanczosGBidiag) -> None:
        """
        Set the bidiagonalization choice to use in the GSVD TRLanczos solver.

        Logically collective.

        Parameters
        ----------
        bidiag
            The bidiagonalization choice.

        See Also
        --------
        getTRLanczosGBidiag, slepc.SVDTRLanczosSetGBidiag
        """
        cdef SlepcSVDTRLanczosGBidiag val = bidiag
        CHKERR( SVDTRLanczosSetGBidiag(self.svd, val) )

    def getTRLanczosGBidiag(self) -> TRLanczosGBidiag:
        """
        Get bidiagonalization choice used in the GSVD TRLanczos solver.

        Not collective.

        Returns
        -------
        TRLanczosGBidiag
            The bidiagonalization choice.

        See Also
        --------
        setTRLanczosGBidiag, slepc.SVDTRLanczosGetGBidiag
        """
        cdef SlepcSVDTRLanczosGBidiag val = SVD_TRLANCZOS_GBIDIAG_LOWER
        CHKERR( SVDTRLanczosGetGBidiag(self.svd, &val) )
        return val

    def setTRLanczosRestart(self, keep: float) -> None:
        """
        Set the restart parameter for the thick-restart Lanczos method.

        Logically collective.

        Set the restart parameter for the thick-restart Lanczos method, in
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
        getTRLanczosRestart, slepc.SVDTRLanczosSetRestart
        """
        cdef PetscReal val = asReal(keep)
        CHKERR( SVDTRLanczosSetRestart(self.svd, val) )

    def getTRLanczosRestart(self) -> float:
        """
        Get the restart parameter used in the thick-restart Lanczos method.

        Not collective.

        Returns
        -------
        float
            The number of vectors to be kept at restart.

        See Also
        --------
        setTRLanczosRestart, slepc.SVDTRLanczosGetRestart
        """
        cdef PetscReal val = 0
        CHKERR( SVDTRLanczosGetRestart(self.svd, &val) )
        return toReal(val)

    def setTRLanczosLocking(self, lock: bool) -> None:
        """
        Toggle between locking and non-locking variants of TRLanczos.

        Logically collective.

        Parameters
        ----------
        lock
            ``True`` if the locking variant must be selected.

        Notes
        -----
        The default is to lock converged singular triplets when the method restarts.
        This behavior can be changed so that all directions are kept in the
        working subspace even if already converged to working accuracy (the
        non-locking variant).

        See Also
        --------
        getTRLanczosLocking, slepc.SVDTRLanczosSetLocking
        """
        cdef PetscBool val = asBool(lock)
        CHKERR( SVDTRLanczosSetLocking(self.svd, val) )

    def getTRLanczosLocking(self) -> bool:
        """
        Get the locking flag used in the thick-restart Lanczos method.

        Not collective.

        Returns
        -------
        bool
            The locking flag.

        See Also
        --------
        setTRLanczosLocking, slepc.SVDTRLanczosGetLocking
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( SVDTRLanczosGetLocking(self.svd, &tval) )
        return toBool(tval)

    def setTRLanczosKSP(self, KSP ksp) -> None:
        """
        Set a linear solver object associated to the SVD solver.

        Collective.

        Parameters
        ----------
        ``ksp``
            The linear solver object.

        See Also
        --------
        getTRLanczosKSP, slepc.SVDTRLanczosSetKSP
        """
        CHKERR( SVDTRLanczosSetKSP(self.svd, ksp.ksp) )

    def getTRLanczosKSP(self) -> KSP:
        """
        Get the linear solver object associated with the SVD solver.

        Collective.

        Returns
        -------
        `petsc4py.PETSc.KSP`
            The linear solver object.

        See Also
        --------
        setTRLanczosKSP, slepc.SVDTRLanczosGetKSP
        """
        cdef KSP ksp = KSP()
        CHKERR( SVDTRLanczosGetKSP(self.svd, &ksp.ksp) )
        CHKERR( PetscINCREF(ksp.obj) )
        return ksp

    def setTRLanczosExplicitMatrix(self, flag: bool = True) -> None:
        """
        Set if the matrix :math:`Z=[A^*,B^*]^*` must be built explicitly.

        Logically collective.

        Parameters
        ----------
        flag
            ``True`` if :math:`Z=[A^*,B^*]^*` is built explicitly.

        Notes
        -----
        This option is relevant for the GSVD case only. :math:`Z` is the
        coefficient matrix of the least-squares solver used internally.

        See Also
        --------
        getTRLanczosExplicitMatrix, slepc.SVDTRLanczosSetExplicitMatrix
        """
        cdef PetscBool tval = asBool(flag)
        CHKERR( SVDTRLanczosSetExplicitMatrix(self.svd, tval) )

    def getTRLanczosExplicitMatrix(self) -> bool:
        """
        Get the flag indicating if :math:`Z=[A^*,B^*]^*` is built explicitly.

        Not collective.

        Returns
        -------
        bool
            ``True`` if :math:`Z=[A^*,B^*]^*` is built explicitly.

        See Also
        --------
        setTRLanczosExplicitMatrix, slepc.SVDTRLanczosGetExplicitMatrix
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( SVDTRLanczosGetExplicitMatrix(self.svd, &tval) )
        return toBool(tval)

    setOperator = setOperators  # backward compatibility

    #

    property problem_type:
        """The type of the eigenvalue problem."""
        def __get__(self) -> SVDProblemType:
            return self.getProblemType()
        def __set__(self, value):
            self.setProblemType(value)

    property transpose_mode:
        """How to handle the transpose of the matrix."""
        def __get__(self) -> bool:
            return self.getTransposeMode()
        def __set__(self, value):
            self.setTransposeMode(value)

    property which:
        """The portion of the spectrum to be sought."""
        def __get__(self) -> SVDWhich:
            return self.getWhichSingularTriplets()
        def __set__(self, value):
            self.setWhichSingularTriplets(value)

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

    property ds:
        """The direct solver (`DS`) object associated."""
        def __get__(self) -> DS:
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
