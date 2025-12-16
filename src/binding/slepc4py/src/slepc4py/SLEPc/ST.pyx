# -----------------------------------------------------------------------------

class STType(object):
    """
    ST type.

    - `SHIFT`:   Shift from origin.
    - `SINVERT`: Shift-and-invert.
    - `CAYLEY`:  Cayley transform.
    - `PRECOND`: Preconditioner.
    - `FILTER`:  Polynomial filter.
    - `SHELL`:   User-defined.

    See Also
    --------
    slepc.STType
    """
    SHIFT   = S_(STSHIFT)
    SINVERT = S_(STSINVERT)
    CAYLEY  = S_(STCAYLEY)
    PRECOND = S_(STPRECOND)
    FILTER  = S_(STFILTER)
    SHELL   = S_(STSHELL)

class STMatMode(object):
    """
    ST matrix mode.

    - `COPY`:    A working copy of the matrix is created.
    - `INPLACE`: The operation is computed in-place.
    - `SHELL`:   The matrix :math:`A - \sigma B` is handled as an
      implicit matrix.

    See Also
    --------
    slepc.STMatMode
    """
    COPY    = ST_MATMODE_COPY
    INPLACE = ST_MATMODE_INPLACE
    SHELL   = ST_MATMODE_SHELL

class STFilterType(object):
    """
    ST filter type.

    - ``FILTLAN``:  An adapted implementation of the Filtered Lanczos Package.
    - ``CHEBYSEV``: A polynomial filter based on a truncated Chebyshev series.

    See Also
    --------
    slepc.STFilterType
    """
    FILTLAN   = ST_FILTER_FILTLAN
    CHEBYSHEV = ST_FILTER_CHEBYSHEV

class STFilterDamping(object):
    """
    ST filter damping.

    - `NONE`:    No damping
    - `JACKSON`: Jackson damping
    - `LANCZOS`: Lanczos damping
    - `FEJER`:   Fejer damping

    See Also
    --------
    slepc.STFilterDamping
    """
    NONE    = ST_FILTER_DAMPING_NONE
    JACKSON = ST_FILTER_DAMPING_JACKSON
    LANCZOS = ST_FILTER_DAMPING_LANCZOS
    FEJER   = ST_FILTER_DAMPING_FEJER

# -----------------------------------------------------------------------------

cdef class ST(Object):

    """
    Spectral Transformation.

    The Spectral Transformation (`ST`) class encapsulates the functionality
    required for acceleration techniques based on the transformation of the
    spectrum. The eigensolvers implemented in `EPS` work by applying an
    operator to a set of vectors and this operator can adopt different forms.
    The `ST` object handles all the different possibilities in a uniform way,
    so that the solver can proceed without knowing which transformation has
    been selected. Polynomial eigensolvers in `PEP` also support spectral
    transformation.
    """

    Type          = STType
    MatMode       = STMatMode
    FilterType    = STFilterType
    FilterDamping = STFilterDamping

    def __cinit__(self):
        self.obj = <PetscObject*> &self.st
        self.st = NULL

    def view(self, Viewer viewer=None) -> None:
        """
        Print the ST data structure.

        Collective.

        Parameters
        ----------
        viewer
            Visualization context; if not provided, the standard
            output is used.

        See Also
        --------
        slepc.STView
        """
        cdef PetscViewer vwr = def_Viewer(viewer)
        CHKERR( STView(self.st, vwr) )

    def destroy(self) -> Self:
        """
        Destroy the ST object.

        Collective.

        See Also
        --------
        slepc.STDestroy
        """
        CHKERR( STDestroy(&self.st) )
        self.st = NULL
        return self

    def reset(self) -> None:
        """
        Reset the ST object.

        Collective.

        See Also
        --------
        slepc.STReset
        """
        CHKERR( STReset(self.st) )

    def create(self, comm: Comm | None = None) -> Self:
        """
        Create the ST object.

        Collective.

        Parameters
        ----------
        comm
            MPI communicator; if not provided, it defaults to all processes.

        See Also
        --------
        slepc.STCreate
        """
        cdef MPI_Comm ccomm = def_Comm(comm, SLEPC_COMM_DEFAULT())
        cdef SlepcST newst = NULL
        CHKERR( STCreate(ccomm, &newst) )
        CHKERR( SlepcCLEAR(self.obj) ); self.st = newst
        return self

    def setType(self, st_type: Type | str) -> None:
        """
        Set the particular spectral transformation to be used.

        Logically collective.

        Parameters
        ----------
        st_type
            The spectral transformation to be used.

        Notes
        -----
        The default is `SHIFT` with a zero shift. Normally, it is best
        to use `setFromOptions()` and then set the ST type from the
        options database rather than by using this routine. Using the
        options database provides the user with maximum flexibility in
        evaluating the different available methods.

        See Also
        --------
        getType, slepc.STSetType
        """
        cdef SlepcSTType cval = NULL
        st_type = str2bytes(st_type, &cval)
        CHKERR( STSetType(self.st, cval) )

    def getType(self) -> str:
        """
        Get the ST type of this object.

        Not collective.

        Returns
        -------
        str
            The spectral transformation currently being used.

        See Also
        --------
        setType, slepc.STGetType
        """
        cdef SlepcSTType st_type = NULL
        CHKERR( STGetType(self.st, &st_type) )
        return bytes2str(st_type)

    def setOptionsPrefix(self, prefix: str | None = None) -> None:
        """
        Set the prefix used for searching for all ST options in the database.

        Logically collective.

        Parameters
        ----------
        prefix
            The prefix string to prepend to all ST option requests.

        Notes
        -----
        A hyphen (``-``) must NOT be given at the beginning of the
        prefix name.  The first character of all runtime options is
        AUTOMATICALLY the hyphen.

        See Also
        --------
        appendOptionsPrefix, getOptionsPrefix, slepc.STGetOptionsPrefix
        """
        cdef const char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( STSetOptionsPrefix(self.st, cval) )

    def getOptionsPrefix(self) -> str:
        """
        Get the prefix used for searching for all ST options in the database.

        Not collective.

        Returns
        -------
        str
            The prefix string set for this ST object.

        See Also
        --------
        setOptionsPrefix, appendOptionsPrefix, slepc.STGetOptionsPrefix
        """
        cdef const char *prefix = NULL
        CHKERR( STGetOptionsPrefix(self.st, &prefix) )
        return bytes2str(prefix)

    def appendOptionsPrefix(self, prefix: str | None = None) -> None:
        """
        Append to the prefix used for searching for all ST options in the database.

        Logically collective.

        Parameters
        ----------
        prefix
            The prefix string to prepend to all ST option requests.

        See Also
        --------
        setOptionsPrefix, getOptionsPrefix, slepc.STAppendOptionsPrefix
        """
        cdef const char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( STAppendOptionsPrefix(self.st, cval) )

    def setFromOptions(self) -> None:
        """
        Set ST options from the options database.

        Collective.

        Notes
        -----
        To see all options, run your program with the ``-help`` option.

        This routine must be called before `setUp()` if the user is to be
        allowed to set the solver type.

        See Also
        --------
        setOptionsPrefix, slepc.STSetFromOptions
        """
        CHKERR( STSetFromOptions(self.st) )

    #

    def setShift(self, shift: Scalar) -> None:
        """
        Set the shift associated with the spectral transformation.

        Collective.

        Parameters
        ----------
        shift
            The value of the shift.

        Notes
        -----
        In some spectral transformations, changing the shift may have
        associated a lot of work, for example recomputing a
        factorization.

        This function is normally not directly called by users, since the
        shift is indirectly set by `EPS.setTarget()`.

        See Also
        --------
        getShift, slepc.STSetShift
        """
        cdef PetscScalar sval = asScalar(shift)
        CHKERR( STSetShift(self.st, sval) )

    def getShift(self) -> Scalar:
        """
        Get the shift associated with the spectral transformation.

        Not collective.

        Returns
        -------
        Scalar
            The value of the shift.

        See Also
        --------
        setShift, slepc.STGetShift
        """
        cdef PetscScalar sval = 0
        CHKERR( STGetShift(self.st, &sval) )
        return toScalar(sval)

    def setTransform(self, flag: bool = True) -> None:
        """
        Set a flag to indicate whether the transformed matrices are computed or not.

        Logically collective.

        Parameters
        ----------
        flag
            This flag is intended for the case of polynomial
            eigenproblems solved via linearization.
            If this flag is ``False`` (default) the spectral transformation
            is applied to the linearization (handled by the eigensolver),
            otherwise it is applied to the original problem.

        See Also
        --------
        getTransform, slepc.STSetTransform
        """
        cdef PetscBool sval = asBool(flag)
        CHKERR( STSetTransform(self.st, sval) )

    def getTransform(self) -> bool:
        """
        Get the flag indicating whether the transformed matrices are computed or not.

        Not collective.

        Returns
        -------
        bool
            This flag is intended for the case of polynomial
            eigenproblems solved via linearization.
            If this flag is ``False`` (default) the spectral transformation
            is applied to the linearization (handled by the eigensolver),
            otherwise it is applied to the original problem.

        See Also
        --------
        setTransform, slepc.STGetTransform
        """
        cdef PetscBool sval = PETSC_FALSE
        CHKERR( STGetTransform(self.st, &sval) )
        return toBool(sval)

    def setMatMode(self, mode: MatMode) -> None:
        """
        Set a flag related to management of transformed matrices.

        Logically collective.

        The flag indicates how the transformed matrices are being
        stored in the spectral transformation.

        Parameters
        ----------
        mode
            The mode flag.

        Notes
        -----
        By default (`ST.MatMode.COPY`), a copy of matrix :math:`A` is made
        and then this copy is modified explicitly, e.g.,
        :math:`A \leftarrow (A - \sigma B)`.

        With `ST.MatMode.INPLACE`, the original matrix :math:`A` is modified at
        `setUp()` and reverted at the end of the computations. With respect to
        the previous one, this mode avoids a copy of matrix :math:`A`. However,
        a backdraw is that the recovered matrix might be slightly different
        from the original one (due to roundoff).

        With `ST.MatMode.SHELL`, the solver works with an implicit shell matrix
        that represents the shifted matrix. This mode is the most efficient in
        creating the transformed matrix but it places serious limitations to the
        linear solves performed in each iteration of the eigensolver
        (typically, only iterative solvers with Jacobi preconditioning can be
        used).

        In the two first modes the efficiency of this computation can be
        controlled with `setMatStructure()`.

        See Also
        --------
        setMatrices, setMatStructure, getMatMode, slepc.STSetMatMode
        """
        cdef SlepcSTMatMode val = mode
        CHKERR( STSetMatMode(self.st, val) )

    def getMatMode(self) -> MatMode:
        """
        Get a flag that indicates how the matrix is being shifted.

        Not collective.

        Get a flag that indicates how the matrix is being shifted in
        the shift-and-invert and Cayley spectral transformations.

        Returns
        -------
        MatMode
            The mode flag.

        See Also
        --------
        setMatMode, slepc.STGetMatMode
        """
        cdef SlepcSTMatMode val = ST_MATMODE_INPLACE
        CHKERR( STGetMatMode(self.st, &val) )
        return val

    def setMatrices(self, operators: list[Mat]) -> None:
        """
        Set the matrices associated with the eigenvalue problem.

        Collective.

        Parameters
        ----------
        operators
            The matrices associated with the eigensystem.

        Notes
        -----
        It must be called before `setUp()`. If it is called again after
        `setUp()` then the `ST` object is reset.

        In standard eigenproblems only one matrix is passed, while in
        generalized problems two matrices are provided. The number of
        matrices is larger in polynomial eigenproblems.

        In normal usage, matrices are provided via the corresponding
        `EPS` of `PEP` interface function.

        See Also
        --------
        getMatrices, setUp, reset, slepc.STSetMatrices
        """
        operators = tuple(operators)
        cdef PetscMat *mats = NULL
        cdef Py_ssize_t k=0, n = len(operators)
        cdef tmp = allocate(<size_t>n*sizeof(PetscMat),<void**>&mats)
        for k from 0 <= k < n: mats[k] = (<Mat?>operators[k]).mat
        CHKERR( STSetMatrices(self.st, <PetscInt>n, mats) )

    def getMatrices(self) -> list[Mat]:
        """
        Get the matrices associated with the eigenvalue problem.

        Collective.

        Returns
        -------
        list of petsc4py.PETSc.Mat
            The matrices associated with the eigensystem.

        See Also
        --------
        setMatrices, slepc.STGetNumMatrices, slepc.STGetMatrix
        """
        cdef Mat A
        cdef PetscMat mat = NULL
        cdef PetscInt k=0, n=0
        CHKERR( STGetNumMatrices(self.st, &n) )
        cdef object operators = []
        for k from 0 <= k < n:
            CHKERR( STGetMatrix(self.st, k, &mat) )
            A = Mat(); A.mat = mat; CHKERR( PetscINCREF(A.obj) )
            operators.append(A)
        return tuple(operators)

    def setMatStructure(self, structure: petsc4py.PETSc.Mat.Structure) -> None:
        """
        Set the matrix structure attribute.

        Logically collective.

        Set an internal `petsc4py.PETSc.Mat.Structure` attribute to indicate
        which is the relation of the sparsity pattern of all the `ST` matrices.

        Parameters
        ----------
        structure
            The matrix structure specification.

        Notes
        -----
        By default, the sparsity patterns are assumed to be
        different. If the patterns are equal or a subset then it is
        recommended to set this attribute for efficiency reasons (in
        particular, for internal ``Mat.axpy()`` operations).

        This function has no effect in the case of standard eigenproblems.

        In case of polynomial eigenproblems, the flag applies to all
        matrices relative to the first one.

        See Also
        --------
        getMatStructure, setMatrices, slepc.STSetMatStructure
        """
        cdef PetscMatStructure val = matstructure(structure)
        CHKERR( STSetMatStructure(self.st, val) )

    def getMatStructure(self) -> petsc4py.PETSc.Mat.Structure:
        """
        Get the internal matrix structure attribute.

        Not collective.

        Get the internal `petsc4py.PETSc.Mat.Structure` attribute to
        indicate which is the relation of the sparsity pattern of the
        matrices.

        Returns
        -------
        petsc4py.PETSc.Mat.Structure
            The structure flag.

        See Also
        --------
        setMatStructure, slepc.STGetMatStructure
        """
        cdef PetscMatStructure val
        CHKERR( STGetMatStructure(self.st, &val) )
        return val

    def setKSP(self, KSP ksp) -> None:
        """
        Set the ``KSP`` object associated with the spectral transformation.

        Collective.

        Parameters
        ----------
        ksp
            The linear solver object.

        See Also
        --------
        getKSP, slepc.STSetKSP
        """
        CHKERR( STSetKSP(self.st, ksp.ksp) )

    def getKSP(self) -> KSP:
        """
        Get the ``KSP`` object associated with the spectral transformation.

        Collective.

        Returns
        -------
        `petsc4py.PETSc.KSP`
            The linear solver object.

        See Also
        --------
        setKSP, slepc.STGetKSP
        """
        cdef KSP ksp = KSP()
        CHKERR( STGetKSP(self.st, &ksp.ksp) )
        CHKERR( PetscINCREF(ksp.obj) )
        return ksp

    def setPreconditionerMat(self, Mat P = None) -> None:
        """
        Set the matrix to be used to build the preconditioner.

        Collective.

        Parameters
        ----------
        P
            The matrix that will be used in constructing the preconditioner.

        Notes
        -----
        This matrix will be passed to the internal ``KSP`` object (via the last
        argument of ``KSP.setOperators()``) as the matrix to be
        used when constructing the preconditioner. If no matrix is set then
        :math:`A-\sigma B` will be used to build the preconditioner, being
        :math:`\sigma` the value set by `setShift()`.

        More precisely, this is relevant for spectral transformations that
        represent a rational matrix function, and use a ``KSP`` object for the
        denominator. It includes also the `PRECOND` case. If the user has a
        good approximation to matrix that can be used to build a cheap
        preconditioner, it can be passed with this function. Note that it
        affects only the ``Pmat`` argument of ``KSP.setOperators()``,
        not the ``Amat`` argument.

        If a preconditioner matrix is set, the default is to use an iterative
        ``KSP`` rather than a direct method.

        An alternative to pass an approximation of :math:`A-\sigma B` with this
        function is to provide approximations of :math:`A` and :math:`B` via
        `setSplitPreconditioner()`. The difference is that when :math:`\sigma`
        changes the preconditioner is recomputed.

        A call with no matrix argument will remove a previously set matrix.

        See Also
        --------
        getPreconditionerMat, slepc.STSetPreconditionerMat
        """
        cdef PetscMat Pmat = P.mat if P is not None else <PetscMat>NULL
        CHKERR( STSetPreconditionerMat(self.st, Pmat) )

    def getPreconditionerMat(self) -> Mat:
        """
        Get the matrix previously set by `setPreconditionerMat()`.

        Not collective.

        Returns
        -------
        petsc4py.PETSc.Mat
            The matrix that will be used in constructing the preconditioner.

        See Also
        --------
        setPreconditionerMat, slepc.STGetPreconditionerMat
        """
        cdef Mat P = Mat()
        CHKERR( STGetPreconditionerMat(self.st, &P.mat) )
        CHKERR( PetscINCREF(P.obj) )
        return P

    def setSplitPreconditioner(self, operators: list[petsc4py.PETSc.Mat], structure: petsc4py.PETSc.Mat.
Structure | None = None) -> None:
        """
        Set the matrices to be used to build the preconditioner.

        Collective.

        Parameters
        ----------
        operators
            The matrices associated with the preconditioner.
        structure
            The matrix structure specification.

        Notes
        -----
        The number of matrices passed here must be the same as in `setMatrices()`.

        For linear eigenproblems, the preconditioner matrix is computed as
        :math:`P(\sigma) = A_0-\sigma B_0`, where :math:`A_0,B_0` are
        approximations of :math:`A,B` (the eigenproblem matrices) provided via the
        ``operators`` argument in this function. Compared to `setPreconditionerMat()`,
        this function allows setting a preconditioner in a way that is independent
        of the shift :math:`\sigma`. Whenever the value of :math:`\sigma` changes
        the preconditioner is recomputed.

        Similarly, for polynomial eigenproblems the matrix for the preconditioner
        is expressed as :math:`P(\sigma) = \sum_i P_i \phi_i(\sigma)`, for
        :math:`i=1,\dots,n`, where :math:`P_i` are given in ``operators`` and the
        :math:`\phi_i`'s are the polynomial basis functions.

        The ``structure`` flag provides information about the relative nonzero
        pattern of the ``operators`` matrices, in the same way as in
        `setMatStructure()`.

        See Also
        --------
        getSplitPreconditioner, setPreconditionerMat, slepc.STSetSplitPreconditioner
        """
        operators = tuple(operators)
        cdef PetscMatStructure cstructure = matstructure(structure)
        cdef PetscMat *mats = NULL
        cdef Py_ssize_t k=0, n = len(operators)
        cdef tmp = allocate(<size_t>n*sizeof(PetscMat),<void**>&mats)
        for k from 0 <= k < n: mats[k] = (<Mat?>operators[k]).mat
        CHKERR( STSetSplitPreconditioner(self.st, <PetscInt>n, mats, cstructure) )

    def getSplitPreconditioner(self) -> tuple[list[petsc4py.PETSc.Mat], petsc4py.PETSc.Mat.Structure]:
        """
        Get the matrices to be used to build the preconditioner.

        Not collective.

        Returns
        -------
        list of petsc4py.PETSc.Mat
            The list of matrices associated with the preconditioner.
        petsc4py.PETSc.Mat.Structure
            The structure flag.

        See Also
        --------
        slepc.STGetSplitPreconditionerInfo, slepc.STGetSplitPreconditionerTerm
        """
        cdef PetscInt k=0,n=0
        cdef PetscMatStructure cstructure
        cdef PetscMat mat = NULL
        CHKERR( STGetSplitPreconditionerInfo(self.st, &n, &cstructure) )
        cdef object operators = []
        for k from 0 <= k < n:
            CHKERR( STGetSplitPreconditionerTerm(self.st, k, &mat) )
            A = Mat(); A.mat = mat; CHKERR( PetscINCREF(A.obj) )
            operators.append(A)
        return tuple(operators, cstructure)

    #

    def setUp(self) -> None:
        """
        Prepare for the use of a spectral transformation.

        Collective.

        See Also
        --------
        apply, slepc.STSetUp
        """
        CHKERR( STSetUp(self.st) )

    def apply(self, Vec x, Vec y) -> None:
        """
        Apply the spectral transformation operator to a vector.

        Collective.

        Apply the spectral transformation operator to a vector, for instance
        :math:`y=(A-\sigma B)^{-1}Bx` in the case of the shift-and-invert
        transformation and generalized eigenproblem.

        Parameters
        ----------
        x
            The input vector.
        y
            The result vector.

        See Also
        --------
        applyTranspose, applyHermitianTranspose, applyMat, slepc.STApply
        """
        CHKERR( STApply(self.st, x.vec, y.vec) )

    def applyTranspose(self, Vec x, Vec y) -> None:
        """
        Apply the transpose of the operator to a vector.

        Collective.

        Apply the transpose of the operator to a vector, for instance
        :math:`y=B^T(A-\sigma B)^{-T}x` in the case of the shift-and-invert
        transformation and generalized eigenproblem.

        Parameters
        ----------
        x
            The input vector.
        y
            The result vector.

        See Also
        --------
        apply, applyHermitianTranspose, slepc.STApplyTranspose
        """
        CHKERR( STApplyTranspose(self.st, x.vec, y.vec) )

    def applyHermitianTranspose(self, Vec x, Vec y) -> None:
        """
        Apply the Hermitian-transpose of the operator to a vector.

        Collective.

        Apply the Hermitian-transpose of the operator to a vector, for instance
        :math:`y=B^*(A - \sigma B)^{-*}x` in the case of the shift-and-invert
        transformation and generalized eigenproblem.

        Parameters
        ----------
        x
            The input vector.
        y
            The result vector.

        See Also
        --------
        apply, applyTranspose, slepc.STApplyHermitianTranspose
        """
        CHKERR( STApplyHermitianTranspose(self.st, x.vec, y.vec) )

    def applyMat(self, Mat X, Mat Y) -> None:
        """
        Apply the spectral transformation operator to a matrix.

        Collective.

        Apply the spectral transformation operator to a matrix, for instance
        :math:`Y=(A-\sigma B)^{-1}BX` in the case of the shift-and-invert
        transformation and generalized eigenproblem.

        Parameters
        ----------
        X
            The input matrix.
        Y
            The result matrix.

        See Also
        --------
        apply, slepc.STApplyMat
        """
        CHKERR( STApplyMat(self.st, X.mat, Y.mat) )

    def getOperator(self) -> Mat:
        """
        Get a shell matrix that represents the operator of the spectral transformation.

        Collective.

        Returns
        -------
        petsc4py.PETSc.Mat
            Operator matrix.

        Notes
        -----
        The operator is defined in linear eigenproblems only, not in
        polynomial ones, so the call will fail if more than 2 matrices
        were passed in `setMatrices()`.

        The returned shell matrix is essentially a wrapper to the `apply()`
        and `applyTranspose()` operations. The operator can often be expressed as

        .. math::

           Op = D K^{-1} M D^{-1}

        where :math:`D` is the balancing matrix, and :math:`M` and :math:`K` are
        two matrices corresponding to the numerator and denominator for spectral
        transformations that represent a rational matrix function.

        The preconditioner matrix :math:`K` typically depends on the value of the
        shift, and its inverse is handled via an internal ``KSP`` object. Normal
        usage does not require explicitly calling `getOperator()`, but it can be
        used to force the creation of :math:`K` and :math:`M`, and then :math:`K`
        is passed to the ``KSP``. This is useful for setting options associated
        with the ``PCFactor`` (to set MUMPS options, for instance).

        The returned matrix must NOT be destroyed by the user. Instead, when no
        longer needed it must be returned with `restoreOperator()`. In particular,
        this is required before modifying the `ST` matrices or the shift.

        See Also
        --------
        apply, setMatrices, setShift, restoreOperator, slepc.STGetOperator
        """
        cdef Mat op = Mat()
        CHKERR( STGetOperator(self.st, &op.mat) )
        CHKERR( PetscINCREF(op.obj) )
        return op

    def restoreOperator(self, Mat op) -> None:
        """
        Restore the previously seized operator matrix.

        Logically collective.

        Parameters
        ----------
        op
            Operator matrix previously obtained with `getOperator()`.

        See Also
        --------
        getOperator, slepc.STRestoreOperator
        """
        CHKERR( PetscObjectDereference(<PetscObject>op.mat) )
        CHKERR( STRestoreOperator(self.st, &op.mat) )

    #

    def setCayleyAntishift(self, mu: Scalar) -> None:
        """
        Set the value of the anti-shift for the Cayley spectral transformation.

        Logically collective.

        Parameters
        ----------
        mu
            The anti-shift.

        Notes
        -----
        In the generalized Cayley transform, the operator can be expressed as
        :math:`(A - \sigma B)^{-1}(A + \mu B)`. This function sets the value
        of :math:`mu`.  Use `setShift()` for setting :math:`\sigma`.

        See Also
        --------
        setShift, getCayleyAntishift, slepc.STCayleySetAntishift
        """
        cdef PetscScalar sval = asScalar(mu)
        CHKERR( STCayleySetAntishift(self.st, sval) )

    def getCayleyAntishift(self) -> Scalar:
        """
        Get the value of the anti-shift for the Cayley spectral transformation.

        Not collective.

        Returns
        -------
        Scalar
            The anti-shift.

        See Also
        --------
        setCayleyAntishift, slepc.STCayleyGetAntishift
        """
        cdef PetscScalar sval = 0
        CHKERR( STCayleyGetAntishift(self.st, &sval) )
        return toScalar(sval)

    def setFilterType(self, filter_type: FilterType) -> None:
        """
        Set the method to be used to build the polynomial filter.

        Logically collective.

        Parameter
        ---------
        filter_type
            The type of filter.

        See Also
        --------
        getFilterType, slepc.STFilterSetType
        """
        cdef SlepcSTFilterType val = filter_type
        CHKERR( STFilterSetType(self.st, val) )

    def getFilterType(self) -> FilterType:
        """
        Get the method to be used to build the polynomial filter.

        Not collective.

        Returns
        -------
        FilterType
            The type of filter.

        See Also
        --------
        setFilterType, slepc.STFilterGetType
        """
        cdef SlepcSTFilterType val = ST_FILTER_FILTLAN
        CHKERR( STFilterGetType(self.st, &val) )
        return val

    def setFilterInterval(self, inta: float, intb: float) -> None:
        """
        Set the interval containing the desired eigenvalues.

        Logically collective.

        Parameters
        ----------
        inta
            The left end of the interval.
        intb
            The right end of the interval.

        Notes
        -----
        The filter will be configured to emphasize eigenvalues contained
        in the given interval, and damp out eigenvalues outside it. If the
        interval is open, then the filter is low- or high-pass, otherwise
        it is mid-pass.

        Common usage is to set the interval in `EPS` with `EPS.setInterval()`.

        The interval must be contained within the numerical range of the
        matrix, see `setFilterRange()`.

        See Also
        --------
        getFilterInterval, setFilterRange, slepc.STFilterSetInterval
        """
        cdef PetscReal rval1 = asReal(inta)
        cdef PetscReal rval2 = asReal(intb)
        CHKERR( STFilterSetInterval(self.st, rval1, rval2) )

    def getFilterInterval(self) -> tuple[float, float]:
        """
        Get the interval containing the desired eigenvalues.

        Not collective.

        Returns
        -------
        inta: float
            The left end of the interval.
        intb: float
            The right end of the interval.

        See Also
        --------
        setFilterInterval, slepc.STFilterGetInterval
        """
        cdef PetscReal inta = 0
        cdef PetscReal intb = 0
        CHKERR( STFilterGetInterval(self.st, &inta, &intb) )
        return (toReal(inta), toReal(intb))

    def setFilterRange(self, left: float, right: float) -> None:
        """
        Set the numerical range (or field of values) of the matrix.

        Logically collective.

        Set the numerical range (or field of values) of the matrix, that is,
        the interval containing all eigenvalues.

        Parameters
        ----------
        left
            The left end of the spectral range.
        right
            The right end of the spectral range.

        Notes
        -----
        The filter will be most effective if the numerical range is tight,
        that is, ``left`` and ``right`` are good approximations to the
        leftmost and rightmost eigenvalues, respectively.

        See Also
        --------
        setFilterInterval, getFilterRange, slepc.STFilterSetRange
        """
        cdef PetscReal rval1 = asReal(left)
        cdef PetscReal rval2 = asReal(right)
        CHKERR( STFilterSetRange(self.st, rval1, rval2) )

    def getFilterRange(self) -> tuple[float, float]:
        """
        Get the interval containing all eigenvalues.

        Not collective.

        Returns
        -------
        left: float
            The left end of the spectral range.
        right: float
            The right end of the spectral range.

        See Also
        --------
        getFilterInterval, slepc.STFilterGetRange
        """
        cdef PetscReal left = 0
        cdef PetscReal right = 0
        CHKERR( STFilterGetRange(self.st, &left, &right) )
        return (toReal(left), toReal(right))

    def setFilterDegree(self, deg: int) -> None:
        """
        Set the degree of the filter polynomial.

        Logically collective.

        Parameters
        ----------
        deg
            The polynomial degree.

        See Also
        --------
        getFilterDegree, slepc.STFilterSetDegree
        """
        cdef PetscInt val = asInt(deg)
        CHKERR( STFilterSetDegree(self.st, val) )

    def getFilterDegree(self) -> int:
        """
        Get the degree of the filter polynomial.

        Not collective.

        Returns
        -------
        int
            The polynomial degree.

        See Also
        --------
        setFilterDegree, slepc.STFilterGetDegree
        """
        cdef PetscInt val = 0
        CHKERR( STFilterGetDegree(self.st, &val) )
        return toInt(val)

    def setFilterDamping(self, damping: FilterDamping) -> None:
        """
        Set the type of damping to be used in the polynomial filter.

        Logically collective.

        Parameter
        ---------
        damping
            The type of damping.

        Notes
        -----
        Only used in `FilterType.CHEBYSHEV` filters.

        See Also
        --------
        getFilterDamping, slepc.STFilterSetDamping
        """
        cdef SlepcSTFilterDamping val = damping
        CHKERR( STFilterSetDamping(self.st, val) )

    def getFilterDamping(self) -> FilterDamping:
        """
        Get the type of damping used in the polynomial filter.

        Not collective.

        Returns
        -------
        FilterDamping
            The type of damping.

        See Also
        --------
        setFilterDamping, slepc.STFilterGetDamping
        """
        cdef SlepcSTFilterDamping val = ST_FILTER_DAMPING_NONE
        CHKERR( STFilterGetDamping(self.st, &val) )
        return val

    #

    property shift:
        """Value of the shift."""
        def __get__(self) -> float:
            return self.getShift()
        def __set__(self, value):
            self.setShift(value)

    property transform:
        """If the transformed matrices are computed."""
        def __get__(self) -> bool:
            return self.getTransform()
        def __set__(self, value):
            self.setTransform(value)

    property mat_mode:
        """How the transformed matrices are being stored in the ST."""
        def __get__(self) -> STMatMode:
            return self.getMatMode()
        def __set__(self, value):
            self.setMatMode(value)

    property mat_structure:
        """Relation of the sparsity pattern of all ST matrices."""
        def __get__(self) -> MatStructure:
            return self.getMatStructure()
        def __set__(self, value):
            self.setMatStructure(value)

    property ksp:
        """KSP object associated with the spectral transformation."""
        def __get__(self) -> KSP:
            return self.getKSP()
        def __set__(self, value):
            self.setKSP(value)

# -----------------------------------------------------------------------------

del STType
del STMatMode
del STFilterType
del STFilterDamping

# -----------------------------------------------------------------------------
