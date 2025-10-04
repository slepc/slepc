# -----------------------------------------------------------------------------

class STType(object):
    """
    ST types.

    - `SHELL`:   User-defined.
    - `SHIFT`:   Shift from origin.
    - `SINVERT`: Shift-and-invert.
    - `CAYLEY`:  Cayley transform.
    - `PRECOND`: Preconditioner.
    - `FILTER`:  Polynomial filter.
    """
    SHELL   = S_(STSHELL)
    SHIFT   = S_(STSHIFT)
    SINVERT = S_(STSINVERT)
    CAYLEY  = S_(STCAYLEY)
    PRECOND = S_(STPRECOND)
    FILTER  = S_(STFILTER)

class STMatMode(object):
    """
    ST matrix mode.

    - `COPY`:    A working copy of the matrix is created.
    - `INPLACE`: The operation is computed in-place.
    - `SHELL`:   The matrix :math:`A - \sigma B` is handled as an
                 implicit matrix.
    """
    COPY    = ST_MATMODE_COPY
    INPLACE = ST_MATMODE_INPLACE
    SHELL   = ST_MATMODE_SHELL

class STFilterType(object):
    """
    ST filter type.

    - ``FILTLAN``:  An adapted implementation of the Filtered Lanczos Package.
    - ``CHEBYSEV``: A polynomial filter based on a truncated Chebyshev series.
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
    """
    NONE    = ST_FILTER_DAMPING_NONE
    JACKSON = ST_FILTER_DAMPING_JACKSON
    LANCZOS = ST_FILTER_DAMPING_LANCZOS
    FEJER   = ST_FILTER_DAMPING_FEJER

# -----------------------------------------------------------------------------

cdef class ST(Object):

    """ST."""

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
        """
        cdef PetscViewer vwr = def_Viewer(viewer)
        CHKERR( STView(self.st, vwr) )

    def destroy(self) -> Self:
        """
        Destroy the ST object.

        Collective.
        """
        CHKERR( STDestroy(&self.st) )
        self.st = NULL
        return self

    def reset(self) -> None:
        """
        Reset the ST object.

        Collective.
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
        See `ST.Type` for available methods. The default is
        `ST.Type.SHIFT` with a zero shift.  Normally, it is best to
        use `setFromOptions()` and then set the ST type from the
        options database rather than by using this routine.  Using the
        options database provides the user with maximum flexibility in
        evaluating the different available methods.
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
        """
        cdef const char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( STAppendOptionsPrefix(self.st, cval) )

    def setFromOptions(self) -> None:
        """
        Set ST options from the options database.

        Collective.

        This routine must be called before `setUp()` if the user is to be
        allowed to set the solver type.

        Notes
        -----
        To see all options, run your program with the -help option.
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
            If this flag is False (default) the spectral transformation
            is applied to the linearization (handled by the eigensolver),
            otherwise it is applied to the original problem.
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
            If this flag is False (default) the spectral transformation
            is applied to the linearization (handled by the eigensolver),
            otherwise it is applied to the original problem.
        """
        cdef PetscBool sval = PETSC_FALSE
        CHKERR( STGetTransform(self.st, &sval) )
        return toBool(sval)

    def setMatMode(self, mode: MatMode) -> None:
        """
        Set a flag to indicate how the matrix is being shifted.

        Logically collective.

        Set a flag to indicate how the matrix is being shifted in the
        shift-and-invert and Cayley spectral transformations.

        Parameters
        ----------
        mode
            The mode flag.

        Notes
        -----
        By default (`ST.MatMode.COPY`), a copy of matrix :math:`A` is made
        and then this copy is shifted explicitly, e.g.
        :math:`A \leftarrow (A - s B)`.

        With `ST.MatMode.INPLACE`, the original matrix :math:`A` is shifted at
        `setUp()` and unshifted at the end of the computations. With respect to
        the previous one, this mode avoids a copy of matrix :math:`A`. However,
        a backdraw is that the recovered matrix might be slightly different
        from the original one (due to roundoff).

        With `ST.MatMode.SHELL`, the solver works with an implicit shell matrix
        that represents the shifted matrix. This mode is the most efficient in
        creating the shifted matrix but it places serious limitations to the
        linear solves performed in each iteration of the eigensolver
        (typically, only iterative solvers with Jacobi preconditioning can be
        used).

        In the case of generalized problems, in the two first modes the matrix
        :math:`A - s B` has to be computed explicitly. The efficiency of
        this computation can be controlled with `setMatStructure()`.
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
        """
        operators = tuple(operators)
        cdef PetscMat *mats = NULL
        cdef Py_ssize_t k=0, n = len(operators)
        cdef tmp = allocate(<size_t>n*sizeof(PetscMat),<void**>&mats)
        for k from 0 <= k < n: mats[k] = (<Mat?>operators[k]).mat
        CHKERR( STSetMatrices(self.st, <PetscInt>n, mats) )

    def getMatrices(self) -> list[petsc4py.PETSc.Mat]:
        """
        Get the matrices associated with the eigenvalue problem.

        Collective.

        Returns
        -------
        list of petsc4py.PETSc.Mat
            The matrices associated with the eigensystem.
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
        Set an internal Mat.Structure attribute.

        Logically collective.

        Set an internal Mat.Structure attribute to indicate which is the
        relation of the sparsity pattern of the two matrices :math:`A` and
        :math:`B` constituting the generalized eigenvalue problem. This
        function has no effect in the case of standard eigenproblems.

        Parameters
        ----------
        structure
            Either same, different, or a subset of the non-zero
            sparsity pattern.

        Notes
        -----
        By default, the sparsity patterns are assumed to be
        different. If the patterns are equal or a subset then it is
        recommended to set this attribute for efficiency reasons (in
        particular, for internal *AXPY()* matrix operations).
        """
        cdef PetscMatStructure val = matstructure(structure)
        CHKERR( STSetMatStructure(self.st, val) )

    def getMatStructure(self) -> petsc4py.PETSc.Mat.Structure:
        """
        Get the internal Mat.Structure attribute.

        Not collective.

        Get the internal Mat.Structure attribute to indicate which is
        the relation of the sparsity pattern of the matrices.

        Returns
        -------
        petsc4py.PETSc.Mat.Structure
            The structure flag.
        """
        cdef PetscMatStructure val
        CHKERR( STGetMatStructure(self.st, &val) )
        return val

    def setKSP(self, KSP ksp: petsc4py.PETSc.KSP) -> None:
        """
        Set the ``KSP`` object associated with the spectral transformation.

        Collective.

        Parameters
        ----------
        `petsc4py.PETSc.KSP`
            The linear solver object.
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

        Notes
        -----
        On output, the internal value of `petsc4py.PETSc.KSP` can be ``NULL`` if the
        combination of eigenproblem type and selected transformation
        does not require to solve a linear system of equations.
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
        """
        cdef PetscMat Pmat = P.mat if P is not None else <PetscMat>NULL
        CHKERR( STSetPreconditionerMat(self.st, Pmat) )

    def getPreconditionerMat(self) -> petsc4py.PETSc.Mat:
        """
        Get the matrix previously set by setPreconditionerMat().

        Not collective.

        Returns
        -------
        petsc4py.PETSc.Mat
            The matrix that will be used in constructing the preconditioner.
        """
        cdef Mat P = Mat()
        CHKERR( STGetPreconditionerMat(self.st, &P.mat) )
        CHKERR( PetscINCREF(P.obj) )
        return P

    def setSplitPreconditioner(self, operators: list[petsc4py.PETSc.Mat], structure: petsc4py.PETSc.Mat.Structure | None = None) -> None:
        """
        Set the matrices to be used to build the preconditioner.

        Collective.

        Parameters
        ----------
        operators
            The matrices associated with the preconditioner.
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
        """
        CHKERR( STSetUp(self.st) )

    def apply(self, Vec x, Vec y) -> None:
        """
        Apply the spectral transformation operator to a vector.

        Collective.

        Apply the spectral transformation operator to a vector, for instance
        :math:`(A - s B)^{-1} B` in the case of the shift-and-invert
        transformation and generalized eigenproblem.

        Parameters
        ----------
        x
            The input vector.
        y
            The result vector.
        """
        CHKERR( STApply(self.st, x.vec, y.vec) )

    def applyTranspose(self, Vec x, Vec y) -> None:
        """
        Apply the transpose of the operator to a vector.

        Collective.

        Apply the transpose of the operator to a vector, for instance
        :math:`B^T(A - s B)^{-T}` in the case of the shift-and-invert
        transformation and generalized eigenproblem.

        Parameters
        ----------
        x
            The input vector.
        y
            The result vector.
        """
        CHKERR( STApplyTranspose(self.st, x.vec, y.vec) )

    def applyHermitianTranspose(self, Vec x, Vec y) -> None:
        """
        Apply the hermitian-transpose of the operator to a vector.

        Collective.

        Apply the hermitian-transpose of the operator to a vector, for instance
        :math:`B^H(A - s B)^{-H}` in the case of the shift-and-invert
        transformation and generalized eigenproblem.

        Parameters
        ----------
        x
            The input vector.
        y
            The result vector.
        """
        CHKERR( STApplyHermitianTranspose(self.st, x.vec, y.vec) )

    def applyMat(self, Mat x, Mat y) -> None:
        """
        Apply the spectral transformation operator to a matrix.

        Collective.

        Apply the spectral transformation operator to a matrix, for instance
        :math:`(A - s B)^{-1} B` in the case of the shift-and-invert
        transformation and generalized eigenproblem.

        Parameters
        ----------
        x
            The input matrix.
        y
            The result matrix.
        """
        CHKERR( STApplyMat(self.st, x.mat, y.mat) )

    def getOperator(self) -> petsc4py.PETSc.Mat:
        """
        Get a shell matrix that represents the operator of the spectral transformation.

        Collective.

        Returns
        -------
        petsc4py.PETSc.Mat
            Operator matrix.
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
            Operator matrix previously obtained with getOperator().
        """
        CHKERR( PetscObjectDereference(<PetscObject>op.mat) )
        CHKERR( STRestoreOperator(self.st, &op.mat) )

    #

    def setCayleyAntishift(self, tau: Scalar) -> None:
        """
        Set the value of the anti-shift for the Cayley spectral transformation.

        Logically collective.

        Parameters
        ----------
        tau
            The anti-shift.

        Notes
        -----
        In the generalized Cayley transform, the operator can be expressed as
        :math:`OP = inv(A - \sigma B) (A + tau B)`. This function sets
        the value of :math:`tau`.  Use `setShift()` for setting
        :math:`\sigma`.
        """
        cdef PetscScalar sval = asScalar(tau)
        CHKERR( STCayleySetAntishift(self.st, sval) )

    def getCayleyAntishift(self) -> Scalar:
        """
        Get the value of the anti-shift for the Cayley spectral transformation.

        Not collective.

        Returns
        -------
        Scalar
            The anti-shift.
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
        matrix, see `ST.setFilterRange()`.
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
            The left end of the interval.
        right
            The right end of the interval.

        Notes
        -----
        The filter will be most effective if the numerical range is tight,
        that is, left and right are good approximations to the leftmost and
        rightmost eigenvalues, respectively.
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
            The left end of the interval.
        right: float
            The right end of the interval.
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

# -----------------------------------------------------------------------------
