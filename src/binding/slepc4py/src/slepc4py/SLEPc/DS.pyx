# -----------------------------------------------------------------------------

class DSType(object):
    """DS type."""
    HEP     = S_(DSHEP)
    NHEP    = S_(DSNHEP)
    GHEP    = S_(DSGHEP)
    GHIEP   = S_(DSGHIEP)
    GNHEP   = S_(DSGNHEP)
    NHEPTS  = S_(DSNHEPTS)
    SVD     = S_(DSSVD)
    HSVD    = S_(DSHSVD)
    GSVD    = S_(DSGSVD)
    PEP     = S_(DSPEP)
    NEP     = S_(DSNEP)

class DSStateType(object):
    """
    DS state types.

    - `RAW`:          Not processed yet.
    - `INTERMEDIATE`: Reduced to Hessenberg or tridiagonal form (or equivalent).
    - `CONDENSED`:    Reduced to Schur or diagonal form (or equivalent).
    - `TRUNCATED`:    Condensed form truncated to a smaller size.
    """
    RAW          = DS_STATE_RAW
    INTERMEDIATE = DS_STATE_INTERMEDIATE
    CONDENSED    = DS_STATE_CONDENSED
    TRUNCATED    = DS_STATE_TRUNCATED

class DSMatType(object):
    """
    To refer to one of the matrices stored internally in DS.

    - `A`:  first matrix of eigenproblem/singular value problem.
    - `B`:  second matrix of a generalized eigenproblem.
    - `C`:  third matrix of a quadratic eigenproblem.
    - `T`:  tridiagonal matrix.
    - `D`:  diagonal matrix.
    - `Q`:  orthogonal matrix of (right) Schur vectors.
    - `Z`:  orthogonal matrix of left Schur vectors.
    - `X`:  right eigenvectors.
    - `Y`:  left eigenvectors.
    - `U`:  left singular vectors.
    - `V`:  right singular vectors.
    - `W`:  workspace matrix.
    """
    A  = DS_MAT_A
    B  = DS_MAT_B
    C  = DS_MAT_C
    T  = DS_MAT_T
    D  = DS_MAT_D
    Q  = DS_MAT_Q
    Z  = DS_MAT_Z
    X  = DS_MAT_X
    Y  = DS_MAT_Y
    U  = DS_MAT_U
    V  = DS_MAT_V
    W  = DS_MAT_W

class DSParallelType(object):
    """
    DS parallel types.

    - `REDUNDANT`:    Every process performs the computation redundantly.
    - `SYNCHRONIZED`: The first process sends the result to the rest.
    - `DISTRIBUTED`:  Used in some cases to distribute the computation among processes.
    """
    REDUNDANT    = DS_PARALLEL_REDUNDANT
    SYNCHRONIZED = DS_PARALLEL_SYNCHRONIZED
    DISTRIBUTED  = DS_PARALLEL_DISTRIBUTED

# -----------------------------------------------------------------------------

cdef class DS(Object):

    """DS."""

    Type         = DSType
    StateType    = DSStateType
    MatType      = DSMatType
    ParallelType = DSParallelType

    def __cinit__(self):
        self.obj = <PetscObject*> &self.ds
        self.ds = NULL

    def view(self, Viewer viewer=None) -> None:
        """
        Prints the DS data structure.

        Parameters
        ----------
        viewer
            Visualization context; if not provided, the standard
            output is used.
        """
        cdef PetscViewer vwr = def_Viewer(viewer)
        CHKERR( DSView(self.ds, vwr) )

    def destroy(self) -> Self:
        """Destroys the DS object."""
        CHKERR( DSDestroy(&self.ds) )
        self.ds = NULL
        return self

    def reset(self) -> None:
        """Resets the DS object."""
        CHKERR( DSReset(self.ds) )

    def create(self, comm: Comm | None = None) -> Self:
        """
        Creates the DS object.

        Parameters
        ----------
        comm
            MPI communicator; if not provided, it defaults to all processes.
        """
        cdef MPI_Comm ccomm = def_Comm(comm, SLEPC_COMM_DEFAULT())
        cdef SlepcDS newds = NULL
        CHKERR( DSCreate(ccomm, &newds) )
        CHKERR( SlepcCLEAR(self.obj) ); self.ds = newds
        return self

    def setType(self, ds_type: Type | str) -> None:
        """
        Selects the type for the DS object.

        Parameters
        ----------
        ds_type
            The direct solver type to be used.
        """
        cdef SlepcDSType cval = NULL
        ds_type = str2bytes(ds_type, &cval)
        CHKERR( DSSetType(self.ds, cval) )

    def getType(self) -> str:
        """
        Gets the DS type of this object.

        Returns
        -------
        str
            The direct solver type currently being used.
        """
        cdef SlepcDSType ds_type = NULL
        CHKERR( DSGetType(self.ds, &ds_type) )
        return bytes2str(ds_type)

    def setOptionsPrefix(self, prefix: str | None = None) -> None:
        """
        Sets the prefix used for searching for all DS options in the
        database.

        Parameters
        ----------
        prefix
            The prefix string to prepend to all DS option requests.

        Notes
        -----
        A hyphen (``-``) must NOT be given at the beginning of the
        prefix name.  The first character of all runtime options is
        AUTOMATICALLY the hyphen.
        """
        cdef const char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( DSSetOptionsPrefix(self.ds, cval) )

    def appendOptionsPrefix(self, prefix: str | None = None) -> None:
        """
        Appends to the prefix used for searching for all DS options
        in the database.

        Parameters
        ----------
        prefix
            The prefix string to prepend to all DS option requests.
        """
        cdef const char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( DSAppendOptionsPrefix(self.ds, cval) )

    def getOptionsPrefix(self) -> str:
        """
        Gets the prefix used for searching for all DS options in the
        database.

        Returns
        -------
        str
            The prefix string set for this DS object.
        """
        cdef const char *prefix = NULL
        CHKERR( DSGetOptionsPrefix(self.ds, &prefix) )
        return bytes2str(prefix)

    def setFromOptions(self) -> None:
        """
        Sets DS options from the options database.

        Notes
        -----
        To see all options, run your program with the ``-help``
        option.
        """
        CHKERR( DSSetFromOptions(self.ds) )

    def duplicate(self) -> DS:
        """Duplicate the DS object with the same type and dimensions."""
        cdef DS ds = type(self)()
        CHKERR( DSDuplicate(self.ds, &ds.ds) )
        return ds

    #

    def allocate(self, ld: int) -> None:
        """
        Allocates memory for internal storage or matrices in DS.

        Parameters
        ----------
        ld
            Leading dimension (maximum allowed dimension for the
            matrices, including the extra row if present).
        """
        cdef PetscInt val = asInt(ld)
        CHKERR( DSAllocate(self.ds, val) )

    def getLeadingDimension(self) -> int:
        """
        Returns the leading dimension of the allocated matrices.

        Returns
        -------
        int
            Leading dimension (maximum allowed dimension for the matrices).
        """
        cdef PetscInt val = 0
        CHKERR( DSGetLeadingDimension(self.ds, &val) )
        return toInt(val)

    def setState(self, state: StateType) -> None:
        """
        Change the state of the DS object.

        Parameters
        ----------
        state
            The new state.

        Notes
        -----
        The state indicates that the dense system is in an initial
        state (raw), in an intermediate state (such as tridiagonal,
        Hessenberg or Hessenberg-triangular), in a condensed state
        (such as diagonal, Schur or generalized Schur), or in a
        truncated state.

        This function is normally used to return to the raw state when
        the condensed structure is destroyed.
        """
        cdef SlepcDSStateType val = state
        CHKERR( DSSetState(self.ds, val) )

    def getState(self) -> StateType:
        """
        Returns the current state.

        Returns
        -------
        StateType
            The current state.
        """
        cdef SlepcDSStateType val = DS_STATE_RAW
        CHKERR( DSGetState(self.ds, &val) )
        return val

    def setParallel(self, pmode: ParallelType) -> None:
        """
        Selects the mode of operation in parallel runs.

        Parameters
        ----------
        pmode
            The parallel mode.
        """
        cdef SlepcDSParallelType val = pmode
        CHKERR( DSSetParallel(self.ds, val) )

    def getParallel(self) -> ParallelType:
        """
        Gets the mode of operation in parallel runs.

        Returns
        -------
        ParallelType
            The parallel mode.
        """
        cdef SlepcDSParallelType val = DS_PARALLEL_REDUNDANT
        CHKERR( DSGetParallel(self.ds, &val) )
        return val

    def setDimensions(self, n: int | None = None, l: int | None = None, k: int | None = None) -> None:
        """
        Resize the matrices in the DS object.

        Parameters
        ----------
        n
            The new size.
        l
            Number of locked (inactive) leading columns.
        k
            Intermediate dimension (e.g., position of arrow).

        Notes
        -----
        The internal arrays are not reallocated.
        """
        cdef PetscInt ival1 = PETSC_CURRENT
        cdef PetscInt ival2 = PETSC_CURRENT
        cdef PetscInt ival3 = PETSC_CURRENT
        if n is not None: ival1 = asInt(n)
        if l is not None: ival2 = asInt(l)
        if k is not None: ival3 = asInt(k)
        CHKERR( DSSetDimensions(self.ds, ival1, ival2, ival3) )

    def getDimensions(self) -> tuple[int, int, int, int]:
        """
        Returns the current dimensions.

        Returns
        -------
        n: int
            The new size.
        l: int
            Number of locked (inactive) leading columns.
        k: int
            Intermediate dimension (e.g., position of arrow).
        t: int
            Truncated length.
        """
        cdef PetscInt ival1 = 0
        cdef PetscInt ival2 = 0
        cdef PetscInt ival3 = 0
        cdef PetscInt ival4 = 0
        CHKERR( DSGetDimensions(self.ds, &ival1, &ival2, &ival3, &ival4) )
        return (toInt(ival1), toInt(ival2), toInt(ival3), toInt(ival4))

    def setBlockSize(self, bs: int) -> None:
        """
        Selects the block size.

        Parameters
        ----------
        bs
            The block size.
        """
        cdef PetscInt val = bs
        CHKERR( DSSetBlockSize(self.ds, val) )

    def getBlockSize(self) -> int:
        """
        Gets the block size.

        Returns
        -------
        int
            The block size.
        """
        cdef PetscInt val = 0
        CHKERR( DSGetBlockSize(self.ds, &val) )
        return val

    def setMethod(self, meth: int) -> None:
        """
        Selects the method to be used to solve the problem.

        Parameters
        ----------
        meth
            An index identifying the method.
        """
        cdef PetscInt val = meth
        CHKERR( DSSetMethod(self.ds, val) )

    def getMethod(self) -> int:
        """
        Gets the method currently used in the DS.

        Returns
        -------
        int
            Identifier of the method.
        """
        cdef PetscInt val = 0
        CHKERR( DSGetMethod(self.ds, &val) )
        return val

    def setCompact(self, comp: bool) -> None:
        """
        Switch to compact storage of matrices.

        Parameters
        ----------
        comp
            True means compact storage.

        Notes
        -----
        Compact storage is used in some `DS` types such as
        `DS.Type.HEP` when the matrix is tridiagonal. This flag
        can be used to indicate whether the user provides the
        matrix entries via the compact form (the tridiagonal
        `DS.MatType.T`) or the non-compact one (`DS.MatType.A`).

        The default is ``False``.
        """
        cdef PetscBool val = asBool(comp)
        CHKERR( DSSetCompact(self.ds, val) )

    def getCompact(self) -> bool:
        """
        Gets the compact storage flag.

        Returns
        -------
        bool
            The flag.
        """
        cdef PetscBool val = PETSC_FALSE
        CHKERR( DSGetCompact(self.ds, &val) )
        return toBool(val)

    def setExtraRow(self, ext: bool) -> None:
        """
        Sets a flag to indicate that the matrix has one extra row.

        Parameters
        ----------
        ext
            True if the matrix has extra row.

        Notes
        -----
        In Krylov methods it is useful that the matrix representing
        the direct solver has one extra row, i.e., has dimension
        (n+1)*n . If this flag is activated, all transformations
        applied to the right of the matrix also affect this additional
        row. In that case, (n+1) must be less or equal than the
        leading dimension.

        The default is ``False``.
        """
        cdef PetscBool val = asBool(ext)
        CHKERR( DSSetExtraRow(self.ds, val) )

    def getExtraRow(self) -> bool:
        """
        Gets the extra row flag.

        Returns
        -------
        bool
            The flag.
        """
        cdef PetscBool val = PETSC_FALSE
        CHKERR( DSGetExtraRow(self.ds, &val) )
        return toBool(val)

    def setRefined(self, ref: bool) -> None:
        """
        Sets a flag to indicate that refined vectors must be computed.

        Parameters
        ----------
        ref
            True if refined vectors must be used.

        Notes
        -----
        Normally the vectors returned in `DS.MatType.X` are eigenvectors
        of the projected matrix. With this flag activated, `vectors()`
        will return the right singular vector of the smallest singular
        value of matrix At-theta*I, where At is the extended (n+1)xn
        matrix and theta is the Ritz value. This is used in the
        refined Ritz approximation.

        The default is ``False``.
        """
        cdef PetscBool val = asBool(ref)
        CHKERR( DSSetRefined(self.ds, val) )

    def getRefined(self) -> bool:
        """
        Gets the refined vectors flag.

        Returns
        -------
        bool
            The flag.
        """
        cdef PetscBool val = PETSC_FALSE
        CHKERR( DSGetRefined(self.ds, &val) )
        return toBool(val)

    def truncate(self, n: int, trim: bool = False) -> None:
        """
        Truncates the system represented in the DS object.

        Parameters
        ----------
        n
            The new size.
        trim
            A flag to indicate if the factorization must be trimmed.
        """
        cdef PetscInt val = asInt(n)
        cdef PetscBool flg = asBool(trim)
        CHKERR( DSTruncate(self.ds, val, flg) )

    def updateExtraRow(self) -> None:
        """Performs all necessary operations so that the extra row gets up-to-date after a call to `solve()`."""
        CHKERR( DSUpdateExtraRow(self.ds) )

    def getMat(self, matname: MatType) -> Mat:
        """
        Returns the requested matrix as a sequential dense Mat object.

        Parameters
        ----------
        matname
            The requested matrix.

        Returns
        -------
        Mat
            The matrix.
        """
        cdef SlepcDSMatType mname = matname
        cdef Mat mat = Mat()
        CHKERR( DSGetMat(self.ds, mname, &mat.mat) )
        CHKERR( PetscINCREF(mat.obj) )
        return mat

    def restoreMat(self, matname: MatType, Mat mat) -> None:
        """
        Restore the previously seized matrix.

        Parameters
        ----------
        matname
            The selected matrix.
        mat
            The matrix previously obtained with `getMat()`.
        """
        cdef SlepcDSMatType mname = matname
        CHKERR( PetscObjectDereference(<PetscObject>mat.mat) )
        CHKERR( DSRestoreMat(self.ds, mname, &mat.mat) )

    def setIdentity(self, matname: MatType) -> None:
        """
        Copy the identity on the active part of a matrix.

        Parameters
        ----------
        matname
            The requested matrix.
        """
        cdef SlepcDSMatType mname = matname
        CHKERR( DSSetIdentity(self.ds, mname) )

    #

    def cond(self) -> float:
        """
        Compute the inf-norm condition number of the first matrix.

        Returns
        -------
        float
            Condition number.
        """
        cdef PetscReal rval = 0
        CHKERR( DSCond(self.ds, &rval) )
        return toReal(rval)

    def solve(self) -> ArrayReal | ArrayComplex:
        """
        Solve the problem.

        Returns
        -------
        ArrayReal | ArrayComplex
            Eigenvalues or singular values.
        """
        n = self.getDimensions()[0]
        cdef PetscScalar *eigr = NULL
        cdef PetscScalar *eigi = NULL
        cdef tmpr = allocate(<size_t>n*sizeof(PetscScalar), <void**>&eigr)
        cdef tmpi = allocate(<size_t>n*sizeof(PetscScalar), <void**>&eigi)
        CHKERR( DSSolve(self.ds, eigr, eigi) )
        cdef object kr = array_s(n, eigr)
        cdef object ki = array_s(n, eigi)
        if self.getType().upper() in ['HEP','GHEP','BSE','SVD','HSVD','GSVD']:
            return kr.real.copy()
        else:
            return kr+1j*ki

    def vectors(self, matname = MatType.X) -> None:
        """
        Compute vectors associated to the dense system such
        as eigenvectors.

        Parameters
        ----------
        matname: `DS.MatType` enumerate
           The matrix, used to indicate which vectors are required.
        """
        cdef SlepcDSMatType mname = matname
        CHKERR( DSVectors(self.ds, mname, NULL, NULL) )

    #

    def setSVDDimensions(self, m: int) -> None:
        """
        Sets the number of columns of a `DS` of type `SVD`.

        Parameters
        ----------
        m
            The number of columns.
        """
        cdef PetscInt val = asInt(m)
        CHKERR( DSSVDSetDimensions(self.ds, val) )

    def getSVDDimensions(self) -> int:
        """
        Gets the number of columns of a `DS` of type `SVD`.

        Returns
        -------
        int
            The number of columns.
        """
        cdef PetscInt val = 0
        CHKERR( DSSVDGetDimensions(self.ds, &val) )
        return toInt(val)

    def setHSVDDimensions(self, m: int) -> None:
        """
        Sets the number of columns of a `DS` of type `HSVD`.

        Parameters
        ----------
        m
            The number of columns.
        """
        cdef PetscInt val = asInt(m)
        CHKERR( DSHSVDSetDimensions(self.ds, val) )

    def getHSVDDimensions(self) -> int:
        """
        Gets the number of columns of a `DS` of type `HSVD`.

        Returns
        -------
        int
            The number of columns.
        """
        cdef PetscInt val = 0
        CHKERR( DSHSVDGetDimensions(self.ds, &val) )
        return toInt(val)

    def setGSVDDimensions(self, m: int, p: int) -> None:
        """
        Sets the number of columns and rows of a `DS` of type `GSVD`.

        Parameters
        ----------
        m
            The number of columns.
        p
            The number of rows for the second matrix.
        """
        cdef PetscInt val1 = asInt(m)
        cdef PetscInt val2 = asInt(p)
        CHKERR( DSGSVDSetDimensions(self.ds, val1, val2) )

    def getGSVDDimensions(self) -> tuple[int, int]:
        """
        Gets the number of columns and rows of a `DS` of type `GSVD`.

        Returns
        -------
        m: int
            The number of columns.
        p: int
            The number of rows for the second matrix.
        """
        cdef PetscInt val1 = 0
        cdef PetscInt val2 = 0
        CHKERR( DSGSVDGetDimensions(self.ds, &val1, &val2) )
        return (toInt(val1), toInt(val2))

    def setPEPDegree(self, deg: int) -> None:
        """
        Sets the polynomial degree of a `DS` of type `PEP`.

        Parameters
        ----------
        deg
            The polynomial degree.
        """
        cdef PetscInt val = asInt(deg)
        CHKERR( DSPEPSetDegree(self.ds, val) )

    def getPEPDegree(self) -> int:
        """
        Gets the polynomial degree of a `DS` of type `PEP`.

        Returns
        -------
        int
            The polynomial degree.
        """
        cdef PetscInt val = 0
        CHKERR( DSPEPGetDegree(self.ds, &val) )
        return toInt(val)

    def setPEPCoefficients(self, pbc: Sequence[float]) -> None:
        """
        Sets the polynomial basis coefficients of a `DS` of type `PEP`.

        Parameters
        ----------
        pbc
            Coefficients.
        """
        cdef PetscInt na = 0
        cdef PetscReal *a = NULL
        cdef object tmp1 = iarray_r(pbc, &na, &a)
        CHKERR( DSPEPSetCoefficients(self.ds, a) )

    def getPEPCoefficients(self) -> ArrayReal:
        """
        Gets the polynomial basis coefficients of a `DS` of type `PEP`.

        Returns
        -------
        ArrayReal
            Coefficients.
        """
        cdef PetscInt np = 0
        cdef PetscReal *coeff = NULL
        CHKERR( DSPEPGetDegree(self.ds, &np) )
        CHKERR( DSPEPGetCoefficients(self.ds, &coeff) )
        cdef object ocoeff = None
        try:
            ocoeff = array_r(3*(np+1), coeff)
        finally:
            CHKERR( PetscFree(coeff) )
        return ocoeff

    #

    property state:
        def __get__(self):
            return self.getState()
        def __set__(self, value):
            self.setState(value)

    property parallel:
        def __get__(self):
            return self.getParallel()
        def __set__(self, value):
            self.setParallel(value)

    property block_size:
        def __get__(self):
            return self.getBlockSize()
        def __set__(self, value):
            self.setBlockSize(value)

    property method:
        def __get__(self):
            return self.getMethod()
        def __set__(self, value):
            self.setMethod(value)

    property compact:
        def __get__(self):
            return self.getCompact()
        def __set__(self, value):
            self.setCompact(value)

    property extra_row:
        def __get__(self):
            return self.getExtraRow()
        def __set__(self, value):
            self.setExtraRow(value)

    property refined:
        def __get__(self):
            return self.getRefined()
        def __set__(self, value):
            self.setRefined(value)

# -----------------------------------------------------------------------------

del DSType
del DSStateType
del DSMatType
del DSParallelType

# -----------------------------------------------------------------------------
