# -----------------------------------------------------------------------------

class BVType(object):
    """
    BV type
    """
    MAT        = S_(BVMAT)
    SVEC       = S_(BVSVEC)
    VECS       = S_(BVVECS)
    CONTIGUOUS = S_(BVCONTIGUOUS)
    TENSOR     = S_(BVTENSOR)

class BVOrthogType(object):
    """
    BV orthogonalization types

    - `CGS`: Classical Gram-Schmidt.
    - `MGS`: Modified Gram-Schmidt.
    """
    CGS = BV_ORTHOG_CGS
    MGS = BV_ORTHOG_MGS

class BVOrthogRefineType(object):
    """
    BV orthogonalization refinement types

    - `IFNEEDED`: Reorthogonalize if a criterion is satisfied.
    - `NEVER`:    Never reorthogonalize.
    - `ALWAYS`:   Always reorthogonalize.
    """
    IFNEEDED = BV_ORTHOG_REFINE_IFNEEDED
    NEVER    = BV_ORTHOG_REFINE_NEVER
    ALWAYS   = BV_ORTHOG_REFINE_ALWAYS

class BVOrthogBlockType(object):
    """
    BV block-orthogonalization types

    - `GS`:       Gram-Schmidt.
    - `CHOL`:     Cholesky.
    - `TSQR`:     Tall-skinny QR.
    - `TSQRCHOL`: Tall-skinny QR with Cholesky.
    - `SVQB`:     SVQB.
    """
    GS       = BV_ORTHOG_BLOCK_GS
    CHOL     = BV_ORTHOG_BLOCK_CHOL
    TSQR     = BV_ORTHOG_BLOCK_TSQR
    TSQRCHOL = BV_ORTHOG_BLOCK_TSQRCHOL
    SVQB     = BV_ORTHOG_BLOCK_SVQB

class BVMatMultType(object):
    """
    BV mat-mult types

    - `VECS`: Perform a matrix-vector multiply per each column.
    - `MAT`:  Carry out a Mat-Mat product with a dense matrix.
    """
    VECS     = BV_MATMULT_VECS
    MAT      = BV_MATMULT_MAT

class BVSVDMethod(object):
    """
    BV methods for computing the SVD

    - `REFINE`: Based on the SVD of the cross product matrix S'*S, with refinement.
    - `QR`:     Based on the SVD of the triangular factor of qr(S).
    - `QR_CAA`: Variant of QR intended for use in cammunication-avoiding Arnoldi.
    """
    REFINE   = BV_SVD_METHOD_REFINE
    QR       = BV_SVD_METHOD_QR
    QR_CAA   = BV_SVD_METHOD_QR_CAA

# -----------------------------------------------------------------------------

cdef class BV(Object):

    """
    BV
    """

    Type             = BVType
    OrthogType       = BVOrthogType
    OrthogRefineType = BVOrthogRefineType
    RefineType       = BVOrthogRefineType
    OrthogBlockType  = BVOrthogBlockType
    BlockType        = BVOrthogBlockType
    MatMultType      = BVMatMultType

    def __cinit__(self):
        self.obj = <PetscObject*> &self.bv
        self.bv = NULL

    # unary operations

    def __pos__(self):
        return bv_pos(self)

    def __neg__(self):
        return bv_neg(self)

    # inplace binary operations

    def __iadd__(self, other):
        return bv_iadd(self, other)

    def __isub__(self, other):
        return bv_isub(self, other)

    def __imul__(self, other):
        return bv_imul(self, other)

    def __idiv__(self, other):
        return bv_idiv(self, other)

    def __itruediv__(self, other):
        return bv_idiv(self, other)

    # binary operations

    def __add__(self, other):
        return bv_add(self, other)

    def __radd__(self, other):
        return bv_radd(self, other)

    def __sub__(self, other):
        return bv_sub(self, other)

    def __rsub__(self, other):
        return bv_rsub(self, other)

    def __mul__(self, other):
        return bv_mul(self, other)

    def __rmul__(self, other):
        return bv_rmul(self, other)

    def __div__(self, other):
        return bv_div(self, other)

    def __rdiv__(self, other):
        return bv_rdiv(self, other)

    def __truediv__(self, other):
        return bv_div(self, other)

    def __rtruediv__(self, other):
        return bv_rdiv(self, other)

    #

    def view(self, Viewer viewer=None) -> None:
        """
        Prints the BV data structure.

        Parameters
        ----------
        viewer
            Visualization context; if not provided, the standard
            output is used.
        """
        cdef PetscViewer vwr = def_Viewer(viewer)
        CHKERR( BVView(self.bv, vwr) )

    def destroy(self) -> Self:
        """
        Destroys the BV object.
        """
        CHKERR( BVDestroy(&self.bv) )
        self.bv = NULL
        return self

    def create(self, comm: Comm | None = None) -> Self:
        """
        Creates the BV object.

        Parameters
        ----------
        comm
            MPI communicator; if not provided, it defaults to all
            processes.
        """
        cdef MPI_Comm ccomm = def_Comm(comm, SLEPC_COMM_DEFAULT())
        cdef SlepcBV newbv = NULL
        CHKERR( BVCreate(ccomm, &newbv) )
        CHKERR( SlepcCLEAR(self.obj) ); self.bv = newbv
        return self

    def createFromMat(self, Mat A) -> Self:
        """
        Creates a basis vectors object from a dense Mat object.

        Parameters
        ----------
        A
            A dense tall-skinny matrix.
        """
        cdef SlepcBV newbv = NULL
        CHKERR( BVCreateFromMat(A.mat, &newbv) )
        CHKERR( SlepcCLEAR(self.obj) ); self.bv = newbv
        return self

    def createMat(self) -> Mat:
        """
        Creates a new Mat object of dense type and copies the contents of the
        BV object.

        Returns
        -------
        Mat
            The new matrix.
        """
        cdef Mat mat = Mat()
        CHKERR( BVCreateMat(self.bv, &mat.mat) )
        return mat

    def duplicate(self) -> BV:
        """
        Duplicate the BV object with the same type and dimensions.
        """
        cdef BV bv = type(self)()
        CHKERR( BVDuplicate(self.bv, &bv.bv) )
        return bv

    def duplicateResize(self, m: int) -> BV:
        """
        Creates a new BV object of the same type and dimensions as
        an existing one, but with possibly different number of columns.

        Parameters
        ----------
        m
            The number of columns.
        """
        cdef BV bv = type(self)()
        cdef PetscInt ival = asInt(m)
        CHKERR( BVDuplicateResize(self.bv, ival, &bv.bv) )
        return bv

    def copy(self, BV result=None) -> BV:
        """
        Copies a basis vector object into another one.

        Parameters
        ----------
        result
            The copy.
        """
        if result is None:
            result = type(self)()
        if result.bv == NULL:
            CHKERR( BVDuplicate(self.bv, &result.bv) )
        CHKERR( BVCopy(self.bv, result.bv) )
        return result

    def setType(self, bv_type: Type | str) -> None:
        """
        Selects the type for the BV object.

        Parameters
        ----------
        bv_type
            The inner product type to be used.
        """
        cdef SlepcBVType cval = NULL
        bv_type = str2bytes(bv_type, &cval)
        CHKERR( BVSetType(self.bv, cval) )

    def getType(self) -> str:
        """
        Gets the BV type of this object.

        Returns
        -------
        str
            The inner product type currently being used.
        """
        cdef SlepcBVType bv_type = NULL
        CHKERR( BVGetType(self.bv, &bv_type) )
        return bytes2str(bv_type)

    def setSizes(self, sizes: LayoutSizeSpec, m: int) -> None:
        """
        Sets the local and global sizes, and the number of columns.

        Parameters
        ----------
        sizes
            The global size ``N`` or a two-tuple ``(n, N)``
            with the local and global sizes.
        m
            The number of columns.

        Notes
        -----
        Either ``n`` or ``N`` (but not both) can be ``PETSc.DECIDE``
        or ``None`` to have it automatically set.
        """
        cdef PetscInt n=0, N=0
        cdef PetscInt ival = asInt(m)
        BV_Sizes(sizes, &n, &N)
        CHKERR( BVSetSizes(self.bv, n, N, ival) )

    def setSizesFromVec(self, Vec w, m: int) -> None:
        """
        Sets the local and global sizes, and the number of columns. Local and
        global sizes are specified indirectly by passing a template vector.

        Parameters
        ----------
        w
            The template vector.
        m
            The number of columns.
        """
        cdef PetscInt ival = asInt(m)
        CHKERR( BVSetSizesFromVec(self.bv, w.vec, ival) )

    def getSizes(self) -> tuple[LayoutSizeSpec, int]:
        """
        Returns the local and global sizes, and the number of columns.

        Returns
        -------
        (n, N): tuple[int, int]
            The local and global sizes
        m: int
            The number of columns.
        """
        cdef PetscInt n=0, N=0, m=0
        CHKERR( BVGetSizes(self.bv, &n, &N, &m) )
        return ((toInt(n), toInt(N)), toInt(m))

    def setLeadingDimension(self, ld: int) -> None:
        """
        Sets the leading dimension.

        Parameters
        ----------
        ld
            The leading dimension.
        """
        cdef PetscInt val = asInt(ld)
        CHKERR( BVSetLeadingDimension(self.bv, val) )

    def getLeadingDimension(self) -> int:
        """
        Gets the leading dimension.

        Returns
        -------
        int
            The leading dimension.
        """
        cdef PetscInt val = 0
        CHKERR( BVGetLeadingDimension(self.bv, &val) )
        return toInt(val)

    def setOptionsPrefix(self, prefix: str | None = None) -> None:
        """
        Sets the prefix used for searching for all BV options in the
        database.

        Parameters
        ----------
        prefix
            The prefix string to prepend to all BV option requests.

        Notes
        -----
        A hyphen (``-``) must NOT be given at the beginning of the
        prefix name.  The first character of all runtime options is
        AUTOMATICALLY the hyphen.
        """
        cdef const char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( BVSetOptionsPrefix(self.bv, cval) )

    def appendOptionsPrefix(self, prefix: str | None = None) -> None:
        """
        Appends to the prefix used for searching for all BV options
        in the database.

        Parameters
        ----------
        prefix
            The prefix string to prepend to all BV option requests.
        """
        cdef const char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( BVAppendOptionsPrefix(self.bv, cval) )

    def getOptionsPrefix(self) -> str:
        """
        Gets the prefix used for searching for all BV options in the
        database.

        Returns
        -------
        str
            The prefix string set for this BV object.
        """
        cdef const char *prefix = NULL
        CHKERR( BVGetOptionsPrefix(self.bv, &prefix) )
        return bytes2str(prefix)

    def setFromOptions(self) -> None:
        """
        Sets BV options from the options database.

        Notes
        -----
        To see all options, run your program with the ``-help``
        option.
        """
        CHKERR( BVSetFromOptions(self.bv) )

    #

    def getOrthogonalization(self) -> tuple[OrthogType, OrthogRefineType, float, OrthogBlockType]:
        """
        Gets the orthogonalization settings from the BV object.

        Returns
        -------
        type: OrthogType
            The type of orthogonalization technique.
        refine: OrthogRefineType
            The type of refinement.
        eta: float
            Parameter for selective refinement (used when the
            refinement type is `BV.OrthogRefineType.IFNEEDED`).
        block: OrthogBlockType
            The type of block orthogonalization .
        """
        cdef SlepcBVOrthogType val1 = BV_ORTHOG_CGS
        cdef SlepcBVOrthogRefineType val2 = BV_ORTHOG_REFINE_IFNEEDED
        cdef SlepcBVOrthogBlockType val3 = BV_ORTHOG_BLOCK_GS
        cdef PetscReal rval = PETSC_DEFAULT
        CHKERR( BVGetOrthogonalization(self.bv, &val1, &val2, &rval, &val3) )
        return (val1, val2, toReal(rval), val3)

    def setOrthogonalization(
        self,
        otype: OrthogType | None = None,
        refine: OrthogRefineType | None = None,
        eta: float | None = None,
        block: OrthogBlockType | None = None,
    ) -> None:
        """
        Specifies the method used for the orthogonalization of vectors
        (classical or modified Gram-Schmidt with or without refinement),
        and for the block-orthogonalization (simultaneous orthogonalization
        of a set of vectors).

        Parameters
        ----------
        otype
            The type of orthogonalization technique.
        refine
            The type of refinement.
        eta
            Parameter for selective refinement.
        block
            The type of block orthogonalization.

        Notes
        -----
        The default settings work well for most problems.

        The parameter `eta` should be a real value between ``0`` and
        ``1`` (or `DETERMINE`).  The value of `eta` is used only when
        the refinement type is `BV.OrthogRefineType.IFNEEDED`.

        When using several processors, `BV.OrthogType.MGS` is likely to
        result in bad scalability.

        If the method set for block orthogonalization is GS, then the
        computation is done column by column with the vector orthogonalization.
        """
        cdef SlepcBVOrthogType val1 = BV_ORTHOG_CGS
        cdef SlepcBVOrthogRefineType val2 = BV_ORTHOG_REFINE_IFNEEDED
        cdef SlepcBVOrthogBlockType val3 = BV_ORTHOG_BLOCK_GS
        cdef PetscReal rval = PETSC_CURRENT
        CHKERR( BVGetOrthogonalization(self.bv, &val1, &val2, NULL, &val3) )
        if otype  is not None: val1 = otype
        if refine is not None: val2 = refine
        if block  is not None: val3 = block
        if eta    is not None: rval = asReal(eta)
        CHKERR( BVSetOrthogonalization(self.bv, val1, val2, rval, val3) )

    def getMatMultMethod(self) -> MatMultType:
        """
        Gets the method used for the `matMult()` operation.

        Returns
        -------
        MatMultType
            The method for the `matMult()` operation.
        """
        cdef SlepcBVMatMultType val = BV_MATMULT_MAT
        CHKERR( BVGetMatMultMethod(self.bv, &val) )
        return val

    def setMatMultMethod(self, method: MatMultType) -> None:
        """
        Specifies the method used for the `matMult()` operation.

        Parameters
        ----------
        method
            The method for the `matMult()` operation.
        """
        cdef SlepcBVMatMultType val = method
        CHKERR( BVSetMatMultMethod(self.bv, val) )

    #

    def getMatrix(self) -> tuple[Mat, bool] | tuple[None, bool]:
        """
        Retrieves the matrix representation of the inner product.

        Returns
        -------
        mat: Mat
            The matrix of the inner product
        indef: bool
            Whether the matrix is indefinite
        """
        cdef Mat mat = Mat()
        cdef PetscBool indef = PETSC_FALSE
        CHKERR( BVGetMatrix(self.bv, &mat.mat, &indef) )
        if mat.mat:
            CHKERR( PetscINCREF(mat.obj) )
            return (mat, toBool(indef))
        else:
            return (None, False)

    def setMatrix(self, Mat mat or None, indef: bool = False) -> None:
        """
        Sets the bilinear form to be used for inner products.

        Parameters
        ----------
        mat
            The matrix of the inner product.
        indef
            Whether the matrix is indefinite
        """
        cdef PetscMat m = <PetscMat>NULL if mat is None else mat.mat
        cdef PetscBool tval = PETSC_TRUE if indef else PETSC_FALSE
        CHKERR( BVSetMatrix(self.bv, m, tval) )

    def applyMatrix(self, Vec x, Vec y) -> None:
        """
        Multiplies a vector with the matrix associated to the bilinear
        form.

        Parameters
        ----------
        x
            The input vector.
        y
            The result vector.

        Notes
        -----
        If the bilinear form has no associated matrix this function
        copies the vector.
        """
        CHKERR( BVApplyMatrix(self.bv, x.vec, y.vec) )

    def setActiveColumns(self, l: int, k: int) -> None:
        """
        Specify the columns that will be involved in operations.

        Parameters
        ----------
        l
            The leading number of columns.
        k
            The active number of columns.
        """
        cdef PetscInt ival1 = asInt(l)
        cdef PetscInt ival2 = asInt(k)
        CHKERR( BVSetActiveColumns(self.bv, ival1, ival2) )

    def getActiveColumns(self) -> tuple[int, int]:
        """
        Returns the current active dimensions.

        Returns
        -------
        l: int
            The leading number of columns.
        k: int
            The active number of columns.
        """
        cdef PetscInt l=0, k=0
        CHKERR( BVGetActiveColumns(self.bv, &l, &k) )
        return (toInt(l), toInt(k))

    def scaleColumn(self, j: int, alpha: Scalar) -> None:
        """
        Scale column j by alpha

        Parameters
        ----------
        j
            column number to be scaled.
        alpha
            scaling factor.
        """
        cdef PetscInt ival = asInt(j)
        cdef PetscScalar sval = asScalar(alpha)
        CHKERR( BVScaleColumn(self.bv, ival, sval) )

    def scale(self, alpha: Scalar) -> None:
        """
        Multiply the entries by a scalar value.

        Parameters
        ----------
        alpha
            scaling factor.

        Notes
        -----
        All active columns (except the leading ones) are scaled.
        """
        cdef PetscScalar sval = asScalar(alpha)
        CHKERR( BVScale(self.bv, sval) )

    def insertVec(self, j: int, Vec w) -> None:
        """
        Insert a vector into the specified column.

        Parameters
        ----------
        j
            The column to be overwritten.
        w
            The vector to be copied.
        """
        cdef PetscInt ival = asInt(j)
        CHKERR( BVInsertVec(self.bv, ival, w.vec) )

    def insertVecs(self, s: int, W: Vec | list[Vec], orth: bool = False) -> int:
        """
        Insert a set of vectors into specified columns.

        Parameters
        ----------
        s
            The first column to be overwritten.
        W
            Set of vectors to be copied.
        orth
            Flag indicating if the vectors must be orthogonalized.

        Returns
        -------
        int
            Number of linearly independent vectors.

        Notes
        -----
        Copies the contents of vectors W into self(:,s:s+n), where n is the
        length of W. If orthogonalization flag is set then the vectors are
        copied one by one then orthogonalized against the previous one.  If any
        are linearly dependent then it is discared and the value of m is
        decreased.
        """
        if isinstance(W, Vec): W = [W]
        cdef PetscInt ival = asInt(s)
        cdef PetscVec *ws = NULL
        cdef Py_ssize_t i = 0, ns = len(W)
        cdef tmp = allocate(<size_t>ns*sizeof(PetscVec),<void**>&ws)
        for i in range(ns): ws[i] = (<Vec?>W[i]).vec
        cdef PetscInt m = <PetscInt>ns
        cdef PetscBool tval = PETSC_TRUE if orth else PETSC_FALSE
        CHKERR( BVInsertVecs(self.bv, ival, &m, ws, tval) )
        return toInt(m)

    def insertConstraints(self, C: Vec | list[Vec]) -> int:
        """
        Insert a set of vectors as constraints.

        Parameters
        ----------
        C
            Set of vectors to be inserted as constraints.

        Returns
        -------
        int
            Number of constraints.

        Notes
        -----
        The constraints are relevant only during orthogonalization. Constraint
        vectors span a subspace that is deflated in every orthogonalization
        operation, so they are intended for removing those directions from the
        orthogonal basis computed in regular BV columns.
        """
        if isinstance(C, Vec): C = [C]
        cdef PetscVec *cs = NULL
        cdef Py_ssize_t i = 0, nc = len(C)
        cdef tmp = allocate(<size_t>nc*sizeof(PetscVec),<void**>&cs)
        for i in range(nc): cs[i] = (<Vec?>C[i]).vec
        cdef PetscInt m = <PetscInt>nc
        CHKERR( BVInsertConstraints(self.bv, &m, cs) )
        return toInt(m)

    def setNumConstraints(self, nc: int) -> None:
        """
        Sets the number of constraints.

        Parameters
        ----------
        nc
            The number of constraints.
        """
        cdef PetscInt val = asInt(nc)
        CHKERR( BVSetNumConstraints(self.bv, val) )

    def getNumConstraints(self) -> int:
        """
        Gets the number of constraints.

        Returns
        -------
        int
            The number of constraints.
        """
        cdef PetscInt val = 0
        CHKERR( BVGetNumConstraints(self.bv, &val) )
        return toInt(val)

    def createVec(self) -> Vec:
        """
        Creates a new Vec object with the same type and dimensions as
        the columns of the basis vectors object.

        Returns
        -------
        Vec
            New vector.
        """
        cdef Vec v = Vec()
        CHKERR( BVCreateVec(self.bv, &v.vec) )
        return v

    def setVecType(self, vec_type: Vec.Type | str) -> None:
        """
        Set the vector type.

        Parameters
        ----------
        vec_type
            Vector type used when creating vectors with `createVec`.
        """
        cdef PetscVecType cval = NULL
        vec_type = str2bytes(vec_type, &cval)
        CHKERR( BVSetVecType(self.bv, cval) )

    def getVecType(self) -> str:
        """
        Return the vector type used by the basis vectors object.
        """
        cdef PetscVecType cval = NULL
        CHKERR( BVGetVecType(self.bv, &cval) )
        return bytes2str(cval)

    def copyVec(self, j: int, Vec v) -> None:
        """
        Copies one of the columns of a basis vectors object into a Vec.

        Parameters
        ----------
        j
            The column number to be copied.
        v
            A vector.
        """
        cdef PetscInt ival = asInt(j)
        CHKERR( BVCopyVec(self.bv, ival, v.vec) )

    def copyColumn(self, j: int, i: int) -> None:
        """
        Copies the values from one of the columns to another one.

        Parameters
        ----------
        j
            The number of the source column.
        i
            The number of the destination column.
        """
        cdef PetscInt ival1 = asInt(j)
        cdef PetscInt ival2 = asInt(i)
        CHKERR( BVCopyColumn(self.bv, ival1, ival2) )

    def setDefiniteTolerance(self, deftol: float) -> None:
        """
        Sets the tolerance to be used when checking a definite inner product.

        Parameters
        ----------
        deftol
            The tolerance.
        """
        cdef PetscReal val = asReal(deftol)
        CHKERR( BVSetDefiniteTolerance(self.bv, val) )

    def getDefiniteTolerance(self) -> float:
        """
        Gets the tolerance to be used when checking a definite inner product.

        Returns
        -------
        float
            The tolerance.
        """
        cdef PetscReal val = 0
        CHKERR( BVGetDefiniteTolerance(self.bv, &val) )
        return toReal(val)

    def dotVec(self, Vec v) -> ArrayScalar:
        """
        Computes multiple dot products of a vector against all the column
        vectors of a BV.

        Parameters
        ----------
        v
            A vector.

        Returns
        -------
        ArrayScalar 
            The computed values.

        Notes
        -----
        This is analogue to VecMDot(), but using BV to represent a collection
        of vectors. The result is m = X^H*y, so m_i is equal to x_j^H y. Note
        that here X is transposed as opposed to BVDot().

        If a non-standard inner product has been specified with BVSetMatrix(),
        then the result is m = X^H*B*y.
        """
        l, k = self.getActiveColumns()
        cdef PetscScalar* mval = NULL
        cdef tmp = allocate(<size_t>(k - l)*sizeof(PetscScalar), <void**>&mval)

        CHKERR( BVDotVec(self.bv, v.vec, mval) )

        cdef object m = None
        m = array_s(k - l, mval)
        return m

    def dotColumn(self, j: int) -> ArrayScalar:
        """
        Computes multiple dot products of a column against all the column
        vectors of a BV.

        Parameters
        ----------
        j
            The index of the column.

        Returns
        -------
        ArrayScalar
            The computed values.
        """
        cdef PetscInt ival = asInt(j)
        l, k = self.getActiveColumns()
        cdef PetscScalar* mval = NULL
        cdef tmp = allocate(<size_t>(k - l)*sizeof(PetscScalar), <void**>&mval)

        CHKERR( BVDotColumn(self.bv, ival, mval) )

        cdef object m = None
        m = array_s(k - l, mval)
        return m

    def getColumn(self, j: int) -> Vec:
        """
        Returns a Vec object that contains the entries of the requested column
        of the basis vectors object.

        Parameters
        ----------
        j
            The index of the requested column.

        Returns
        -------
        Vec
            The vector containing the jth column.

        Notes
        -----
        Modifying the returned Vec will change the BV entries as well.
        """
        cdef Vec v = Vec()
        cdef PetscInt ival = asInt(j)
        CHKERR( BVGetColumn(self.bv, j, &v.vec) )
        CHKERR( PetscINCREF(v.obj) )
        return v

    def restoreColumn(self, j: int, Vec v) -> None:
        """
        Restore a column obtained with `getColumn()`.

        Parameters
        ----------
        j
            The index of the requested column.
        v
            The vector obtained with `getColumn()`.

        Notes
        -----
        The arguments must match the corresponding call to `getColumn()`.
        """
        cdef PetscInt ival = asInt(j)
        CHKERR( PetscObjectDereference(<PetscObject>v.vec) )
        CHKERR( BVRestoreColumn(self.bv, ival, &v.vec) )

    def getMat(self) -> Mat:
        """
        Returns a Mat object of dense type that shares the memory
        of the basis vectors object.

        Returns
        -------
        Mat
            The matrix.

        Notes
        -----
        The returned matrix contains only the active columns. If the content
        of the Mat is modified, these changes are also done in the BV object.
        The user must call `restoreMat()` when no longer needed.
        """
        cdef Mat A = Mat()
        CHKERR( BVGetMat(self.bv, &A.mat) )
        CHKERR( PetscINCREF(A.obj) )
        return A

    def restoreMat(self, Mat A) -> None:
        """
        Restores the Mat obtained with `getMat()`.

        Parameters
        ----------
        A
            The matrix obtained with `getMat()`.

        Notes
        -----
        A call to this function must match a previous call of `getMat()`.
        The effect is that the contents of the Mat are copied back to the
        BV internal data structures.
        """
        CHKERR( PetscObjectDereference(<PetscObject>A.mat) )
        CHKERR( BVRestoreMat(self.bv, &A.mat) )

    def dot(self, BV Y) -> Mat:
        """
        Computes the 'block-dot' product of two basis vectors objects.
            M = Y^H*X (m_ij = y_i^H x_j) or M = Y^H*B*X

        Parameters
        ----------
        Y
            Left basis vectors, can be the same as self, giving M = X^H X.

        Returns
        -------
        Mat
            The resulting matrix.

        Notes
        -----
        This is the generalization of VecDot() for a collection of vectors, M =
        Y^H*X. The result is a matrix M whose entry m_ij is equal to y_i^H x_j
        (where y_i^H denotes the conjugate transpose of y_i).

        X and Y can be the same object.

        If a non-standard inner product has been specified with setMatrix(),
        then the result is M = Y^H*B*X. In this case, both X and Y must have
        the same associated matrix.

        Only rows (resp. columns) of M starting from ly (resp. lx) are
        computed, where ly (resp. lx) is the number of leading columns of Y
        (resp. X).
        """
        cdef BV X = self
        cdef PetscInt ky=0, kx=0
        CHKERR( BVGetActiveColumns(Y.bv, NULL, &ky) )
        CHKERR( BVGetActiveColumns(X.bv, NULL, &kx) )
        cdef Mat M = Mat().createDense((ky, kx), comm=COMM_SELF).setUp()
        CHKERR( BVDot(X.bv, Y.bv, M.mat) )
        return M

    def matProject(self, Mat A or None, BV Y) -> Mat:
        """
        Computes the projection of a matrix onto a subspace.

        M = Y^H A X

        Parameters
        ----------
        A
            Matrix to be projected.
        Y
            Left basis vectors, can be the same as self, giving M = X^H A X.

        Returns
        -------
        Mat
            Projection of the matrix A onto the subspace.
        """
        cdef BV X = self
        cdef PetscInt  kx=0, ky=0
        CHKERR( BVGetActiveColumns(X.bv, NULL, &kx) )
        CHKERR( BVGetActiveColumns(Y.bv, NULL, &ky) )
        cdef PetscMat Amat = <PetscMat>NULL if A is None else A.mat
        cdef Mat M = Mat().createDense((ky, kx), comm=COMM_SELF).setUp()
        CHKERR( BVMatProject(X.bv, Amat, Y.bv, M.mat) )
        return M

    def matMult(self, Mat A, BV Y=None) -> BV:
        """
        Computes the matrix-vector product for each column, Y = A*V.

        Parameters
        ----------
        A
            The matrix.

        Returns
        -------
        BV
            The result.

        Notes
        -----
        Only active columns (excluding the leading ones) are processed.

        It is possible to choose whether the computation is done column by column
        or using dense matrices using the options database keys:

            -bv_matmult_vecs
            -bv_matmult_mat

        The default is bv_matmult_mat.
        """
        cdef MPI_Comm comm = PetscObjectComm(<PetscObject>self.bv)
        cdef SlepcBVType bv_type = NULL
        cdef PetscInt n=0, N=0, m=0
        cdef SlepcBVOrthogType val1 = BV_ORTHOG_CGS
        cdef SlepcBVOrthogRefineType val2 = BV_ORTHOG_REFINE_IFNEEDED
        cdef SlepcBVOrthogBlockType val3 = BV_ORTHOG_BLOCK_GS
        cdef PetscReal rval = PETSC_DEFAULT
        if Y is None: Y = BV()
        if Y.bv == NULL:
            CHKERR( BVGetType(self.bv, &bv_type) )
            CHKERR( MatGetLocalSize(A.mat, &n, NULL) )
            CHKERR( MatGetSize(A.mat, &N, NULL) )
            CHKERR( BVGetSizes(self.bv, NULL, NULL, &m) )
            CHKERR( BVGetOrthogonalization(self.bv, &val1, &val2, &rval, &val3) )
            CHKERR( BVCreate(comm, &Y.bv) )
            CHKERR( BVSetType(Y.bv, bv_type) )
            CHKERR( BVSetSizes(Y.bv, n, N, m) )
            CHKERR( BVSetOrthogonalization(Y.bv, val1, val2, rval, val3) )
        CHKERR( BVMatMult(self.bv, A.mat, Y.bv) )
        return Y

    def matMultHermitianTranspose(self, Mat A, BV Y=None) -> BV:
        """
        Computes the matrix-vector product with the conjugate transpose of a
        matrix for each column, Y=A^H*V.

        Parameters
        ----------
        A
            The matrix.

        Returns
        -------
        BV
            The result.

        Notes
        -----
        Only active columns (excluding the leading ones) are processed.

        As opoosed to matMult(), this operation is always done by column by
        column, with a sequence of calls to MatMultHermitianTranspose().
        """
        cdef MPI_Comm comm = PetscObjectComm(<PetscObject>self.bv)
        cdef SlepcBVType bv_type = NULL
        cdef PetscInt n=0, N=0, m=0
        cdef SlepcBVOrthogType val1 = BV_ORTHOG_CGS
        cdef SlepcBVOrthogRefineType val2 = BV_ORTHOG_REFINE_IFNEEDED
        cdef SlepcBVOrthogBlockType val3 = BV_ORTHOG_BLOCK_GS
        cdef PetscReal rval = PETSC_DEFAULT
        if Y is None: Y = BV()
        if Y.bv == NULL:
            CHKERR( BVGetType(self.bv, &bv_type) )
            CHKERR( MatGetLocalSize(A.mat, NULL, &n) )
            CHKERR( MatGetSize(A.mat, NULL, &N) )
            CHKERR( BVGetSizes(self.bv, NULL, NULL, &m) )
            CHKERR( BVGetOrthogonalization(self.bv, &val1, &val2, &rval, &val3) )
            CHKERR( BVCreate(comm, &Y.bv) )
            CHKERR( BVSetType(Y.bv, bv_type) )
            CHKERR( BVSetSizes(Y.bv, n, N, m) )
            CHKERR( BVSetOrthogonalization(Y.bv, val1, val2, rval, val3) )
        CHKERR( BVMatMultHermitianTranspose(self.bv, A.mat, Y.bv) )
        return Y

    def matMultColumn(self, Mat A, j: int) -> None:
        """
        Computes the matrix-vector product for a specified column, storing
        the result in the next column: v_{j+1}=A*v_j.

        Parameters
        ----------
        A
            The matrix.
        j
            Index of column.
        """
        cdef PetscInt ival = asInt(j)
        CHKERR( BVMatMultColumn(self.bv, A.mat, ival) )

    def matMultTransposeColumn(self, Mat A, j: int) -> None:
        """
        Computes the transpose matrix-vector product for a specified column,
        storing the result in the next column: v_{j+1}=A^T*v_j.

        Parameters
        ----------
        A
            The matrix.
        j
            Index of column.
        """
        cdef PetscInt ival = asInt(j)
        CHKERR( BVMatMultTransposeColumn(self.bv, A.mat, ival) )

    def matMultHermitianTransposeColumn(self, Mat A, j: int) -> None:
        """
        Computes the conjugate-transpose matrix-vector product for a specified column,
        storing the result in the next column: v_{j+1}=A^H*v_j.

        Parameters
        ----------
        A
            The matrix.
        j
            Index of column.
        """
        cdef PetscInt ival = asInt(j)
        CHKERR( BVMatMultHermitianTransposeColumn(self.bv, A.mat, ival) )

    def mult(self, alpha: Scalar, beta: Scalar, BV X, Mat Q or None) -> None:
        """
        Computes Y = beta*Y + alpha*X*Q.

        Parameters
        ----------
        alpha
            Coefficient that multiplies X.
        beta
            Coefficient that multiplies Y.
        X
            Input basis vectors.
        Q
            Input matrix, if not given the identity matrix is assumed.
        """
        cdef PetscScalar sval1 = asScalar(alpha)
        cdef PetscScalar sval2 = asScalar(beta)
        cdef PetscMat Qmat = <PetscMat>NULL if Q is None else Q.mat
        CHKERR( BVMult(self.bv, sval1, sval2, X.bv, Qmat) )

    def multInPlace(self, Mat Q, s: int, e: int) -> None:
        """
        Update a set of vectors as V(:,s:e-1) = V*Q(:,s:e-1).

        Parameters
        ----------
        Q
            A sequential dense matrix.
        s
            First column to be overwritten.
        e
            Last column to be overwritten.
        """
        cdef PetscInt ival1 = asInt(s)
        cdef PetscInt ival2 = asInt(e)
        CHKERR( BVMultInPlace(self.bv, Q.mat, ival1, ival2) )

    def multColumn(self, alpha: Scalar, beta: Scalar, j: int, q: Sequence[Scalar]) -> None:
        """
        Computes y = beta*y + alpha*X*q, where y is the j-th column.

        Parameters
        ----------
        alpha
            Coefficient that multiplies X.
        beta
            Coefficient that multiplies y.
        j
            The column index.
        q
            Input coefficients.
        """
        cdef PetscScalar sval1 = asScalar(alpha)
        cdef PetscScalar sval2 = asScalar(beta)
        cdef PetscInt ival = asInt(j)
        cdef PetscInt nq = 0
        cdef PetscScalar* qval = NULL
        cdef tmp = iarray_s(q, &nq, &qval)
        cdef PetscInt l=0, k=0
        CHKERR( BVGetActiveColumns(self.bv, &l, &k) )
        assert nq == k-l
        CHKERR( BVMultColumn(self.bv, sval1, sval2, ival, qval) )

    def multVec(self, alpha: Scalar, beta: Scalar, Vec y, q: Sequence[Scalar]) -> None:
        """
        Computes y = beta*y + alpha*X*q.

        Parameters
        ----------
        alpha
            Coefficient that multiplies X.
        beta
            Coefficient that multiplies y.
        y
            Input/output vector.
        q
            Input coefficients.
        """
        cdef PetscScalar sval1 = asScalar(alpha)
        cdef PetscScalar sval2 = asScalar(beta)
        cdef PetscInt nq = 0
        cdef PetscScalar* qval = NULL
        cdef tmp = iarray_s(q, &nq, &qval)
        cdef PetscInt l=0, k=0
        CHKERR( BVGetActiveColumns(self.bv, &l, &k) )
        assert nq == k-l
        CHKERR( BVMultVec(self.bv, sval1, sval2, y.vec, qval) )

    def normColumn(self, j: int, norm_type: NormType | None = None) -> float:
        """
        Computes the matrix norm of the BV.

        Parameters
        ----------
        j
            Index of column.
        norm_type
            The norm type.

        Returns
        -------
        float
            The norm.

        Notes
        -----
        The norm of V[j] is computed (NORM_1, NORM_2, or NORM_INFINITY).

        If a non-standard inner product has been specified with BVSetMatrix(),
        then the returned value is ``sqrt(V[j]'* B*V[j])``, where B is the inner
        product matrix (argument 'type' is ignored).
        """
        cdef PetscNormType ntype = PETSC_NORM_2
        if norm_type is not None: ntype = norm_type
        cdef PetscReal norm = 0
        CHKERR( BVNormColumn(self.bv, j, ntype, &norm) )
        return toReal(norm)

    def norm(self, norm_type: NormType | None = None) -> float:
        """
        Computes the matrix norm of the BV.

        Parameters
        ----------
        norm_type
            The norm type.

        Returns
        -------
        float
            The norm.

        Notes
        -----
        All active columns (except the leading ones) are considered as a
        matrix. The allowed norms are NORM_1, NORM_FROBENIUS, and
        NORM_INFINITY.

        This operation fails if a non-standard inner product has been specified
        with BVSetMatrix().
        """
        cdef PetscNormType ntype = PETSC_NORM_FROBENIUS
        if norm_type is not None: ntype = norm_type
        cdef PetscReal norm = 0
        CHKERR( BVNorm(self.bv, ntype, &norm) )
        return toReal(norm)

    def resize(self, m: int, copy: bool = True) -> None:
        """
        Change the number of columns.

        Parameters
        ----------
        m
            The new number of columns.
        copy
            A flag indicating whether current values should be kept.

        Notes
        -----
        Internal storage is reallocated. If copy is True, then the contents are
        copied to the leading part of the new space.
        """
        cdef PetscInt ival = asInt(m)
        cdef PetscBool tval = PETSC_TRUE if copy else PETSC_FALSE
        CHKERR( BVResize(self.bv, ival, tval) )

    def setRandom(self) -> None:
        """
        Set the active columns of the BV to random numbers.

        Notes
        -----
        All active columns (except the leading ones) are modified.
        """
        CHKERR( BVSetRandom(self.bv) )

    def setRandomNormal(self) -> None:
        """
        Set the active columns of the BV to random numbers (with normal
        distribution).

        Notes
        -----
        All active columns (except the leading ones) are modified.
        """
        CHKERR( BVSetRandomNormal(self.bv) )

    def setRandomSign(self) -> None:
        """
        Set the entries of a BV to values 1 or -1 with equal probability.

        Notes
        -----
        All active columns (except the leading ones) are modified.
        """
        CHKERR( BVSetRandomSign(self.bv) )

    def setRandomColumn(self, j: int) -> None:
        """
        Set one column of the BV to random numbers.

        Parameters
        ----------
        j
            Column number to be set.
        """
        cdef PetscInt ival = asInt(j)
        CHKERR( BVSetRandomColumn(self.bv, ival) )

    def setRandomCond(self, condn: float) -> None:
        """
        Set the columns of a BV to random numbers, in a way that the generated
        matrix has a given condition number.

        Parameters
        ----------
        condn
            Condition number.
        """
        cdef PetscReal rval = asReal(condn)
        CHKERR( BVSetRandomCond(self.bv, rval) )

    def setRandomContext(self, Random rnd) -> None:
        """
        Sets the `PETSc.Random` object associated with the BV, to be used
        in operations that need random numbers.

        Parameters
        ----------
        rnd
            The random number generator context.
        """
        CHKERR( BVSetRandomContext(self.bv, rnd.rnd) )

    def getRandomContext(self) -> Random:
        """
        Gets the `PETSc.Random` object associated with the BV.

        Returns
        -------
        PETSc.Random
            The random number generator context.
        """
        cdef Random rnd = Random()
        CHKERR( BVGetRandomContext(self.bv, &rnd.rnd) )
        CHKERR( PetscINCREF(rnd.obj) )
        return rnd

    def orthogonalizeVec(self, Vec v) -> tuple[float, bool]:
        """
        Orthogonalize a vector with respect to a set of vectors.

        Parameters
        ----------
        v
            Vector to be orthogonalized, modified on return.

        Returns
        -------
        norm: float
            The norm of the resulting vector.
        lindep: bool
            Flag indicating that refinement did not improve the
            quality of orthogonalization.

        Notes
        -----
        This function applies an orthogonal projector to project
        vector ``v`` onto the orthogonal complement of the span of the
        columns of the BV.

        This routine does not normalize the resulting vector.
        """
        cdef PetscReal norm = 0
        cdef PetscBool ldep = PETSC_FALSE
        CHKERR( BVOrthogonalizeVec(self.bv, v.vec, NULL, &norm, &ldep) )
        return (toReal(norm), toBool(ldep))

    def orthogonalizeColumn(self, j: int) -> tuple[float, bool]:
        """
        Orthogonalize one of the column vectors with respect to the previous ones.

        Parameters
        ----------
        j
            Index of the column to be orthogonalized.

        Returns
        -------
        norm: float
            The norm of the resulting vector.
        lindep: bool
            Flag indicating that refinement did not improve the
            quality of orthogonalization.

        Notes
        -----
        This function applies an orthogonal projector to project
        vector ``V[j]`` onto the orthogonal complement of the span of the
        columns ``V[0..j-1]``, where ``V[.]`` are the vectors of the BV.
        The columns ``V[0..j-1]`` are assumed to be mutually orthonormal.

        This routine does not normalize the resulting vector.
        """
        cdef PetscInt ival = asInt(j)
        cdef PetscReal norm = 0
        cdef PetscBool ldep = PETSC_FALSE
        CHKERR( BVOrthogonalizeColumn(self.bv, ival, NULL, &norm, &ldep) )
        return (toReal(norm), toBool(ldep))

    def orthonormalizeColumn(self, j: int, replace: bool = False) -> tuple[float, bool]:
        """
        Orthonormalize one of the column vectors with respect to the previous
        ones.  This is equivalent to a call to `orthogonalizeColumn()`
        followed by a call to `scaleColumn()` with the reciprocal of the norm.

        Parameters
        ----------
        j
            Index of the column to be orthonormalized.
        replace
            Whether it is allowed to set the vector randomly.

        Returns
        -------
        norm: float
            The norm of the resulting vector.
        lindep: bool
            Flag indicating that refinement did not improve the
            quality of orthogonalization.
        """
        cdef PetscInt ival = asInt(j)
        cdef PetscBool bval = PETSC_FALSE
        if replace is not None: bval = asBool(replace)
        cdef PetscReal norm = 0
        cdef PetscBool ldep = PETSC_FALSE
        CHKERR( BVOrthonormalizeColumn(self.bv, ival, bval, &norm, &ldep) )
        return (toReal(norm), toBool(ldep))

    def orthogonalize(self, Mat R=None, **kargs: Any) -> None:
        """
        Orthogonalize all columns (except leading ones),
        that is, compute the QR decomposition.

        Parameters
        ----------
        R
            A sequential dense matrix.

        Notes
        -----
        The output satisfies ``V0 = V*R`` (where V0 represent the input V) and ``V'*V = I``.
        """
        if kargs: self.setOrthogonalization(**kargs)
        cdef PetscMat Rmat = <PetscMat>NULL if R is None else R.mat
        CHKERR( BVOrthogonalize(self.bv, Rmat) )

    #

    property sizes:
        """Basis vectors local and global sizes, and the number of columns."""
        def __get__(self) -> tuple[LayoutSizeSpec, int]:
            return self.getSizes()

    property size:
        """Basis vectors global size."""
        def __get__(self) -> tuple[int, int]:
            return self.getSizes()[0][0]

    property local_size:
        """Basis vectors local size."""
        def __get__(self) -> int:
            return self.getSizes()[0][1]

    property column_size:
        """Basis vectors column size."""
        def __get__(self) -> int:
            return self.getSizes()[1]

# -----------------------------------------------------------------------------

del BVType
del BVOrthogType
del BVOrthogRefineType
del BVOrthogBlockType
del BVMatMultType

# -----------------------------------------------------------------------------
