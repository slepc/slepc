# -----------------------------------------------------------------------------

class BVType(object):
    """
    BV type.

    - `MAT`: A `BV` stored as a dense `petsc.Mat`.
    - `SVEC`: A `BV` stored as a single `petsc.Vec`.
    - `VECS`: A `BV` stored as an array of independent `petsc.Vec`.
    - `CONTIGUOUS`: A `BV` stored as an array of `petsc.Vec`
      sharing a contiguous array of scalars.
    - `TENSOR`: A special `BV` represented in compact form as
      :math:`V = (I \otimes U) S`.

    See Also
    --------
    slepc.BVType
    """
    MAT        = S_(BVMAT)
    SVEC       = S_(BVSVEC)
    VECS       = S_(BVVECS)
    CONTIGUOUS = S_(BVCONTIGUOUS)
    TENSOR     = S_(BVTENSOR)

class BVOrthogType(object):
    """
    BV orthogonalization types.

    - `CGS`: Classical Gram-Schmidt.
    - `MGS`: Modified Gram-Schmidt.

    See Also
    --------
    slepc.BVOrthogType
    """
    CGS = BV_ORTHOG_CGS
    MGS = BV_ORTHOG_MGS

class BVOrthogRefineType(object):
    """
    BV orthogonalization refinement types.

    - `IFNEEDED`: Reorthogonalize if a criterion is satisfied.
    - `NEVER`:    Never reorthogonalize.
    - `ALWAYS`:   Always reorthogonalize.

    See Also
    --------
    slepc.BVOrthogRefineType
    """
    IFNEEDED = BV_ORTHOG_REFINE_IFNEEDED
    NEVER    = BV_ORTHOG_REFINE_NEVER
    ALWAYS   = BV_ORTHOG_REFINE_ALWAYS

class BVOrthogBlockType(object):
    """
    BV block-orthogonalization types.

    - `GS`:       Gram-Schmidt, column by column.
    - `CHOL`:     Cholesky QR method.
    - `TSQR`:     Tall-skinny QR method.
    - `TSQRCHOL`: Tall-skinny QR, but computing the triangular factor only.
    - `SVQB`:     SVQB method.

    See Also
    --------
    slepc.BVOrthogBlockType
    """
    GS       = BV_ORTHOG_BLOCK_GS
    CHOL     = BV_ORTHOG_BLOCK_CHOL
    TSQR     = BV_ORTHOG_BLOCK_TSQR
    TSQRCHOL = BV_ORTHOG_BLOCK_TSQRCHOL
    SVQB     = BV_ORTHOG_BLOCK_SVQB

class BVMatMultType(object):
    """
    BV mat-mult types.

    - `VECS`: Perform a matrix-vector multiply per each column.
    - `MAT`:  Carry out a Mat-Mat product with a dense matrix.

    See Also
    --------
    slepc.BVMatMultType
    """
    VECS     = BV_MATMULT_VECS
    MAT      = BV_MATMULT_MAT

class BVSVDMethod(object):
    """
    BV methods for computing the SVD.

    - `REFINE`: Based on the SVD of the cross product matrix :math:`S^* S`,
      with refinement.
    - `QR`:     Based on the SVD of the triangular factor of qr(S).
    - `QR_CAA`: Variant of QR intended for use in communication-avoiding.
      Arnoldi.

    See Also
    --------
    slepc.BVSVDMethod
    """
    REFINE   = BV_SVD_METHOD_REFINE
    QR       = BV_SVD_METHOD_QR
    QR_CAA   = BV_SVD_METHOD_QR_CAA

# -----------------------------------------------------------------------------

cdef class BV(Object):

    """
    Basis Vectors.

    The `BV` package provides the concept of a block of vectors that
    represent the basis of a subspace. It is a convenient way of handling
    a collection of vectors that often operate together, rather than
    working with an array of `petsc4py.PETSc.Vec`.
    """

    Type             = BVType
    OrthogType       = BVOrthogType
    OrthogRefineType = BVOrthogRefineType
    RefineType       = BVOrthogRefineType
    OrthogBlockType  = BVOrthogBlockType
    BlockType        = BVOrthogBlockType
    MatMultType      = BVMatMultType
    SVDMethod        = BVSVDMethod

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
        Print the BV data structure.

        Collective.

        Parameters
        ----------
        viewer
            Visualization context; if not provided, the standard
            output is used.

        See Also
        --------
        slepc.BVView
        """
        cdef PetscViewer vwr = def_Viewer(viewer)
        CHKERR( BVView(self.bv, vwr) )

    def destroy(self) -> Self:
        """
        Destroy the BV object.

        Collective.

        See Also
        --------
        slepc.BVDestroy
        """
        CHKERR( BVDestroy(&self.bv) )
        self.bv = NULL
        return self

    def create(self, comm: Comm | None = None) -> Self:
        """
        Create the BV object.

        Collective.

        Parameters
        ----------
        comm
            MPI communicator; if not provided, it defaults to all
            processes.

        See Also
        --------
        createFromMat, slepc.BVCreate
        """
        cdef MPI_Comm ccomm = def_Comm(comm, SLEPC_COMM_DEFAULT())
        cdef SlepcBV newbv = NULL
        CHKERR( BVCreate(ccomm, &newbv) )
        CHKERR( SlepcCLEAR(self.obj) ); self.bv = newbv
        return self

    def createFromMat(self, Mat A) -> Self:
        """
        Create a basis vectors object from a dense matrix.

        Collective.

        Parameters
        ----------
        A
            A dense tall-skinny matrix.

        Notes
        -----
        The matrix values are copied to the `BV` data storage, memory is not
        shared.

        The communicator of the `BV` object will be the same as `A`, and so
        will be the dimensions.

        See Also
        --------
        create, createMat, slepc.BVCreateFromMat
        """
        cdef SlepcBV newbv = NULL
        CHKERR( BVCreateFromMat(A.mat, &newbv) )
        CHKERR( SlepcCLEAR(self.obj) ); self.bv = newbv
        return self

    def createMat(self) -> Mat:
        """
        Create a new dense matrix and copy the contents of the BV.

        Collective.

        Returns
        -------
        petsc4py.PETSc.Mat
            The new matrix.

        Notes
        -----
        The matrix contains all columns of the `BV`, not just the active
        columns.

        See Also
        --------
        createFromMat, createVec, getMat, slepc.BVCreateMat
        """
        cdef Mat mat = Mat()
        CHKERR( BVCreateMat(self.bv, &mat.mat) )
        return mat

    def duplicate(self) -> BV:
        """
        Duplicate the BV object with the same type and dimensions.

        Collective.

        Returns
        -------
        BV
            The new object.

        Notes
        -----
        This function does not copy the entries, it just allocates the
        storage for the new `BV`. Use `copy()` to copy the content.

        See Also
        --------
        duplicateResize, slepc.BVDuplicate
        """
        cdef BV bv = type(self)()
        CHKERR( BVDuplicate(self.bv, &bv.bv) )
        return bv

    def duplicateResize(self, m: int) -> BV:
        """
        Create a BV object of the same type and dimensions as an existing one.

        Collective.

        Parameters
        ----------
        m
            The number of columns.

        Returns
        -------
        BV
            The new object.

        Notes
        -----
        This is equivalent to a call to `duplicate()` followed by `resize()`
        with possibly different number of columns.
        The contents of this `BV` are not copied to the new one.

        See Also
        --------
        duplicate, resize, slepc.BVDuplicateResize
        """
        cdef BV bv = type(self)()
        cdef PetscInt ival = asInt(m)
        CHKERR( BVDuplicateResize(self.bv, ival, &bv.bv) )
        return bv

    def copy(self, BV result=None) -> BV:
        """
        Copy a basis vector object into another one.

        Logically collective.

        Returns
        -------
        BV
            The copy.

        Parameters
        ----------
        result
            The copy.

        Notes
        -----
        Both objects must be distributed in the same manner; local copies are
        done. Only active columns (excluding the leading ones) are copied.
        In the destination BV, columns are overwritten starting from the
        leading ones. Constraints are not copied.

        See Also
        --------
        slepc.BVCopy
        """
        if result is None:
            result = type(self)()
        if result.bv == NULL:
            CHKERR( BVDuplicate(self.bv, &result.bv) )
        CHKERR( BVCopy(self.bv, result.bv) )
        return result

    def setType(self, bv_type: Type | str) -> None:
        """
        Set the type for the BV object.

        Logically collective.

        Parameters
        ----------
        bv_type
            The basis vectors type to be used.

        See Also
        --------
        getType, slepc.BVSetType
        """
        cdef SlepcBVType cval = NULL
        bv_type = str2bytes(bv_type, &cval)
        CHKERR( BVSetType(self.bv, cval) )

    def getType(self) -> str:
        """
        Get the BV type of this object.

        Not collective.

        Returns
        -------
        str
            The basis vectors type currently being used.

        See Also
        --------
        setType, slepc.BVGetType
        """
        cdef SlepcBVType bv_type = NULL
        CHKERR( BVGetType(self.bv, &bv_type) )
        return bytes2str(bv_type)

    def setSizes(self, sizes: LayoutSizeSpec, m: int) -> None:
        """
        Set the local and global sizes, and the number of columns.

        Collective.

        Parameters
        ----------
        sizes
            The global size ``N`` or a two-tuple ``(n, N)``
            with the local and global sizes.
        m
            The number of columns.

        Notes
        -----
        Either ``n`` or ``N`` (but not both) can be `DETERMINE`
        or ``None`` to have it automatically set.

        See Also
        --------
        setSizesFromVec, getSizes, slepc.BVSetSizes
        """
        cdef PetscInt n=0, N=0
        cdef PetscInt ival = asInt(m)
        BV_Sizes(sizes, &n, &N)
        CHKERR( BVSetSizes(self.bv, n, N, ival) )

    def setSizesFromVec(self, Vec w, m: int) -> None:
        """
        Set the local and global sizes, and the number of columns.

        Collective.

        Local and global sizes are specified indirectly by passing a template
        vector.

        Parameters
        ----------
        w
            The template vector.
        m
            The number of columns.

        See Also
        --------
        setSizes, getSizes, slepc.BVSetSizesFromVec
        """
        cdef PetscInt ival = asInt(m)
        CHKERR( BVSetSizesFromVec(self.bv, w.vec, ival) )

    def getSizes(self) -> tuple[LayoutSizeSpec, int]:
        """
        Get the local and global sizes, and the number of columns.

        Not collective.

        Returns
        -------
        (n, N): tuple of int
            The local and global sizes.
        m: int
            The number of columns.

        See Also
        --------
        setSizes, setSizesFromVec, slepc.BVGetSizes
        """
        cdef PetscInt n=0, N=0, m=0
        CHKERR( BVGetSizes(self.bv, &n, &N, &m) )
        return ((toInt(n), toInt(N)), toInt(m))

    def setLeadingDimension(self, ld: int) -> None:
        """
        Set the leading dimension.

        Not collective.

        Parameters
        ----------
        ld
            The leading dimension.

        Notes
        -----
        This parameter is relevant for a BV of `BV.Type.MAT`.

        See Also
        --------
        getLeadingDimension, slepc.BVSetLeadingDimension
        """
        cdef PetscInt val = asInt(ld)
        CHKERR( BVSetLeadingDimension(self.bv, val) )

    def getLeadingDimension(self) -> int:
        """
        Get the leading dimension.

        Not collective.

        Returns
        -------
        int
            The leading dimension.

        Notes
        -----
        The returned value may be different in different processes.

        The leading dimension must be used when accessing the internal
        array via `getArray()`.

        See Also
        --------
        setLeadingDimension, slepc.BVGetLeadingDimension
        """
        cdef PetscInt val = 0
        CHKERR( BVGetLeadingDimension(self.bv, &val) )
        return toInt(val)

    def getArray(self, readonly: bool = False) -> ArrayScalar:
        """
        Return the array where the data is stored.

        Not collective.

        Parameters
        ----------
        readonly
            Enable to obtain a read only array.

        Returns
        -------
        ArrayScalar
            The array.

        See Also
        --------
        slepc.BVGetArray, slepc.BVGetArrayRead
        """
        cdef PetscInt m=0, N=0, lda=0, k=0, l=0
        cdef PetscScalar *data = NULL
        CHKERR(BVGetSizes(self.bv, NULL, &N, NULL))
        CHKERR(BVGetLeadingDimension(self.bv, &lda))
        CHKERR(BVGetActiveColumns(self.bv, &l, &k))
        m = k-l
        if readonly:
            CHKERR(BVGetArrayRead(self.bv, <const PetscScalar**>&data))
        else:
            CHKERR(BVGetArray(self.bv, &data))
        cdef int typenum = NPY_PETSC_SCALAR
        cdef int itemsize = <int>sizeof(PetscScalar)
        cdef int flags = NPY_ARRAY_FARRAY_RO if readonly else NPY_ARRAY_FARRAY
        cdef npy_intp dims[2], strides[2]
        dims[0] = <npy_intp>N; strides[0] = <npy_intp>sizeof(PetscScalar)
        dims[1] = <npy_intp>m; strides[1] = <npy_intp>(lda*sizeof(PetscScalar))
        cdef ndarray array = PyArray_New(<PyTypeObject*>ndarray, 2,
                                         dims, typenum, strides,
                                         data, itemsize, flags, NULL)
        Py_INCREF(self)
        PyArray_SetBaseObject(array, self)
        if readonly:
            CHKERR(BVRestoreArrayRead(self.bv, <const PetscScalar**>&data))
        else:
            CHKERR(BVRestoreArray(self.bv, &data))
        return array

    def setOptionsPrefix(self, prefix: str | None = None) -> None:
        """
        Set the prefix used for searching for all BV options in the database.

        Logically collective.

        Parameters
        ----------
        prefix
            The prefix string to prepend to all BV option requests.

        Notes
        -----
        A hyphen (``-``) must NOT be given at the beginning of the
        prefix name.  The first character of all runtime options is
        AUTOMATICALLY the hyphen.

        See Also
        --------
        appendOptionsPrefix, getOptionsPrefix, slepc.BVGetOptionsPrefix
        """
        cdef const char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( BVSetOptionsPrefix(self.bv, cval) )

    def appendOptionsPrefix(self, prefix: str | None = None) -> None:
        """
        Append to the prefix used for searching for all BV options in the database.

        Logically collective.

        Parameters
        ----------
        prefix
            The prefix string to prepend to all BV option requests.

        See Also
        --------
        setOptionsPrefix, getOptionsPrefix, slepc.BVAppendOptionsPrefix
        """
        cdef const char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( BVAppendOptionsPrefix(self.bv, cval) )

    def getOptionsPrefix(self) -> str:
        """
        Get the prefix used for searching for all BV options in the database.

        Not collective.

        Returns
        -------
        str
            The prefix string set for this BV object.

        See Also
        --------
        setOptionsPrefix, appendOptionsPrefix, slepc.BVGetOptionsPrefix
        """
        cdef const char *prefix = NULL
        CHKERR( BVGetOptionsPrefix(self.bv, &prefix) )
        return bytes2str(prefix)

    def setFromOptions(self) -> None:
        """
        Set BV options from the options database.

        Collective.

        Notes
        -----
        To see all options, run your program with the ``-help``
        option.

        See Also
        --------
        setOptionsPrefix, slepc.BVSetFromOptions
        """
        CHKERR( BVSetFromOptions(self.bv) )

    #

    def getOrthogonalization(self) -> tuple[OrthogType, OrthogRefineType, float, OrthogBlockType]:
        """
        Get the orthogonalization settings from the BV object.

        Not collective.

        Returns
        -------
        type: OrthogType
            The type of orthogonalization technique.
        refine: OrthogRefineType
            The type of refinement.
        eta: float
            Parameter for selective refinement (used when the
            refinement type is `IFNEEDED`).
        block: OrthogBlockType
            The type of block orthogonalization.

        See Also
        --------
        setOrthogonalization, slepc.BVGetOrthogonalization
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
        Set the method used for the (block-)orthogonalization of vectors.

        Logically collective.

        Ortogonalization of vectors (classical or modified Gram-Schmidt
        with or without refinement), and for the block-orthogonalization
        (simultaneous orthogonalization of a set of vectors).

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

        The parameter ``eta`` should be a real value between ``0`` and
        ``1`` (or `DETERMINE`).  The value of ``eta`` is used only when
        the refinement type is `IFNEEDED`.

        When using several processes, `MGS` is likely to result in bad
        scalability.

        If the method set for block orthogonalization is `GS`, then the
        computation is done column by column with the vector orthogonalization.

        See Also
        --------
        getOrthogonalization, slepc.BVSetOrthogonalization
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
        Get the method used for the `matMult()` operation.

        Not collective.

        Returns
        -------
        MatMultType
            The method for the `matMult()` operation.

        See Also
        --------
        matMult, setMatMultMethod, slepc.BVGetMatMultMethod
        """
        cdef SlepcBVMatMultType val = BV_MATMULT_MAT
        CHKERR( BVGetMatMultMethod(self.bv, &val) )
        return val

    def setMatMultMethod(self, method: MatMultType) -> None:
        """
        Set the method used for the `matMult()` operation.

        Logically collective.

        Parameters
        ----------
        method
            The method for the `matMult()` operation.

        See Also
        --------
        matMult, getMatMultMethod, slepc.BVSetMatMultMethod
        """
        cdef SlepcBVMatMultType val = method
        CHKERR( BVSetMatMultMethod(self.bv, val) )

    #

    def getMatrix(self) -> tuple[Mat, bool] | tuple[None, bool]:
        """
        Get the matrix representation of the inner product.

        Not collective.

        Returns
        -------
        B: petsc4py.PETSc.Mat
            The matrix of the inner product.
        indef: bool
            Whether the matrix is indefinite.

        See Also
        --------
        setMatrix, slepc.BVGetMatrix
        """
        cdef Mat B = Mat()
        cdef PetscBool indef = PETSC_FALSE
        CHKERR( BVGetMatrix(self.bv, &B.mat, &indef) )
        if B.mat:
            CHKERR( PetscINCREF(B.obj) )
            return (B, toBool(indef))
        else:
            return (None, False)

    def setMatrix(self, Mat B or None: Mat | None, indef: bool = False) -> None:
        """
        Set the bilinear form to be used for inner products.

        Collective.

        Parameters
        ----------
        B
            The matrix of the inner product.
        indef
            Whether the matrix is indefinite.

        Notes
        -----
        This is used to specify a non-standard inner product, whose matrix
        representation is given by ``B``. Then, all inner products required
        during orthogonalization are computed as :math:`(x,y)_B=y^*Bx` rather
        than the standard form :math:`(x,y)=y^*x`.

        Matrix ``B`` must be real symmetric (or complex Hermitian). A genuine
        inner product requires that ``B`` is also positive (semi-)definite.
        However, we also allow for an indefinite ``B`` (setting ``indef=True``),
        in which case the orthogonalization uses an indefinite inner product.

        This affects operations `dot()`, `norm()`, `orthogonalize()`, and
        variants.

        Omitting ``B`` has the same effect as if the identity matrix was passed.

        See Also
        --------
        getMatrix, slepc.BVSetMatrix
        """
        cdef PetscMat m = <PetscMat>NULL if B is None else B.mat
        cdef PetscBool tval = PETSC_TRUE if indef else PETSC_FALSE
        CHKERR( BVSetMatrix(self.bv, m, tval) )

    def applyMatrix(self, Vec x, Vec y) -> None:
        """
        Multiply a vector with the matrix associated to the bilinear form.

        Neighbor-wise collective.

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

        See Also
        --------
        setMatrix, slepc.BVApplyMatrix
        """
        CHKERR( BVApplyMatrix(self.bv, x.vec, y.vec) )

    def setActiveColumns(self, l: int, k: int) -> None:
        """
        Set the columns that will be involved in operations.

        Logically collective.

        Parameters
        ----------
        l
            The leading number of columns.
        k
            The active number of columns.

        Notes
        -----
        In operations such as `mult()` or `dot()`, only the first ``k`` columns
        are considered. This is useful when the BV is filled from left to right,
        so the last ``m-k`` columns do not have relevant information.

        Also in operations such as `mult()` or `dot()`, the first ``l`` columns
        are normally not included in the computation.

        In orthogonalization operations, the first ``l`` columns are treated
        differently, they participate in the orthogonalization but the computed
        coefficients are not stored.

        Use `CURRENT` to leave any of the values unchanged. Use `DETERMINE`
        to set ``l`` to the minimum value (``0``) and ``k`` to the maximum (``m``).

        See Also
        --------
        getActiveColumns, setSizes, slepc.BVSetActiveColumns
        """
        cdef PetscInt ival1 = asInt(l)
        cdef PetscInt ival2 = asInt(k)
        CHKERR( BVSetActiveColumns(self.bv, ival1, ival2) )

    def getActiveColumns(self) -> tuple[int, int]:
        """
        Get the current active dimensions.

        Not collective.

        Returns
        -------
        l: int
            The leading number of columns.
        k: int
            The active number of columns.

        See Also
        --------
        setActiveColumns, slepc.BVGetActiveColumns
        """
        cdef PetscInt l=0, k=0
        CHKERR( BVGetActiveColumns(self.bv, &l, &k) )
        return (toInt(l), toInt(k))

    def scaleColumn(self, j: int, alpha: Scalar) -> None:
        """
        Scale a column of a BV.

        Logically collective.

        Parameters
        ----------
        j
            column index to be scaled.
        alpha
            scaling factor.

        See Also
        --------
        scale, slepc.BVScaleColumn
        """
        cdef PetscInt ival = asInt(j)
        cdef PetscScalar sval = asScalar(alpha)
        CHKERR( BVScaleColumn(self.bv, ival, sval) )

    def scale(self, alpha: Scalar) -> None:
        """
        Multiply the entries by a scalar value.

        Logically collective.

        Parameters
        ----------
        alpha
            scaling factor.

        Notes
        -----
        All active columns (except the leading ones) are scaled.

        See Also
        --------
        scaleColumn, setActiveColumns, slepc.BVScale
        """
        cdef PetscScalar sval = asScalar(alpha)
        CHKERR( BVScale(self.bv, sval) )

    def insertVec(self, j: int, Vec w) -> None:
        """
        Insert a vector into the specified column.

        Logically collective.

        Parameters
        ----------
        j
            The column to be overwritten.
        w
            The vector to be copied.

        See Also
        --------
        insertVecs, slepc.BVInsertVec
        """
        cdef PetscInt ival = asInt(j)
        CHKERR( BVInsertVec(self.bv, ival, w.vec) )

    def insertVecs(self, s: int, W: Vec | list[Vec], orth: bool = False) -> int:
        """
        Insert a set of vectors into the specified columns.

        Collective.

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
        Copies the contents of vectors ``W`` into the BV columns ``s:s+n``,
        where ``n`` is the length of ``W``. If ``orth`` is set, then the
        vectors are copied one by one and then orthogonalized against the
        previous one. If any of them is linearly dependent then it is
        discarded and the not counted in the return value.

        See Also
        --------
        insertVec, orthogonalizeColumn, slepc.BVInsertVecs
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

        Collective.

        Parameters
        ----------
        C
            Set of vectors to be inserted as constraints.

        Returns
        -------
        int
            Number of linearly independent constraints.

        Notes
        -----
        The constraints are relevant only during orthogonalization. Constraint
        vectors span a subspace that is deflated in every orthogonalization
        operation, so they are intended for removing those directions from the
        orthogonal basis computed in regular BV columns.

        Constraints are not stored in regular columns, but in a special part of
        the storage. They can be accessed with negative indices in
        `getColumn()`.

        This operation is DESTRUCTIVE, meaning that all data contained in the
        columns of the BV is lost. This is typically invoked just after creating
        the BV. Once a set of constraints has been set, it is not allowed to
        call this function again.

        The vectors are copied one by one and then orthogonalized against the
        previous ones. If any of them is linearly dependent then it is discarded
        and not counted in the return value. The behavior is similar to
        `insertVecs()`.

        See Also
        --------
        insertVecs, setNumConstraints, slepc.BVInsertConstraints
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
        Set the number of constraints.

        Logically collective.

        Parameters
        ----------
        nc
            The number of constraints.
        Notes
        -----
        This function sets the number of constraints to ``nc`` and marks all
        remaining columns as regular. Normal usage would be to call
        `insertConstraints()` instead.

        If ``nc`` is smaller than the previously set value, then some of the
        constraints are discarded. In particular, using ``nc=0`` removes all
        constraints preserving the content of regular columns.

        See Also
        --------
        insertConstraints, getNumConstraints, slepc.BVSetNumConstraints
        """
        cdef PetscInt val = asInt(nc)
        CHKERR( BVSetNumConstraints(self.bv, val) )

    def getNumConstraints(self) -> int:
        """
        Get the number of constraints.

        Not collective.

        Returns
        -------
        int
            The number of constraints.

        See Also
        --------
        insertConstraints, setNumConstraints, slepc.BVGetNumConstraints
        """
        cdef PetscInt val = 0
        CHKERR( BVGetNumConstraints(self.bv, &val) )
        return toInt(val)

    def createVec(self) -> Vec:
        """
        Create a vector with the type and dimensions of the columns of the BV.

        Collective.

        Returns
        -------
        petsc4py.PETSc.Vec
            New vector.

        See Also
        --------
        createMat, setVecType, slepc.BVCreateVec
        """
        cdef Vec v = Vec()
        CHKERR( BVCreateVec(self.bv, &v.vec) )
        return v

    def setVecType(self, vec_type: petsc4py.PETSc.Vec.Type | str) -> None:
        """
        Set the vector type to be used when creating vectors via `createVec()`.

        Collective.

        Parameters
        ----------
        vec_type
            Vector type used when creating vectors with `createVec`.

        Notes
        -----
        This is not needed if the BV object is set up with `setSizesFromVec()`,
        but may be required in the case of `setSizes()` if one wants to work
        with non-standard vectors.

        See Also
        --------
        createVec, getVecType, setSizes, setSizesFromVec, slepc.BVSetVecType
        """
        cdef PetscVecType cval = NULL
        vec_type = str2bytes(vec_type, &cval)
        CHKERR( BVSetVecType(self.bv, cval) )

    def getVecType(self) -> str:
        """
        Get the vector type used when creating vectors via `createVec()`.

        Not collective.

        Returns
        -------
        str
            The vector type.

        See Also
        --------
        createVec, setVecType, slepc.BVGetVecType
        """
        cdef PetscVecType cval = NULL
        CHKERR( BVGetVecType(self.bv, &cval) )
        return bytes2str(cval)

    def copyVec(self, j: int, Vec v) -> None:
        """
        Copy one of the columns of a basis vectors object into a vector.

        Logically collective.

        Parameters
        ----------
        j
            The column index to be copied.
        v
            A vector.

        Notes
        -----
        The BV and ``v`` must be distributed in the same manner; local copies
        are done.

        See Also
        --------
        copy, copyColumn, slepc.BVCopyVec
        """
        cdef PetscInt ival = asInt(j)
        CHKERR( BVCopyVec(self.bv, ival, v.vec) )

    def copyColumn(self, j: int, i: int) -> None:
        """
        Copy the values from one of the columns to another one.

        Logically collective.

        Parameters
        ----------
        j
            The index of the source column.
        i
            The index of the destination column.

        See Also
        --------
        copy, copyVec, slepc.BVCopyColumn
        """
        cdef PetscInt ival1 = asInt(j)
        cdef PetscInt ival2 = asInt(i)
        CHKERR( BVCopyColumn(self.bv, ival1, ival2) )

    def setDefiniteTolerance(self, deftol: float) -> None:
        """
        Set the tolerance to be used when checking a definite inner product.

        Logically collective.

        Parameters
        ----------
        deftol
            The tolerance.

        Notes
        -----
        When using a non-standard inner product, see `setMatrix()`, the solver
        needs to compute :math:`\sqrt{z^*B z}` for various vectors :math:`z`.
        If the inner product has not been declared indefinite, the value
        :math:`z^*B z` must be positive, but due to rounding error a tiny value
        may become negative. A tolerance is used to detect this situation.
        Likewise, in complex arithmetic :math:`z^*B z` should be real, and we
        use the same tolerance to check whether a nonzero imaginary part can be
        considered negligible.

        See Also
        --------
        setMatrix, getDefiniteTolerance, slepc.BVSetDefiniteTolerance
        """
        cdef PetscReal val = asReal(deftol)
        CHKERR( BVSetDefiniteTolerance(self.bv, val) )

    def getDefiniteTolerance(self) -> float:
        """
        Get the tolerance to be used when checking a definite inner product.

        Not collective.

        Returns
        -------
        float
            The tolerance.

        See Also
        --------
        setDefiniteTolerance, slepc.BVGetDefiniteTolerance
        """
        cdef PetscReal val = 0
        CHKERR( BVGetDefiniteTolerance(self.bv, &val) )
        return toReal(val)

    def dotVec(self, Vec v) -> ArrayScalar:
        """
        Dot products of a vector against all the column vectors of the BV.

        Collective.

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
        This is analogue to ``Vec.mDot()``, but using `BV` to represent a
        collection of vectors ``X``. The result is :math:`m = X^* v`, so
        :math:`m_i` is equal to :math:`x_j^* v`. Note that here :math:`X`
        is transposed as opposed to `dot()`.

        If a non-standard inner product has been specified with `setMatrix()`,
        then the result is :math:`m = X^* B v`.

        See Also
        --------
        dot, dotColumn, setMatrix, slepc.BVDotVec
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
        Dot products of a column against all the column vectors of a BV.

        Collective.

        Parameters
        ----------
        j
            The index of the column.

        Returns
        -------
        ArrayScalar
            The computed values.

        Notes
        -----
        This operation is equivalent to `dotVec()` but it uses column ``j`` of
        the BV rather than taking a vector as an argument. The number of active
        columns of the BV is set to ``j`` before the computation, and restored
        afterwards. If the BV has leading columns specified, then these columns
        do not participate in the computation. Therefore, the length of the
        returned array will be ``j`` minus the number of leading columns.

        See Also
        --------
        dot, dotVec, slepc.BVDotColumn
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
        Get a vector with the entries of the column of the BV object.

        Logically collective.

        Parameters
        ----------
        j
            The index of the requested column.

        Returns
        -------
        petsc4py.PETSc.Vec
            The vector containing the ``j``-th column.

        Notes
        -----
        Modifying the returned vector will change the BV entries as well.

        The returned vector must not be destroyed, `restoreColumn()` must be
        called when it is no longer needed. At most, two columns can be
        fetched, that is, this function can only be called twice before the
        corresponding `restoreColumn()` is invoked.

        A negative index ``j`` selects the ``i``-th constraint, where
        ``i=-j``. Constraints should not be modified.

        See Also
        --------
        restoreColumn, insertConstraints, slepc.BVGetColumn
        """
        cdef Vec v = Vec()
        cdef PetscInt ival = asInt(j)
        CHKERR( BVGetColumn(self.bv, j, &v.vec) )
        CHKERR( PetscINCREF(v.obj) )
        return v

    def restoreColumn(self, j: int, Vec v) -> None:
        """
        Restore a column obtained with `getColumn()`.

        Logically collective.

        Parameters
        ----------
        j
            The index of the requested column.
        v
            The vector obtained with `getColumn()`.

        Notes
        -----
        The arguments must match the corresponding call to `getColumn()`.

        See Also
        --------
        getColumn, slepc.BVRestoreColumn
        """
        cdef PetscInt ival = asInt(j)
        CHKERR( PetscObjectDereference(<PetscObject>v.vec) )
        CHKERR( BVRestoreColumn(self.bv, ival, &v.vec) )

    def getMat(self) -> Mat:
        """
        Get a matrix of dense type that shares the memory of the BV object.

        Collective.

        Returns
        -------
        petsc4py.PETSc.Mat
            The matrix.

        Notes
        -----
        The returned matrix contains only the active columns. If the content
        of the matrix is modified, these changes are also done in the BV
        object. The user must call `restoreMat()` when no longer needed.

        This operation implies a call to `getArray()`, which may result in
        data copies.

        See Also
        --------
        restoreMat, createMat, getArray, slepc.BVGetMat
        """
        cdef Mat A = Mat()
        CHKERR( BVGetMat(self.bv, &A.mat) )
        CHKERR( PetscINCREF(A.obj) )
        return A

    def restoreMat(self, Mat A) -> None:
        """
        Restore the matrix obtained with `getMat()`.

        Logically collective.

        Parameters
        ----------
        A
            The matrix obtained with `getMat()`.

        Notes
        -----
        A call to this function must match a previous call of `getMat()`.
        The effect is that the contents of the matrix are copied back to the
        BV internal data structures.

        See Also
        --------
        getMat, slepc.BVRestoreMat
        """
        CHKERR( PetscObjectDereference(<PetscObject>A.mat) )
        CHKERR( BVRestoreMat(self.bv, &A.mat) )

    def dot(self, BV Y) -> Mat:
        """
        Compute the 'block-dot' product of two basis vectors objects.

        Collective.

        :math:`M = Y^* X` :math:`(m_{ij} = y_i^* x_j)` or
        :math:`M = Y^* B X`

        Parameters
        ----------
        Y
            Left basis vectors, can be the same as self, giving
            :math:`M = X^* X`.

        Returns
        -------
        petsc4py.PETSc.Mat
            The resulting matrix.

        Notes
        -----
        This is the generalization of ``Vec.dot()`` for a collection of
        vectors, :math:`M = Y^* X`. The result is a matrix :math:`M` whose
        entry :math:`m_{ij}` is equal to :math:`y_i^* x_j`
        (where :math:`y_i^*` denotes the conjugate transpose of :math:`y_i`).

        :math:`X` and :math:`Y` can be the same object.

        If a non-standard inner product has been specified with `setMatrix()`,
        then the result is :math:`M = Y^* B X`. In this case, both
        :math:`X` and :math:`Y` must have the same associated matrix.

        Only rows (resp. columns) of :math:`M` starting from :math:`l_y` (resp.
        :math:`l_x`) are computed, where :math:`l_y` (resp. :math:`l_x`) is the
        number of leading columns of :math:`Y` (resp. :math:`X`).

        See Also
        --------
        dotVec, dotColumn, setActiveColumns, setMatrix, slepc.BVDot
        """
        cdef BV X = self
        cdef PetscInt ky=0, kx=0
        CHKERR( BVGetActiveColumns(Y.bv, NULL, &ky) )
        CHKERR( BVGetActiveColumns(X.bv, NULL, &kx) )
        cdef Mat M = Mat().createDense((ky, kx), comm=COMM_SELF).setUp()
        CHKERR( BVDot(X.bv, Y.bv, M.mat) )
        return M

    def matProject(self, Mat A: Mat | None, BV Y) -> Mat:
        """
        Compute the projection of a matrix onto a subspace.

        Collective.

        :math:`M = Y^* A X`

        Parameters
        ----------
        A
            Matrix to be projected.
        Y
            Left basis vectors, can be the same as self, giving
            :math:`M = X^* A X`.

        Returns
        -------
        petsc4py.PETSc.Mat
            Projection of the matrix ``A`` onto the subspace.

        Notes
        -----
        If ``A`` is ``None``, then it is assumed that the BV already
        contains :math:`AX`.

        This operation is similar to `dot()`, with important differences.
        The goal is to compute the matrix resulting from the orthogonal
        projection of ``A`` onto the subspace spanned by the columns of
        the BV, :math:`M = X^*AX`, or the oblique projection onto the BV
        along the second one ``Y``, :math:`M = Y^*AX`.

        A difference with respect to `dot()` is that the standard inner
        product is always used, regardless of a non-standard inner product
        being specified with `setMatrix()`.

        See Also
        --------
        dot, setActiveColumns, setMatrix, slepc.BVMatProject
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
        Compute the matrix-vector product for each column, :math:`Y = A V`.

        Neighbor-wise collective.

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
        If ``Y`` is ``None`` a new BV is created.

        It is possible to choose whether the computation is done column by column
        or as a dense matrix-matrix product with `setMatMultMethod()`.

        See Also
        --------
        copy, matMultColumn, matMultTranspose, setMatMultMethod, slepc.BVMatMult
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

    def matMultTranspose(self, Mat A, BV Y=None) -> BV:
        """
        Pre-multiplication with the transpose of a matrix.

        Neighbor-wise collective.

        :math:`Y = A^T V`.

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
        If ``Y`` is ``None`` a new BV is created.

        See Also
        --------
        matMult, matMultTransposeColumn, slepc.BVMatMultTranspose
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
        CHKERR( BVMatMultTranspose(self.bv, A.mat, Y.bv) )
        return Y

    def matMultHermitianTranspose(self, Mat A, BV Y=None) -> BV:
        """
        Pre-multiplication with the conjugate transpose of a matrix.

        Neighbor-wise collective.

        :math:`Y = A^* V`.

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
        If ``Y`` is ``None`` a new BV is created.

        See Also
        --------
        matMult, matMultHermitianTransposeColumn, slepc.BVMatMultHermitianTranspose
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
        Mat-vec product for a column, storing the result in the next column.

        Neighbor-wise collective.

        :math:`v_{j+1} = A v_j`.

        Parameters
        ----------
        A
            The matrix.
        j
            Index of column.

        See Also
        --------
        matMult, slepc.BVMatMultColumn
        """
        cdef PetscInt ival = asInt(j)
        CHKERR( BVMatMultColumn(self.bv, A.mat, ival) )

    def matMultTransposeColumn(self, Mat A, j: int) -> None:
        """
        Transpose matrix-vector product for a specified column.

        Neighbor-wise collective.

        Store the result in the next column: :math:`v_{j+1} = A^T v_j`.

        Parameters
        ----------
        A
            The matrix.
        j
            Index of column.

        See Also
        --------
        matMultColumn, slepc.BVMatMultTransposeColumn
        """
        cdef PetscInt ival = asInt(j)
        CHKERR( BVMatMultTransposeColumn(self.bv, A.mat, ival) )

    def matMultHermitianTransposeColumn(self, Mat A, j: int) -> None:
        """
        Conjugate-transpose matrix-vector product for a specified column.

        Neighbor-wise collective.

        Store the result in the next column: :math:`v_{j+1} = A^* v_j`.

        Parameters
        ----------
        A
            The matrix.
        j
            Index of column.

        See Also
        --------
        matMultColumn, slepc.BVMatMultHermitianTransposeColumn
        """
        cdef PetscInt ival = asInt(j)
        CHKERR( BVMatMultHermitianTransposeColumn(self.bv, A.mat, ival) )

    def mult(self, delta: Scalar, gamma: Scalar, BV X, Mat Q or None: Mat | None) -> None:
        """
        Compute :math:`Y = \gamma Y + \delta X Q`.

        Logically collective.

        Parameters
        ----------
        delta
            Coefficient that multiplies ``X``.
        gamma
            Coefficient that multiplies self (``Y``).
        X
            Input basis vectors.
        Q
            Input matrix, if not given the identity matrix is assumed.

        Notes
        -----
        ``X`` must be different from self (``Y``). The case ``X=Y`` can be
        addressed with `multInPlace()`.

        See Also
        --------
        multVec, multColumn, multInPlace, slepc.BVMult
        """
        cdef PetscScalar sval1 = asScalar(delta)
        cdef PetscScalar sval2 = asScalar(gamma)
        cdef PetscMat Qmat = <PetscMat>NULL if Q is None else Q.mat
        CHKERR( BVMult(self.bv, sval1, sval2, X.bv, Qmat) )

    def multInPlace(self, Mat Q, s: int, e: int) -> None:
        """
        Update a set of vectors as :math:`V(:,s:e-1) = V Q(:,s:e-1)`.

        Logically collective.

        Parameters
        ----------
        Q
            A sequential dense matrix.
        s
            First column to be overwritten.
        e
            Last column to be overwritten.

        See Also
        --------
        mult, multVec, slepc.BVMultInPlace
        """
        cdef PetscInt ival1 = asInt(s)
        cdef PetscInt ival2 = asInt(e)
        CHKERR( BVMultInPlace(self.bv, Q.mat, ival1, ival2) )

    def multColumn(self, delta: Scalar, gamma: Scalar, j: int, q: Sequence[Scalar]) -> None:
        """
        Compute :math:`y = \gamma y + \delta X q`.

        Logically collective.

        Compute :math:`y = \gamma y + \delta X q`, where
        :math:`y` is the ``j``-th column.

        Parameters
        ----------
        delta
            Coefficient that multiplies self (``X``).
        gamma
            Coefficient that multiplies :math:`y`.
        j
            The column index.
        q
            Input coefficients.

        See Also
        --------
        mult, multVec, multInPlace, slepc.BVMultColumn
        """
        cdef PetscScalar sval1 = asScalar(delta)
        cdef PetscScalar sval2 = asScalar(gamma)
        cdef PetscInt ival = asInt(j)
        cdef PetscInt nq = 0
        cdef PetscScalar* qval = NULL
        cdef tmp = iarray_s(q, &nq, &qval)
        cdef PetscInt l=0, k=0
        CHKERR( BVGetActiveColumns(self.bv, &l, &k) )
        assert nq == k-l
        CHKERR( BVMultColumn(self.bv, sval1, sval2, ival, qval) )

    def multVec(self, delta: Scalar, gamma: Scalar, Vec y, q: Sequence[Scalar]) -> None:
        """
        Compute :math:`y = \gamma y + \delta X q`.

        Logically collective.

        Parameters
        ----------
        delta
            Coefficient that multiplies self (``X``).
        gamma
            Coefficient that multiplies ``y``.
        y
            Input/output vector.
        q
            Input coefficients.

        See Also
        --------
        mult, multColumn, multInPlace, slepc.BVMultVec
        """
        cdef PetscScalar sval1 = asScalar(delta)
        cdef PetscScalar sval2 = asScalar(gamma)
        cdef PetscInt nq = 0
        cdef PetscScalar* qval = NULL
        cdef tmp = iarray_s(q, &nq, &qval)
        cdef PetscInt l=0, k=0
        CHKERR( BVGetActiveColumns(self.bv, &l, &k) )
        assert nq == k-l
        CHKERR( BVMultVec(self.bv, sval1, sval2, y.vec, qval) )

    def normColumn(self, j: int, norm_type: NormType | None = None) -> float:
        """
        Compute the vector norm of a selected column.

        Collective.

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
        The norm of :math:`v_j` is computed (``NORM_1``, ``NORM_2``, or
        ``NORM_INFINITY``).

        If a non-standard inner product has been specified with `setMatrix()`,
        then the returned value is :math:`\sqrt{v_j^* B v_j}`,
        where :math:`B` is the inner product matrix (argument 'norm_type' is
        ignored).

        See Also
        --------
        norm, setMatrix, slepc.BVNormColumn
        """
        cdef PetscNormType ntype = PETSC_NORM_2
        if norm_type is not None: ntype = norm_type
        cdef PetscReal norm = 0
        CHKERR( BVNormColumn(self.bv, j, ntype, &norm) )
        return toReal(norm)

    def norm(self, norm_type: NormType | None = None) -> float:
        """
        Compute the matrix norm of the BV.

        Collective.

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
        matrix. The allowed norms are ``NORM_1``, ``NORM_FROBENIUS``, and
        ``NORM_INFINITY``.

        This operation fails if a non-standard inner product has been specified
        with `setMatrix()`.

        See Also
        --------
        normColumn, setMatrix, slepc.BVNorm
        """
        cdef PetscNormType ntype = PETSC_NORM_FROBENIUS
        if norm_type is not None: ntype = norm_type
        cdef PetscReal norm = 0
        CHKERR( BVNorm(self.bv, ntype, &norm) )
        return toReal(norm)

    def resize(self, m: int, copy: bool = True) -> None:
        """
        Change the number of columns.

        Collective.

        Parameters
        ----------
        m
            The new number of columns.
        copy
            A flag indicating whether current values should be kept.

        Notes
        -----
        Internal storage is reallocated. If ``copy`` is ``True``, then the
        contents are copied to the leading part of the new space.

        See Also
        --------
        setSizes, setSizesFromVec, slepc.BVResize
        """
        cdef PetscInt ival = asInt(m)
        cdef PetscBool tval = PETSC_TRUE if copy else PETSC_FALSE
        CHKERR( BVResize(self.bv, ival, tval) )

    def setRandom(self) -> None:
        """
        Set the active columns of the BV to random numbers.

        Logically collective.

        Notes
        -----
        All active columns (except the leading ones) are modified.

        See Also
        --------
        setRandomContext, setRandomColumn, setRandomNormal, slepc.BVSetRandom
        """
        CHKERR( BVSetRandom(self.bv) )

    def setRandomNormal(self) -> None:
        """
        Set the active columns of the BV to normal random numbers.

        Logically collective.

        Notes
        -----
        All active columns (except the leading ones) are modified.

        See Also
        --------
        setRandomContext, setRandom, setRandomSign, slepc.BVSetRandomNormal
        """
        CHKERR( BVSetRandomNormal(self.bv) )

    def setRandomSign(self) -> None:
        """
        Set the entries of a BV to values 1 or -1 with equal probability.

        Logically collective.

        Notes
        -----
        All active columns (except the leading ones) are modified.

        See Also
        --------
        setRandomContext, setRandom, setRandomNormal, slepc.BVSetRandomSign
        """
        CHKERR( BVSetRandomSign(self.bv) )

    def setRandomColumn(self, j: int) -> None:
        """
        Set one column of the BV to random numbers.

        Logically collective.

        Parameters
        ----------
        j
            Column index to be set.

        See Also
        --------
        setRandomContext, setRandom, setRandomNormal, slepc.BVSetRandomColumn
        """
        cdef PetscInt ival = asInt(j)
        CHKERR( BVSetRandomColumn(self.bv, ival) )

    def setRandomCond(self, condn: float) -> None:
        """
        Set the columns of a BV to random numbers.

        Logically collective.

        The generated matrix has a prescribed condition number.

        Parameters
        ----------
        condn
            Condition number.

        See Also
        --------
        setRandomContext, setRandomSign, setRandomNormal, slepc.BVSetRandomCond
        """
        cdef PetscReal rval = asReal(condn)
        CHKERR( BVSetRandomCond(self.bv, rval) )

    def setRandomContext(self, Random rnd) -> None:
        """
        Set the `petsc4py.PETSc.Random` object associated with the BV.

        Collective.

        To be used in operations that need random numbers.

        Parameters
        ----------
        rnd
            The random number generator context.

        See Also
        --------
        getRandomContext, setRandom, setRandomColumn, slepc.BVSetRandomContext
        """
        CHKERR( BVSetRandomContext(self.bv, rnd.rnd) )

    def getRandomContext(self) -> Random:
        """
        Get the `petsc4py.PETSc.Random` object associated with the BV.

        Collective.

        Returns
        -------
        petsc4py.PETSc.Random
            The random number generator context.

        See Also
        --------
        setRandomContext, slepc.BVGetRandomContext
        """
        cdef Random rnd = Random()
        CHKERR( BVGetRandomContext(self.bv, &rnd.rnd) )
        CHKERR( PetscINCREF(rnd.obj) )
        return rnd

    def orthogonalizeVec(self, Vec v) -> tuple[float, bool]:
        """
        Orthogonalize a vector with respect to all active columns.

        Collective.

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
        This function applies an orthogonal projector to project vector
        :math:`v` onto the orthogonal complement of the span of the columns
        of the BV.

        This routine does not normalize the resulting vector.

        See Also
        --------
        orthogonalizeColumn, setOrthogonalization slepc.BVOrthogonalizeVec
        """
        cdef PetscReal norm = 0
        cdef PetscBool ldep = PETSC_FALSE
        CHKERR( BVOrthogonalizeVec(self.bv, v.vec, NULL, &norm, &ldep) )
        return (toReal(norm), toBool(ldep))

    def orthogonalizeColumn(self, j: int) -> tuple[float, bool]:
        """
        Orthogonalize a column vector with respect to the previous ones.

        Collective.

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
        This function applies an orthogonal projector to project vector
        :math:`v_j` onto the orthogonal complement of the span of the columns
        :math:`V[0..j-1]`, where :math:`V[.]` are the vectors of the BV.
        The columns :math:`V[0..j-1]` are assumed to be mutually orthonormal.

        This routine does not normalize the resulting vector.

        See Also
        --------
        orthogonalizeVec, setOrthogonalization slepc.BVOrthogonalizeColumn
        """
        cdef PetscInt ival = asInt(j)
        cdef PetscReal norm = 0
        cdef PetscBool ldep = PETSC_FALSE
        CHKERR( BVOrthogonalizeColumn(self.bv, ival, NULL, &norm, &ldep) )
        return (toReal(norm), toBool(ldep))

    def orthonormalizeColumn(self, j: int, replace: bool = False) -> tuple[float, bool]:
        """
        Orthonormalize a column vector with respect to the previous ones.

        Collective.

        This is equivalent to a call to `orthogonalizeColumn()` followed by a
        call to `scaleColumn()` with the reciprocal of the norm.

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

        See Also
        --------
        orthogonalizeColumn, setOrthogonalization slepc.BVOrthonormalizeColumn
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
        Orthogonalize all columns (except leading ones) (QR decomposition).

        Collective.

        Parameters
        ----------
        R
            A sequential dense matrix.

        Notes
        -----
        The output satisfies :math:`V_0 = V R` (where :math:`V_0` represent the
        input :math:`V`) and :math:`V^* V = I` (or :math:`V^*BV=I` if an inner
        product matrix :math:`B` has been specified with `setMatrix()`).

        See Also
        --------
        orthogonalizeColumn, setMatrix, setOrthogonalization, slepc.BVOrthogonalize
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
del BVSVDMethod

# -----------------------------------------------------------------------------
