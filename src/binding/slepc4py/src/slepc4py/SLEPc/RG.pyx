# -----------------------------------------------------------------------------

class RGType(object):
    """
    RG type.

    - `INTERVAL`: A (generalized) interval.
    - `POLYGON`: A polygonal region defined by its vertices.
    - `ELLIPSE`: An ellipse defined by its center, radius and vertical scale.
    - `RING`: A ring region.

    See Also
    --------
    slepc.RGType
    """
    INTERVAL   = S_(RGINTERVAL)
    POLYGON    = S_(RGPOLYGON)
    ELLIPSE    = S_(RGELLIPSE)
    RING       = S_(RGRING)

class RGQuadRule(object):
    """
    RG quadrature rule for contour integral methods.

    - `TRAPEZOIDAL`: Trapezoidal rule.
    - `CHEBYSHEV`:   Chebyshev points.

    See Also
    --------
    slepc.RGQuadRule
    """
    TRAPEZOIDAL = EPS_CISS_QUADRULE_TRAPEZOIDAL
    CHEBYSHEV   = EPS_CISS_QUADRULE_CHEBYSHEV

# -----------------------------------------------------------------------------

cdef class RG(Object):

    """
    Region.

    The `RG` package provides a way to define a region of the complex plane.
    This is used in various eigensolvers to specify where the wanted
    eigenvalues are located.
    """

    Type     = RGType
    QuadRule = RGQuadRule

    def __cinit__(self):
        self.obj = <PetscObject*> &self.rg
        self.rg = NULL

    def view(self, Viewer viewer=None) -> None:
        """
        Print the RG data structure.

        Collective.

        Parameters
        ----------
        viewer
            Visualization context; if not provided, the standard
            output is used.

        See Also
        --------
        slepc.RGView
        """
        cdef PetscViewer vwr = def_Viewer(viewer)
        CHKERR( RGView(self.rg, vwr) )

    def destroy(self) -> Self:
        """
        Destroy the RG object.

        Collective.

        See Also
        --------
        slepc.RGDestroy
        """
        CHKERR( RGDestroy(&self.rg) )
        self.rg = NULL
        return self

    def create(self, comm: Comm | None = None) -> Self:
        """
        Create the RG object.

        Collective.

        Parameters
        ----------
        comm
            MPI communicator; if not provided, it defaults to all processes.

        See Also
        --------
        slepc.RGCreate
        """
        cdef MPI_Comm ccomm = def_Comm(comm, SLEPC_COMM_DEFAULT())
        cdef SlepcRG newrg = NULL
        CHKERR( RGCreate(ccomm, &newrg) )
        CHKERR( SlepcCLEAR(self.obj) ); self.rg = newrg
        return self

    def setType(self, rg_type: Type | str) -> None:
        """
        Set the type for the RG object.

        Logically collective.

        Parameters
        ----------
        rg_type
            The region type to be used.

        See Also
        --------
        getType, slepc.RGSetType
        """
        cdef SlepcRGType cval = NULL
        rg_type = str2bytes(rg_type, &cval)
        CHKERR( RGSetType(self.rg, cval) )

    def getType(self) -> str:
        """
        Get the RG type of this object.

        Not collective.

        Returns
        -------
        str
            The region type currently being used.

        See Also
        --------
        setType, slepc.RGGetType
        """
        cdef SlepcRGType rg_type = NULL
        CHKERR( RGGetType(self.rg, &rg_type) )
        return bytes2str(rg_type)

    def setOptionsPrefix(self, prefix: str | None = None) -> None:
        """
        Set the prefix used for searching for all RG options in the database.

        Logically collective.

        Parameters
        ----------
        prefix
            The prefix string to prepend to all RG option requests.

        Notes
        -----
        A hyphen (``-``) must NOT be given at the beginning of the
        prefix name.  The first character of all runtime options is
        AUTOMATICALLY the hyphen.

        See Also
        --------
        appendOptionsPrefix, getOptionsPrefix, slepc.RGGetOptionsPrefix
        """
        cdef const char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( RGSetOptionsPrefix(self.rg, cval) )

    def getOptionsPrefix(self) -> str:
        """
        Get the prefix used for searching for all RG options in the database.

        Not collective.

        Returns
        -------
        str
            The prefix string set for this RG object.

        See Also
        --------
        setOptionsPrefix, appendOptionsPrefix, slepc.RGGetOptionsPrefix
        """
        cdef const char *prefix = NULL
        CHKERR( RGGetOptionsPrefix(self.rg, &prefix) )
        return bytes2str(prefix)

    def appendOptionsPrefix(self, prefix: str | None = None) -> None:
        """
        Append to the prefix used for searching for all RG options in the database.

        Logically collective.

        Parameters
        ----------
        prefix
            The prefix string to prepend to all RG option requests.

        See Also
        --------
        setOptionsPrefix, getOptionsPrefix, slepc.RGAppendOptionsPrefix
        """
        cdef const char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( RGAppendOptionsPrefix(self.rg, cval) )

    def setFromOptions(self) -> None:
        """
        Set RG options from the options database.

        Collective.

        Notes
        -----
        To see all options, run your program with the ``-help``
        option.

        See Also
        --------
        setOptionsPrefix, slepc.RGSetFromOptions
        """
        CHKERR( RGSetFromOptions(self.rg) )

    #

    def isTrivial(self) -> bool:
        """
        Tell whether it is the trivial region (whole complex plane).

        Not collective.

        Returns
        -------
        bool
            ``True`` if the region is equal to the whole complex plane, e.g.,
            an interval region with all four endpoints unbounded or an
            ellipse with infinite radius.

        See Also
        --------
        checkInside, slepc.RGIsTrivial
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( RGIsTrivial(self.rg, &tval) )
        return toBool(tval)

    def isAxisymmetric(self, vertical: bool = False) -> bool:
        """
        Determine if the region is axisymmetric.

        Not collective.

        Determine if the region is symmetric with respect to the real or
        imaginary axis.

        Parameters
        ----------
        vertical
            ``True`` if symmetry must be checked against the vertical axis.

        Returns
        -------
        bool
            ``True`` if the region is axisymmetric.

        See Also
        --------
        canUseConjugates, slepc.RGIsAxisymmetric
        """
        cdef PetscBool val = asBool(vertical)
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( RGIsAxisymmetric(self.rg, val, &tval) )
        return toBool(tval)

    def getComplement(self) -> bool:
        """
        Get the flag indicating whether the region is complemented or not.

        Not collective.

        Returns
        -------
        bool
            Whether the region is complemented or not.

        See Also
        --------
        setComplement, slepc.RGGetComplement
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( RGGetComplement(self.rg, &tval) )
        return toBool(tval)

    def setComplement(self, comp: bool = True) -> None:
        """
        Set a flag to indicate that the region is the complement of the specified one.

        Logically collective.

        Parameters
        ----------
        comp
            Activate/deactivate the complementation of the region.

        See Also
        --------
        getComplement, slepc.RGSetComplement
        """
        cdef PetscBool tval = asBool(comp)
        CHKERR( RGSetComplement(self.rg, tval) )

    def setScale(self, sfactor: float = None) -> None:
        """
        Set the scaling factor to be used.

        Logically collective.

        Set the scaling factor to be used when checking that a
        point is inside the region and when computing the contour.

        Parameters
        ----------
        sfactor
            The scaling factor (default=1).

        See Also
        --------
        getScale, checkInside, computeContour, slepc.RGSetScale
        """
        cdef PetscReal rval = 1.0
        if sfactor is not None: rval = asReal(sfactor)
        CHKERR( RGSetScale(self.rg, rval) )

    def getScale(self) -> float:
        """
        Get the scaling factor.

        Not collective.

        Returns
        -------
        float
            The scaling factor.

        See Also
        --------
        setScale, slepc.RGGetScale
        """
        cdef PetscReal rval = 0
        CHKERR( RGGetScale(self.rg, &rval) )
        return toReal(rval)

    def checkInside(self, a: Sequence[complex]) -> ArrayInt:
        """
        Determine if a set of given points are inside the region or not.

        Not collective.

        Parameters
        ----------
        a
            The coordinates of the points.

        Returns
        -------
        ArrayInt
            Computed result for each point (1=inside, 0=on the contour, -1=outside).

        Notes
        -----
        If a scaling factor was set, the points are scaled before checking.

        See Also
        --------
        setScale, setComplement, slepc.RGCheckInside
        """
        cdef Py_ssize_t i = 0, n = len(a)
        cdef PetscScalar *ar = NULL, *ai = NULL
        cdef PetscInt *inside = NULL
        cdef tmp1 = allocate(<size_t>n*sizeof(PetscScalar),<void**>&ar)
        cdef tmp2
        if sizeof(PetscScalar) == sizeof(PetscReal):
            tmp2 = allocate(<size_t>n*sizeof(PetscScalar),<void**>&ai)
            for i in range(n):
                ar[i] = asComplexReal(a[i])
                ai[i] = asComplexImag(a[i])
        else:
            for i in range(n): ar[i] = asScalar(a[i])
        cdef tmp3 = allocate(<size_t>n*sizeof(PetscInt),<void**>&inside)
        CHKERR( RGCheckInside(self.rg, <PetscInt>n, ar, ai, inside) )
        return array_i(<PetscInt>n, inside)

    def computeContour(self, n: int) -> list[complex]:
        """
        Compute points on the contour of the region.

        Not collective.

        Compute the coordinates of several points lying on the contour
        of the region.

        Parameters
        ----------
        n
            The number of points to compute.

        Returns
        -------
        list of complex
            Computed points.

        See Also
        --------
        computeBoundingBox, setScale, slepc.RGComputeContour
        """
        cdef PetscInt k = asInt(n), i = 0
        cdef PetscScalar *cr = NULL, *ci = NULL
        cdef tmp1 = allocate(<size_t>k*sizeof(PetscScalar),<void**>&cr)
        cdef tmp2
        if sizeof(PetscScalar) == sizeof(PetscReal):
            tmp2 = allocate(<size_t>k*sizeof(PetscScalar),<void**>&ci)
        CHKERR( RGComputeContour(self.rg, k, cr, ci) )
        if sizeof(PetscScalar) == sizeof(PetscReal):
            return [toComplex(cr[i],ci[i]) for i from 0 <= i <k]
        else:
            return [toScalar(cr[i]) for i from 0 <= i <k]

    def computeBoundingBox(self) -> tuple[float, float, float, float]:
        """
        Compute box containing the region.

        Not collective.

        Determine the endpoints of a rectangle in the complex plane that
        contains the region.

        Returns
        -------
        a: float
            The left endpoint of the bounding box in the real axis.
        b: float
            The right endpoint of the bounding box in the real axis.
        c: float
            The bottom endpoint of the bounding box in the imaginary axis.
        d: float
            The top endpoint of the bounding box in the imaginary axis.

        See Also
        --------
        computeContour, setScale, slepc.RGComputeBoundingBox
        """
        cdef PetscReal a = 0, b = 0, c = 0, d = 0
        CHKERR( RGComputeBoundingBox(self.rg, &a, &b, &c, &d) )
        return (toReal(a), toReal(b), toReal(c), toReal(d))

    def canUseConjugates(self, realmats: bool = True) -> bool:
        """
        Half of integration points can be avoided (use their conjugates).

        Not collective.

        Used in contour integral methods to determine whether half of
        integration points can be avoided (use their conjugates).

        Parameters
        ----------
        realmats
            ``True`` if the problem matrices are real.

        Returns
        -------
        bool
            Whether it is possible to use conjugates.

        Notes
        -----
        If some integration points are the conjugates of other points, then the
        associated computational cost can be saved. This depends on the problem
        matrices being real and also the region being symmetric with respect to
        the horizontal axis. The result is ``false`` if using real arithmetic or
        in the case of a flat region (height equal to zero).

        See Also
        --------
        isAxisymmetric, slepc.RGCanUseConjugates
        """
        cdef PetscBool bval = asBool(realmats)
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( RGCanUseConjugates(self.rg, bval, &tval) )
        return toBool(tval)

    def computeQuadrature(self, quad: QuadRule, n: int) -> tuple[ArrayScalar, ArrayScalar, ArrayScalar]:
        """
        Compute the values of the parameters used in a quadrature rule.

        Not collective.

        Compute the values of the parameters used in a quadrature rule for a
        contour integral around the boundary of the region.

        Parameters
        ----------
        quad
            The type of quadrature.
        n
            The number of quadrature points to compute.

        Returns
        -------
        z: ArrayScalar
            Quadrature points.
        zn: ArrayScalar
            Normalized quadrature points.
        w: ArrayScalar
            Quadrature weights.

        Notes
        -----
        In complex scalars, the values returned in ``z`` are often the same as
        those computed by `computeContour()`, but this is not the case in real
        scalars where all output arguments are real.

        The computed values change for different quadrature rules.

        See Also
        --------
        computeContour, slepc.RGComputeQuadrature
        """
        cdef SlepcRGQuadRule val = quad
        cdef PetscInt k = asInt(n), i = 0
        cdef PetscScalar *z = NULL, *zn = NULL, *w = NULL
        cdef tmp1 = allocate(<size_t>k*sizeof(PetscScalar),<void**>&z)
        cdef tmp2 = allocate(<size_t>k*sizeof(PetscScalar),<void**>&zn)
        cdef tmp3 = allocate(<size_t>k*sizeof(PetscScalar),<void**>&w)
        CHKERR( RGComputeQuadrature(self.rg, val, k, z, zn, w) )
        return (array_s(k, z), array_s(k, zn), array_s(k, w))

    #

    def setEllipseParameters(self, center: Scalar, radius: float, vscale: float | None = None) -> None:
        """
        Set the parameters defining the ellipse region.

        Logically collective.

        Parameters
        ----------
        center
            The center.
        radius
            The radius.
        vscale
            The vertical scale.

        Notes
        -----
        When PETSc is built with real scalars, the center is restricted to a
        real value.

        See Also
        --------
        getEllipseParameters, slepc.RGEllipseSetParameters
        """
        cdef PetscScalar sval = asScalar(center)
        cdef PetscReal val1 = asReal(radius)
        cdef PetscReal val2 = 1.0
        if vscale is not None: val2 = asReal(vscale)
        CHKERR( RGEllipseSetParameters(self.rg, sval, val1, val2) )

    def getEllipseParameters(self) -> tuple[Scalar, float, float]:
        """
        Get the parameters that define the ellipse region.

        Not collective.

        Returns
        -------
        center: Scalar
            The center.
        radius: float
            The radius.
        vscale: float
            The vertical scale.

        See Also
        --------
        setEllipseParameters, slepc.RGEllipseGetParameters
        """
        cdef PetscScalar sval = 0
        cdef PetscReal val1 = 0
        cdef PetscReal val2 = 0
        CHKERR( RGEllipseGetParameters(self.rg, &sval, &val1, &val2) )
        return (toScalar(sval), toReal(val1), toReal(val2))

    def setIntervalEndpoints(self, a: float, b: float, c: float, d: float) -> None:
        """
        Set the parameters defining the interval region.

        Logically collective.

        Parameters
        ----------
        a
            The left endpoint in the real axis.
        b
            The right endpoint in the real axis.
        c
            The bottom endpoint in the imaginary axis.
        d
            The top endpoint in the imaginary axis.

        Notes
        -----
        The region is defined as :math:`[a,b] x [c,d]`. Particular cases are an
        interval on the real axis (:math:`c=d=0`), similarly for the imaginary
        axis (:math:`a=b=0`), the whole complex plane
        (:math:`a=-\infty,b=\infty,c=-\infty,d=\infty`), and so on.

        When PETSc is built with real scalars, the region must be symmetric with
        respect to the real axis.

        See Also
        --------
        getIntervalEndpoints, slepc.RGIntervalSetEndpoints
        """
        cdef PetscReal va = asReal(a)
        cdef PetscReal vb = asReal(b)
        cdef PetscReal vc = asReal(c)
        cdef PetscReal vd = asReal(d)
        CHKERR( RGIntervalSetEndpoints(self.rg, va, vb, vc, vd) )

    def getIntervalEndpoints(self) -> tuple[float, float, float, float]:
        """
        Get the parameters that define the interval region.

        Not collective.

        Returns
        -------
        a: float
            The left endpoint in the real axis.
        b: float
            The right endpoint in the real axis.
        c: float
            The bottom endpoint in the imaginary axis.
        d: float
            The top endpoint in the imaginary axis.

        See Also
        --------
        setIntervalEndpoints, slepc.RGIntervalGetEndpoints
        """
        cdef PetscReal va = 0
        cdef PetscReal vb = 0
        cdef PetscReal vc = 0
        cdef PetscReal vd = 0
        CHKERR( RGIntervalGetEndpoints(self.rg, &va, &vb, &vc, &vd) )
        return (toReal(va), toReal(vb), toReal(vc), toReal(vd))

    def setPolygonVertices(self, v: Sequence[float] | Sequence[Scalar]) -> None:
        """
        Set the vertices that define the polygon region.

        Logically collective.

        Parameters
        ----------
        v
            The vertices.

        See Also
        --------
        getPolygonVertices, slepc.RGPolygonSetVertices
        """
        cdef Py_ssize_t i = 0, n = len(v)
        cdef PetscScalar *vr = NULL, *vi = NULL
        cdef tmp1 = allocate(<size_t>n*sizeof(PetscScalar),<void**>&vr)
        cdef tmp2
        if sizeof(PetscScalar) == sizeof(PetscReal):
            tmp2 = allocate(<size_t>n*sizeof(PetscScalar),<void**>&vi)
            for i in range(n):
                vr[i] = asComplexReal(v[i])
                vi[i] = asComplexImag(v[i])
        else:
            for i in range(n): vr[i] = asScalar(v[i])
        CHKERR( RGPolygonSetVertices(self.rg, <PetscInt>n, vr, vi) )

    def getPolygonVertices(self) -> ArrayComplex:
        """
        Get the parameters that define the interval region.

        Not collective.

        Returns
        -------
        ArrayComplex
            The vertices.

        See Also
        --------
        setPolygonVertices, slepc.RGPolygonGetVertices
        """
        cdef PetscInt n = 0
        cdef PetscScalar *vr = NULL, *vi = NULL
        CHKERR( RGPolygonGetVertices(self.rg, &n, &vr, &vi) )
        if sizeof(PetscScalar) == sizeof(PetscReal):
            v = [toComplex(vr[i],vi[i]) for i from 0 <= i <n]
            CHKERR( PetscFree(vi) )
        else:
            v = [toScalar(vr[i]) for i from 0 <= i <n]
        CHKERR( PetscFree(vr) )
        return v

    def setRingParameters(
        self,
        center: Scalar,
        radius: float,
        vscale: float,
        start_ang: float,
        end_ang: float,
        width: float,
    ) -> None:
        """
        Set the parameters defining the ring region.

        Logically collective.

        Parameters
        ----------
        center
            The center.
        radius
            The radius.
        vscale
            The vertical scale.
        start_ang
            The right-hand side angle.
        end_ang
            The left-hand side angle.
        width
            The width of the ring.

        Notes
        -----
        The values of ``center``, ``radius`` and ``vscale`` have the same
        meaning as in the ellipse region. The ``start_ang`` and ``end_ang``
        define the span of the ring (by default it is the whole ring), while
        the ``width`` is the separation between the two concentric ellipses
        (above and below the radius by ``width/2``).

        The start and end angles are expressed as a fraction of the
        circumference. The allowed range is :math:`[0,\dots,1]`, with ``0``
        corresponding to 0 radians, ``0.25`` to :math:`\pi/2` radians, and so
        on. It is allowed to have ``start_ang`` > ``end_ang``, in which case
        the ring region crosses over the zero angle.

        When PETSc is built with real scalars, the center is restricted to a
        real value, and the start and end angles must be such that the region
        is symmetric with respect to the real axis.

        See Also
        --------
        getRingParameters, slepc.RGRingSetParameters
        """
        cdef PetscScalar sval = asScalar(center)
        cdef PetscReal val1 = asReal(radius)
        cdef PetscReal val2 = asReal(vscale)
        cdef PetscReal val3 = asReal(start_ang)
        cdef PetscReal val4 = asReal(end_ang)
        cdef PetscReal val5 = asReal(width)
        CHKERR( RGRingSetParameters(self.rg, sval, val1, val2, val3, val4, val5) )

    def getRingParameters(self) -> tuple[Scalar, float, float, float, float, float]:
        """
        Get the parameters that define the ring region.

        Not collective.

        Returns
        -------
        center: Scalar
            The center.
        radius: float
            The radius.
        vscale: float
            The vertical scale.
        start_ang: float
            The right-hand side angle.
        end_ang: float
            The left-hand side angle.
        width: float
            The width of the ring.

        See Also
        --------
        setRingParameters, slepc.RGRingGetParameters
        """
        cdef PetscScalar sval = 0
        cdef PetscReal val1 = 0
        cdef PetscReal val2 = 0
        cdef PetscReal val3 = 0
        cdef PetscReal val4 = 0
        cdef PetscReal val5 = 0
        CHKERR( RGRingGetParameters(self.rg, &sval, &val1, &val2, &val3, &val4, &val5) )
        return (toScalar(sval), toReal(val1), toReal(val2), toReal(val3), toReal(val4), toReal(val5))

    #

    property complement:
        """If the region is the complement of the specified one."""
        def __get__(self) -> bool:
            return self.getComplement()
        def __set__(self, value):
            self.setComplement(value)

    property scale:
        """The scaling factor to be used."""
        def __get__(self) -> float:
            return self.getScale()
        def __set__(self, value):
            self.setScale(value)

# -----------------------------------------------------------------------------

del RGType
del RGQuadRule

# -----------------------------------------------------------------------------
