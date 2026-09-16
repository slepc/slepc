cdef extern from * nogil:

    ctypedef char* SlepcRGType "const char*"
    SlepcRGType RGINTERVAL
    SlepcRGType RGPOLYGON
    SlepcRGType RGELLIPSE
    SlepcRGType RGRING

    ctypedef enum SlepcRGQuadRule "RGQuadRule":
        RG_QUADRULE_TRAPEZOIDAL
        RG_QUADRULE_CHEBYSHEV

    PetscErrorCode RGCreate(MPI_Comm, SlepcRG*)
    PetscErrorCode RGView(SlepcRG, PetscViewer)
    PetscErrorCode RGDestroy(SlepcRG*)
    PetscErrorCode RGSetType(SlepcRG, SlepcRGType)
    PetscErrorCode RGGetType(SlepcRG, SlepcRGType*)

    PetscErrorCode RGSetOptionsPrefix(SlepcRG, char[])
    PetscErrorCode RGGetOptionsPrefix(SlepcRG, char*[])
    PetscErrorCode RGAppendOptionsPrefix(SlepcRG, char[])
    PetscErrorCode RGSetFromOptions(SlepcRG)

    PetscErrorCode RGIsTrivial(SlepcRG, PetscBool*)
    PetscErrorCode RGIsAxisymmetric(SlepcRG, PetscBool, PetscBool*)
    PetscErrorCode RGSetComplement(SlepcRG, PetscBool)
    PetscErrorCode RGGetComplement(SlepcRG, PetscBool*)

    PetscErrorCode RGSetScale(SlepcRG, PetscReal)
    PetscErrorCode RGGetScale(SlepcRG, PetscReal*)
    PetscErrorCode RGPushScale(SlepcRG, PetscReal)
    PetscErrorCode RGPopScale(SlepcRG)

    PetscErrorCode RGCheckInside(SlepcRG, PetscInt, PetscScalar*, PetscScalar*, PetscInt*)
    PetscErrorCode RGComputeContour(SlepcRG, PetscInt, PetscScalar*, PetscScalar*)
    PetscErrorCode RGComputeBoundingBox(SlepcRG, PetscReal*, PetscReal*, PetscReal*, PetscReal*)
    PetscErrorCode RGCanUseConjugates(SlepcRG, PetscBool, PetscBool*)
    PetscErrorCode RGComputeQuadrature(SlepcRG, SlepcRGQuadRule, PetscInt, PetscScalar*, PetscScalar*, PetscScalar*)

    PetscErrorCode RGEllipseSetParameters(SlepcRG, PetscScalar, PetscReal, PetscReal)
    PetscErrorCode RGEllipseGetParameters(SlepcRG, PetscScalar*, PetscReal*, PetscReal*)
    PetscErrorCode RGIntervalSetEndpoints(SlepcRG, PetscReal, PetscReal, PetscReal, PetscReal)
    PetscErrorCode RGIntervalGetEndpoints(SlepcRG, PetscReal*, PetscReal*, PetscReal*, PetscReal*)
    PetscErrorCode RGPolygonSetVertices(SlepcRG, PetscInt, PetscScalar*, PetscScalar*)
    PetscErrorCode RGPolygonGetVertices(SlepcRG, PetscInt*, PetscScalar**, PetscScalar**)
    PetscErrorCode RGRingSetParameters(SlepcRG, PetscScalar, PetscReal, PetscReal, PetscReal, PetscReal, PetscReal)
    PetscErrorCode RGRingGetParameters(SlepcRG, PetscScalar*, PetscReal*, PetscReal*, PetscReal*, PetscReal*, PetscReal*)
