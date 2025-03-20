/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   User interface for the region object in SLEPc
*/

#pragma once

#include <slepcsys.h>
#include <slepcrgtypes.h>

/* SUBMANSEC = RG */

SLEPC_EXTERN PetscErrorCode RGInitializePackage(void);
SLEPC_EXTERN PetscErrorCode RGFinalizePackage(void);

/*J
   RGType - String with the name of the region.

   Level: beginner

.seealso: RGSetType(), RG
J*/
typedef const char *RGType;
#define RGINTERVAL  "interval"
#define RGPOLYGON   "polygon"
#define RGELLIPSE   "ellipse"
#define RGRING      "ring"

/* Logging support */
SLEPC_EXTERN PetscClassId RG_CLASSID;

SLEPC_EXTERN PetscErrorCode RGCreate(MPI_Comm,RG*);
SLEPC_EXTERN PetscErrorCode RGSetType(RG,RGType);
SLEPC_EXTERN PetscErrorCode RGGetType(RG,RGType*);
SLEPC_EXTERN PetscErrorCode RGSetOptionsPrefix(RG,const char *);
SLEPC_EXTERN PetscErrorCode RGAppendOptionsPrefix(RG,const char *);
SLEPC_EXTERN PetscErrorCode RGGetOptionsPrefix(RG,const char *[]);
SLEPC_EXTERN PetscErrorCode RGSetFromOptions(RG);
SLEPC_EXTERN PetscErrorCode RGView(RG,PetscViewer);
SLEPC_EXTERN PetscErrorCode RGViewFromOptions(RG,PetscObject,const char[]);
SLEPC_EXTERN PetscErrorCode RGDestroy(RG*);

/*E
    RGQuadRule - determines how to compute the quadrature parameters

    Level: advanced

.seealso: RGComputeQuadrature()
E*/
typedef enum { RG_QUADRULE_TRAPEZOIDAL = 1,
               RG_QUADRULE_CHEBYSHEV   = 2 } RGQuadRule;

SLEPC_EXTERN PetscErrorCode RGIsTrivial(RG,PetscBool*);
SLEPC_EXTERN PetscErrorCode RGSetComplement(RG,PetscBool);
SLEPC_EXTERN PetscErrorCode RGGetComplement(RG,PetscBool*);
SLEPC_EXTERN PetscErrorCode RGSetScale(RG,PetscReal);
SLEPC_EXTERN PetscErrorCode RGGetScale(RG,PetscReal*);
SLEPC_EXTERN PetscErrorCode RGPushScale(RG,PetscReal);
SLEPC_EXTERN PetscErrorCode RGPopScale(RG);
SLEPC_EXTERN PetscErrorCode RGCheckInside(RG,PetscInt,PetscScalar*,PetscScalar*,PetscInt*);
SLEPC_EXTERN PetscErrorCode RGIsAxisymmetric(RG,PetscBool,PetscBool*);
SLEPC_EXTERN PetscErrorCode RGCanUseConjugates(RG,PetscBool,PetscBool*);
SLEPC_EXTERN PetscErrorCode RGComputeContour(RG,PetscInt,PetscScalar*,PetscScalar*);
SLEPC_EXTERN PetscErrorCode RGComputeBoundingBox(RG,PetscReal*,PetscReal*,PetscReal*,PetscReal*);
SLEPC_EXTERN PetscErrorCode RGComputeQuadrature(RG,RGQuadRule,PetscInt,PetscScalar*,PetscScalar*,PetscScalar*);

SLEPC_EXTERN PetscFunctionList RGList;
SLEPC_EXTERN PetscErrorCode RGRegister(const char[],PetscErrorCode(*)(RG));

/* --------- options specific to particular regions -------- */

SLEPC_EXTERN PetscErrorCode RGEllipseSetParameters(RG,PetscScalar,PetscReal,PetscReal);
SLEPC_EXTERN PetscErrorCode RGEllipseGetParameters(RG,PetscScalar*,PetscReal*,PetscReal*);

SLEPC_EXTERN PetscErrorCode RGIntervalSetEndpoints(RG,PetscReal,PetscReal,PetscReal,PetscReal);
SLEPC_EXTERN PetscErrorCode RGIntervalGetEndpoints(RG,PetscReal*,PetscReal*,PetscReal*,PetscReal*);

SLEPC_EXTERN PetscErrorCode RGPolygonSetVertices(RG,PetscInt,PetscScalar*,PetscScalar*);
SLEPC_EXTERN PetscErrorCode RGPolygonGetVertices(RG,PetscInt*,PetscScalar*[],PetscScalar*[]);

SLEPC_EXTERN PetscErrorCode RGRingSetParameters(RG,PetscScalar,PetscReal,PetscReal,PetscReal,PetscReal,PetscReal);
SLEPC_EXTERN PetscErrorCode RGRingGetParameters(RG,PetscScalar*,PetscReal*,PetscReal*,PetscReal*,PetscReal*,PetscReal*);
