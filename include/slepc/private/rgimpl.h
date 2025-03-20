/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#pragma once

#include <slepcrg.h>
#include <slepc/private/slepcimpl.h>

/* SUBMANSEC = RG */

SLEPC_EXTERN PetscBool RGRegisterAllCalled;
SLEPC_EXTERN PetscErrorCode RGRegisterAll(void);

typedef struct _RGOps *RGOps;

struct _RGOps {
  PetscErrorCode (*istrivial)(RG,PetscBool*);
  PetscErrorCode (*computecontour)(RG,PetscInt,PetscScalar*,PetscScalar*);
  PetscErrorCode (*computebbox)(RG,PetscReal*,PetscReal*,PetscReal*,PetscReal*);
  PetscErrorCode (*computequadrature)(RG,RGQuadRule,PetscInt,PetscScalar*,PetscScalar*,PetscScalar*);
  PetscErrorCode (*checkinside)(RG,PetscReal,PetscReal,PetscInt*);
  PetscErrorCode (*isaxisymmetric)(RG,PetscBool,PetscBool*);
  PetscErrorCode (*setfromoptions)(RG,PetscOptionItems);
  PetscErrorCode (*view)(RG,PetscViewer);
  PetscErrorCode (*destroy)(RG);
};

struct _p_RG {
  PETSCHEADER(struct _RGOps);
  PetscBool   complement;    /* region is the complement of the specified one */
  PetscReal   sfactor;       /* scaling factor */
  PetscReal   osfactor;      /* old scaling factor, before RGPushScale */
  void        *data;
};

/* show an inf instead of PETSC_MAX_REAL */
#define RGShowReal(r) (double)((PetscAbsReal(r)>=PETSC_MAX_REAL)?10*(r):(r))
