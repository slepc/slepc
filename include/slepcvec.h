/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   User interface for various vector operations added in SLEPc
*/

#pragma once

#include <slepcsys.h>

/* SUBMANSEC = sys */

/* VecComp: Vec composed of several smaller Vecs */
#define VECCOMP  "comp"

SLEPC_EXTERN PetscErrorCode VecCreateComp(MPI_Comm,PetscInt*,PetscInt,VecType,Vec,Vec*);
SLEPC_EXTERN PetscErrorCode VecCreateCompWithVecs(Vec*,PetscInt,Vec,Vec*);
SLEPC_EXTERN PetscErrorCode VecCompGetSubVecs(Vec,PetscInt*,const Vec*[]);
SLEPC_EXTERN PetscErrorCode VecCompSetSubVecs(Vec,PetscInt,Vec*);

/* Some auxiliary functions */
SLEPC_EXTERN PetscErrorCode VecNormalizeComplex(Vec,Vec,PetscBool,PetscReal*);
SLEPC_EXTERN PetscErrorCode VecCheckOrthogonality(Vec[],PetscInt,Vec[],PetscInt,Mat,PetscViewer,PetscReal*);
SLEPC_EXTERN PetscErrorCode VecCheckOrthonormality(Vec[],PetscInt,Vec[],PetscInt,Mat,PetscViewer,PetscReal*);
SLEPC_EXTERN PetscErrorCode VecDuplicateEmpty(Vec,Vec*);
SLEPC_EXTERN PetscErrorCode VecSetRandomNormal(Vec,PetscRandom,Vec,Vec);

/* Deprecated functions */
PETSC_DEPRECATED_FUNCTION(3, 8, 0, "VecNormalizeComplex()", ) static inline PetscErrorCode SlepcVecNormalize(Vec xr,Vec xi,PetscBool c,PetscReal *nrm) {return VecNormalizeComplex(xr,xi,c,nrm);}
PETSC_DEPRECATED_FUNCTION(3, 8, 0, "VecCheckOrthogonality()", ) static inline PetscErrorCode SlepcCheckOrthogonality(Vec *V,PetscInt nv,Vec *W,PetscInt nw,Mat B,PetscViewer viewer,PetscReal *lev) {return VecCheckOrthogonality(V,nv,W,nw,B,viewer,lev);}
