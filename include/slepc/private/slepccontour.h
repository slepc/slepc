/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#pragma once

#include <slepc/private/slepcimpl.h>
#include <petscksp.h>

/*
  CISS_BlockHankel - Builds a block Hankel matrix from the contents of Mu.
*/
static inline PetscErrorCode CISS_BlockHankel(PetscScalar *Mu,PetscInt s,PetscInt L,PetscInt M,PetscScalar *H)
{
  PetscInt i,j,k;

  PetscFunctionBegin;
  for (k=0;k<L*M;k++)
    for (j=0;j<M;j++)
      for (i=0;i<L;i++)
        H[j*L+i+k*L*M] = Mu[i+k*L+(j+s)*L*L];
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_SINGLE_LIBRARY_INTERN PetscErrorCode SlepcCISS_isGhost(Mat,PetscInt,PetscReal*,PetscReal,PetscBool*);
SLEPC_SINGLE_LIBRARY_INTERN PetscErrorCode SlepcCISS_BH_SVD(PetscScalar*,PetscInt,PetscReal,PetscReal*,PetscInt*);

/* Data structures and functions for contour integral methods (used in several classes) */
struct _n_SlepcContourData {
  PetscObject  parent;     /* parent object */
  PetscSubcomm subcomm;    /* subcommunicator for top level parallelization */
  PetscInt     npoints;    /* number of integration points assigned to the local subcomm */
  KSP          *ksp;       /* ksp array for storing factorizations at integration points */
  PetscInt     nksp;       /* number of ksp objects in array */
  Mat          *pA;        /* redundant copies of the matrices in the local subcomm */
  Mat          *pP;        /* redundant copies of the matrices (preconditioner) */
  PetscInt     nmat;       /* number of matrices in pA */
  Vec          xsub;       /* aux vector with parallel layout as redundant Mat */
  Vec          xdup;       /* aux vector with parallel layout as original Mat (with contiguous order) */
  VecScatter   scatterin;  /* to scatter from regular vector to xdup */
};
typedef struct _n_SlepcContourData* SlepcContourData;

SLEPC_SINGLE_LIBRARY_INTERN PetscErrorCode SlepcContourDataCreate(PetscInt,PetscInt,PetscObject,SlepcContourData*);
SLEPC_SINGLE_LIBRARY_INTERN PetscErrorCode SlepcContourDataReset(SlepcContourData);
SLEPC_SINGLE_LIBRARY_INTERN PetscErrorCode SlepcContourDataDestroy(SlepcContourData*);

SLEPC_SINGLE_LIBRARY_INTERN PetscErrorCode SlepcContourRedundantMat(SlepcContourData,PetscInt,Mat*,Mat*);
SLEPC_SINGLE_LIBRARY_INTERN PetscErrorCode SlepcContourScatterCreate(SlepcContourData,Vec);
