/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   User interface for various matrix operations added in SLEPc
*/

#pragma once

#include <petscmat.h>

/* SUBMANSEC = Sys */

SLEPC_EXTERN PetscErrorCode MatCreateTile(PetscScalar,Mat,PetscScalar,Mat,PetscScalar,Mat,PetscScalar,Mat,Mat*);
SLEPC_EXTERN PetscErrorCode MatCreateVecsEmpty(Mat,Vec*,Vec*);
SLEPC_EXTERN PetscErrorCode MatNormEstimate(Mat,Vec,Vec,PetscReal*);

/* Matrices for structured eigenproblems */
SLEPC_EXTERN PetscErrorCode MatCreateBSE(Mat,Mat,Mat*);
SLEPC_EXTERN PetscErrorCode MatCreateHamiltonian(Mat,Mat,Mat,Mat*);

/* Deprecated functions */
PETSC_DEPRECATED_FUNCTION(3, 6, 0, "MatCreateRedundantMatrix() followed by MatConvert()", ) static inline PetscErrorCode SlepcMatConvertSeqDense(Mat mat,Mat *newmat)
{
  Mat Ar;

  PetscFunctionBegin;
  PetscCall(MatCreateRedundantMatrix(mat,0,PETSC_COMM_SELF,MAT_INITIAL_MATRIX,&Ar));
  PetscCall(MatConvert(Ar,MATSEQDENSE,MAT_INITIAL_MATRIX,newmat));
  PetscCall(MatDestroy(&Ar));
  PetscFunctionReturn(PETSC_SUCCESS);
}
PETSC_DEPRECATED_FUNCTION(3, 8, 0, "MatCreateTile()", ) static inline PetscErrorCode SlepcMatTile(PetscScalar a,Mat A,PetscScalar b,Mat B,PetscScalar c,Mat C,PetscScalar d,Mat D,Mat *G) {return MatCreateTile(a,A,b,B,c,C,d,D,G);}
