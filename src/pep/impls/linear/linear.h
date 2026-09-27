/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Private header for PEPLINEAR
*/

#pragma once

typedef struct {
  PetscBool  explicitmatrix;
  PEP        pep;
  PetscReal  sfactor,dsfactor; /* scaling factors */
  Mat        A,B;              /* matrices of generalized eigenproblem */
  EPS        eps;              /* linear eigensolver for Az=lBz */
  PetscBool  usereps;          /* eps provided by user */
  Mat        M,C,K;            /* copy of PEP coefficient matrices */
  Vec        w[6];             /* work vectors */
  PetscReal  alpha,beta;       /* coefficients defining the linearization */
} PEP_LINEAR;

/* General case for implicit matrices of degree d */
SLEPC_INTERN PetscErrorCode MatMult_Linear(Mat,Vec,Vec);

/* N */
SLEPC_INTERN PetscErrorCode MatCreateExplicit_Linear_NA(MPI_Comm,PEP_LINEAR*,Mat*);
SLEPC_INTERN PetscErrorCode MatCreateExplicit_Linear_NB(MPI_Comm,PEP_LINEAR*,Mat*);

/* S */
SLEPC_INTERN PetscErrorCode MatCreateExplicit_Linear_SA(MPI_Comm,PEP_LINEAR*,Mat*);
SLEPC_INTERN PetscErrorCode MatCreateExplicit_Linear_SB(MPI_Comm,PEP_LINEAR*,Mat*);

/* H */
SLEPC_INTERN PetscErrorCode MatCreateExplicit_Linear_HA(MPI_Comm,PEP_LINEAR*,Mat*);
SLEPC_INTERN PetscErrorCode MatCreateExplicit_Linear_HB(MPI_Comm,PEP_LINEAR*,Mat*);
