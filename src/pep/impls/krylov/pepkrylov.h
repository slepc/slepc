/*
   Private header for TOAR and STOAR.

   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2014, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.

   SLEPc is free software: you can redistribute it and/or modify it under  the
   terms of version 3 of the GNU Lesser General Public License as published by
   the Free Software Foundation.

   SLEPc  is  distributed in the hope that it will be useful, but WITHOUT  ANY
   WARRANTY;  without even the implied warranty of MERCHANTABILITY or  FITNESS
   FOR  A  PARTICULAR PURPOSE. See the GNU Lesser General Public  License  for
   more details.

   You  should have received a copy of the GNU Lesser General  Public  License
   along with SLEPc. If not, see <http://www.gnu.org/licenses/>.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#if !defined(__TOAR_H)
#define __TOAR_H

typedef struct {
  PetscReal   keep;         /* restart parameter */
  PetscBool   lock;         /* locking/non-locking variant */
  PetscReal   dtol;         /* tolerance for deflation */
  PetscInt    d;            /* polynomial degree */
  PetscInt    ld;           /* leading dimension of auxiliary matrices */
  PetscScalar *S,*qB;       /* auxiliary matrices */
} PEP_TOAR;

#endif


PETSC_INTERN PetscErrorCode PEPExtractVectors_TOAR(PEP);
PETSC_INTERN PetscErrorCode PEPReset_TOAR(PEP);

