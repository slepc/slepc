/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Utility subroutines common to several impls
*/

#pragma once

#include <slepc/private/slepcimpl.h>

#if PetscDefined(HAVE_CUDA)
#include <slepccupmblas.h>

SLEPC_INTERN __host__ PetscErrorCode set_diagonal(PetscInt,PetscScalar*,PetscInt,PetscScalar);
__global__ void set_Cdiagonal_kernel(PetscInt,PetscComplex*,PetscInt,PetscReal,PetscReal,PetscInt);
SLEPC_INTERN __host__ PetscErrorCode set_Cdiagonal(PetscInt,PetscComplex*,PetscInt,PetscReal,PetscReal);
__global__ void shift_diagonal_kernel(PetscInt,PetscScalar*,PetscInt,PetscScalar,PetscInt);
SLEPC_INTERN __host__ PetscErrorCode shift_diagonal(PetscInt,PetscScalar*,PetscInt,PetscScalar);
__global__ void shift_Cdiagonal_kernel(PetscInt,PetscComplex*,PetscInt,PetscComplex,PetscInt);
SLEPC_INTERN __host__ PetscErrorCode shift_Cdiagonal(PetscInt,PetscComplex*,PetscInt,PetscReal,PetscReal);
__global__ void copy_array2D_S2C_kernel(PetscInt,PetscInt,PetscComplex*,PetscInt,PetscScalar*,PetscInt,PetscInt,PetscInt);
SLEPC_INTERN __host__ PetscErrorCode copy_array2D_S2C(PetscInt,PetscInt,PetscComplex*,PetscInt,PetscScalar*,PetscInt);
__global__ void copy_array2D_C2S_kernel(PetscInt,PetscInt,PetscScalar*,PetscInt,PetscComplex*,PetscInt,PetscInt,PetscInt);
SLEPC_INTERN __host__ PetscErrorCode copy_array2D_C2S(PetscInt,PetscInt,PetscScalar*,PetscInt,PetscComplex*,PetscInt);
__global__ void add_array2D_Conj_kernel(PetscInt,PetscInt,PetscComplex*,PetscInt,PetscInt,PetscInt);
SLEPC_INTERN __host__ PetscErrorCode add_array2D_Conj(PetscInt,PetscInt,PetscComplex*,PetscInt);
__global__ void getisreal_array2D_kernel(PetscInt,PetscInt,PetscComplex*,PetscInt,PetscBool*,PetscInt,PetscInt);
SLEPC_INTERN __host__ PetscErrorCode getisreal_array2D(PetscInt,PetscInt,PetscComplex*,PetscInt,PetscBool*);
__global__ void mult_diagonal_kernel(PetscInt,PetscScalar*,PetscInt,PetscScalar*,PetscInt);
SLEPC_INTERN __host__ PetscErrorCode mult_diagonal(PetscInt,PetscScalar*,PetscInt,PetscScalar*);

#endif /* PETSC_HAVE_CUDA */
