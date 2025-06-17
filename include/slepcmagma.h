/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Macro definitions to use MAGMA functionality
*/

#pragma once

/* MANSEC = Sys */

#if defined(PETSC_HAVE_MAGMA)

#include <magma_v2.h>

SLEPC_EXTERN PetscErrorCode SlepcMagmaInit(void);

#define PetscCallMAGMA(func, ...) do { \
    magma_int_t magma_ierr_; \
    PetscStackPushExternal(PetscStringize(func)); \
    func(__VA_ARGS__,&magma_ierr_); \
    PetscStackPop; \
    PetscCheck(magma_ierr_ == MAGMA_SUCCESS,PETSC_COMM_SELF,PETSC_ERR_LIB,"Error calling %s: error code %d",PetscStringize(func(__VA_ARGS__,&magma_ierr)),(int)magma_ierr_); \
  } while (0)
#define CHKERRMAGMA(...) PetscCall(__VA_ARGS__)

#if defined(PETSC_USE_COMPLEX)
#if defined(PETSC_USE_REAL_SINGLE)
#define magma_xgeev(a,b,c,d,e,f,g,h,i,j,k,l,m,n) magma_cgeev((a),(b),(c),(magmaFloatComplex*)(d),(e),(magmaFloatComplex*)(f),(magmaFloatComplex*)(g),(h),(magmaFloatComplex*)(i),(j),(magmaFloatComplex*)(k),(l),(m),(n))
#define magma_xgesv_gpu(a,b,c,d,e,f,g,h)         magma_cgesv_gpu((a),(b),(magmaFloatComplex_ptr)(c),(d),(e),(magmaFloatComplex_ptr)(f),(g),(h))
#define magma_xgetrf_gpu(a,b,c,d,e,f)   magma_cgetrf_gpu((a),(b),(magmaFloatComplex_ptr)(c),(d),(e),(f))
#define magma_xgetri_gpu(a,b,c,d,e,f,g) magma_cgetri_gpu((a),(magmaFloatComplex_ptr)(b),(c),(d),(magmaFloatComplex_ptr)(e),(f),(g))
#define magma_get_xgetri_nb             magma_get_cgetri_nb
#else
#define magma_xgeev(a,b,c,d,e,f,g,h,i,j,k,l,m,n) magma_zgeev((a),(b),(c),(magmaDoubleComplex*)(d),(e),(magmaDoubleComplex*)(f),(magmaDoubleComplex*)(g),(h),(magmaDoubleComplex*)(i),(j),(magmaDoubleComplex*)(k),(l),(m),(n))
#define magma_xgesv_gpu(a,b,c,d,e,f,g,h)         magma_zgesv_gpu((a),(b),(magmaDoubleComplex_ptr)(c),(d),(e),(magmaDoubleComplex_ptr)(f),(g),(h))
#define magma_xgetrf_gpu(a,b,c,d,e,f)   magma_zgetrf_gpu((a),(b),(magmaDoubleComplex_ptr)(c),(d),(e),(f))
#define magma_xgetri_gpu(a,b,c,d,e,f,g) magma_zgetri_gpu((a),(magmaDoubleComplex_ptr)(b),(c),(d),(magmaDoubleComplex_ptr)(e),(f),(g))
#define magma_get_xgetri_nb             magma_get_zgetri_nb
#endif
#else
#if defined(PETSC_USE_REAL_SINGLE)
#define magma_xgeev                     magma_sgeev
#define magma_xgesv_gpu                 magma_sgesv_gpu
#define magma_xgetrf_gpu                magma_sgetrf_gpu
#define magma_xgetri_gpu                magma_sgetri_gpu
#define magma_get_xgetri_nb             magma_get_sgetri_nb
#else
#define magma_xgeev                     magma_dgeev
#define magma_xgesv_gpu                 magma_dgesv_gpu
#define magma_xgetrf_gpu                magma_dgetrf_gpu
#define magma_xgetri_gpu                magma_dgetri_gpu
#define magma_get_xgetri_nb             magma_get_dgetri_nb
#endif
#endif

#if defined(PETSC_USE_REAL_SINGLE)
#define magma_Cgesv_gpu(a,b,c,d,e,f,g,h)         magma_cgesv_gpu((a),(b),(magmaFloatComplex_ptr)(c),(d),(e),(magmaFloatComplex_ptr)(f),(g),(h))
#else
#define magma_Cgesv_gpu(a,b,c,d,e,f,g,h)         magma_zgesv_gpu((a),(b),(magmaDoubleComplex_ptr)(c),(d),(e),(magmaDoubleComplex_ptr)(f),(g),(h))
#endif

#endif
