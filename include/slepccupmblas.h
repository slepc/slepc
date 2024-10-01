/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Macro definitions to use cuBLAS and hipBLAS functionality
*/

#pragma once

#include <petscdevice.h>
#include <petsc/private/petsclegacycupmblas.h>

#if defined(PETSC_HAVE_CUDA)

/* complex single */
#if defined(PETSC_USE_COMPLEX)
#if defined(PETSC_USE_REAL_SINGLE)
#define cublasXgetrfBatched(a,b,c,d,e,f,g) cublasCgetrfBatched((a),(b),(cuComplex**)(c),(d),(e),(f),(g))
#define cublasXgetrsBatched(a,b,c,d,e,f,g,h,i,j,k) cublasCgetrsBatched((a),(b),(c),(d),(const cuComplex**)(e),(f),(g),(cuComplex**)(h),(i),(j),(k))
#else /* complex double */
#define cublasXgetrfBatched(a,b,c,d,e,f,g) cublasZgetrfBatched((a),(b),(cuDoubleComplex**)(c),(d),(e),(f),(g))
#define cublasXgetrsBatched(a,b,c,d,e,f,g,h,i,j,k) cublasZgetrsBatched((a),(b),(c),(d),(const cuDoubleComplex**)(e),(f),(g),(cuDoubleComplex**)(h),(i),(j),(k))
#endif
#else /* real single */
#if defined(PETSC_USE_REAL_SINGLE)
#define cublasXgetrfBatched cublasSgetrfBatched
#define cublasXgetrsBatched cublasSgetrsBatched
#else /* real double */
#define cublasXgetrfBatched cublasDgetrfBatched
#define cublasXgetrsBatched cublasDgetrsBatched
#endif
#endif

/* the following ones are used for PetscComplex in both real and complex scalars */
#if defined(PETSC_USE_REAL_SINGLE)
#define cublasXCaxpy(a,b,c,d,e,f,g)                cublasCaxpy((a),(b),(const cuComplex *)(c),(const cuComplex *)(d),(e),(cuComplex *)(f),(g))
#define cublasXCgemm(a,b,c,d,e,f,g,h,i,j,k,l,m,n)  cublasCgemm((a),(b),(c),(d),(e),(f),(const cuComplex *)(g),(const cuComplex *)(h),(i),(const cuComplex *)(j),(k),(const cuComplex *)(l),(cuComplex *)(m),(n))
#define cublasXCscal(a,b,c,d,e)                    cublasCscal((a),(b),(const cuComplex *)(c),(cuComplex *)(d),(e))
#else
#define cublasXCaxpy(a,b,c,d,e,f,g)                cublasZaxpy((a),(b),(const cuDoubleComplex *)(c),(const cuDoubleComplex *)(d),(e),(cuDoubleComplex *)(f),(g))
#define cublasXCgemm(a,b,c,d,e,f,g,h,i,j,k,l,m,n)  cublasZgemm((a),(b),(c),(d),(e),(f),(const cuDoubleComplex *)(g),(const cuDoubleComplex *)(h),(i),(const cuDoubleComplex *)(j),(k),(const cuDoubleComplex *)(l),(cuDoubleComplex *)(m),(n))
#define cublasXCscal(a,b,c,d,e)                    cublasZscal((a),(b),(const cuDoubleComplex *)(c),(cuDoubleComplex *)(d),(e))
#endif

#endif // PETSC_HAVE_CUDA

#if defined(PETSC_HAVE_HIP)

/* complex single */
#if defined(PETSC_USE_COMPLEX)
#if defined(PETSC_USE_REAL_SINGLE)
#else /* complex double */
#endif
#else /* real single */
#if defined(PETSC_USE_REAL_SINGLE)
#else /* real double */
#endif
#endif

#endif // PETSC_HAVE_HIP
