/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Necessary routines in BLAS and LAPACK not included in petscblaslapack.h
*/

#pragma once

#include <petscblaslapack.h>

/* MANSEC = Sys */

/* Macro to check nonzero info after LAPACK call */
#define SlepcCheckLapackInfo(routine,info) \
  do { \
    PetscCheck(!info,PETSC_COMM_SELF,PETSC_ERR_LIB,"Error in LAPACK subroutine %s: info=%" PetscBLASInt_FMT,routine,info); \
  } while (0)

/* LAPACK return type: we assume slange, etc. behave in the same way as snrm2 */
#if defined(PETSC_USE_REAL_SINGLE) && defined(PETSC_BLASLAPACK_SNRM2_RETURNS_DOUBLE)
#define SlepcLRT double
#else
#define SlepcLRT PetscReal
#endif

/* Special macro for srot, csrot, drot, zdrot (BLASMIXEDrot_) */
#if !defined(PETSC_USE_COMPLEX)
# define PETSC_BLASLAPACK_MIXEDPREFIX_ PETSC_BLASLAPACK_PREFIX_
#else
# if defined(PETSC_BLASLAPACK_CAPS)
#  if defined(PETSC_USE_REAL_SINGLE)
#   define PETSC_BLASLAPACK_MIXEDPREFIX_ CS
#  elif defined(PETSC_USE_REAL_DOUBLE)
#   define PETSC_BLASLAPACK_MIXEDPREFIX_ ZD
#  elif defined(PETSC_USE_REAL___FLOAT128)
#   define PETSC_BLASLAPACK_MIXEDPREFIX_ WQ
#  else
#   define PETSC_BLASLAPACK_MIXEDPREFIX_ KH
#  endif
# else
#  if defined(PETSC_USE_REAL_SINGLE)
#   define PETSC_BLASLAPACK_MIXEDPREFIX_ cs
#  elif defined(PETSC_USE_REAL_DOUBLE)
#   define PETSC_BLASLAPACK_MIXEDPREFIX_ zd
#  elif defined(PETSC_USE_REAL___FLOAT128)
#   define PETSC_BLASLAPACK_MIXEDPREFIX_ wq
#  else
#   define PETSC_BLASLAPACK_MIXEDPREFIX_ kh
#  endif
# endif
#endif
#if defined(PETSC_BLASLAPACK_CAPS)
#  define PETSCBLASMIXED(x,X) PETSC_PASTE3(PETSC_BLASLAPACK_MIXEDPREFIX_, X, PETSC_BLASLAPACK_SUFFIX_)
#else
#  define PETSCBLASMIXED(x,X) PETSC_PASTE3(PETSC_BLASLAPACK_MIXEDPREFIX_, x, PETSC_BLASLAPACK_SUFFIX_)
#endif

#include <slepcblaslapack_mangle.h>

/* LAPACK functions without string parameters */
BLAS_EXTERN void     BLASrot_(PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscReal*,PetscScalar*);
BLAS_EXTERN void     BLASMIXEDrot_(PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscReal*,PetscReal*);
#if !defined(SLEPC_MISSING_LAPACK_LAEV2)
BLAS_EXTERN void     LAPACKlaev2_(const PetscScalar*,const PetscScalar*,const PetscScalar*,PetscReal*,PetscReal*,PetscReal*,PetscScalar*);
#else
#define LAPACKlaev2_(a,b,c,d,e,f,g) PetscMissingLapack("LAEV2",a,b,c,d,e,f,g);
#endif
#if !defined(SLEPC_MISSING_LAPACK_GEHRD)
BLAS_EXTERN void     LAPACKgehrd_(const PetscBLASInt*,const PetscBLASInt*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,PetscScalar*,const PetscBLASInt*,PetscBLASInt*);
#else
#define LAPACKgehrd_(a,b,c,d,e,f,g,h,i) PetscMissingLapack("GEHRD",a,b,c,d,e,f,g,h,i);
#endif
#if !defined(SLEPC_MISSING_LAPACK_GEBRD)
BLAS_EXTERN void     LAPACKgebrd_(const PetscBLASInt*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscReal*,PetscReal*,PetscScalar*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
#else
#define LAPACKgebrd_(a,b,c,d,e,f,g,h,i,j,k) PetscMissingLapack("GEBRD",a,b,c,d,e,f,g,h,i,j,k);
#endif
#if !defined(SLEPC_MISSING_LAPACK_LARFG)
BLAS_EXTERN void     LAPACKlarfg_(const PetscBLASInt*,PetscScalar*,PetscScalar*,const PetscBLASInt*,PetscScalar*);
#else
#define LAPACKlarfg_(a,b,c,d,e) PetscMissingLapack("LARFG",a,b,c,d,e);
#endif
#if !defined(SLEPC_MISSING_LAPACK_LAG2)
BLAS_EXTERN void     LAPACKlag2_(const PetscReal*,const PetscBLASInt*,const PetscReal*,const PetscBLASInt*,const PetscReal*,PetscReal*,PetscReal*,PetscReal*,PetscReal*,PetscReal*);
#else
#define LAPACKlag2_(a,b,c,d,e,f,g,h,i,j) PetscMissingLapack("LAG2",a,b,c,d,e,f,g,h,i,j);
#endif
#if !defined(SLEPC_MISSING_LAPACK_LASV2)
BLAS_EXTERN void     LAPACKlasv2_(const PetscReal*,const PetscReal*,const PetscReal*,PetscReal*,PetscReal*,PetscReal*,PetscReal*,PetscReal*,PetscReal*);
#else
#define LAPACKlasv2_(a,b,c,d,e,f,g,h,i) PetscMissingLapack("LASV2",a,b,c,d,e,f,g,h,i);
#endif
#if !defined(SLEPC_MISSING_LAPACK_LARTG)
BLAS_EXTERN void     LAPACKlartg_(const PetscScalar*,const PetscScalar*,PetscReal*,PetscScalar*,PetscScalar*);
BLAS_EXTERN void     LAPACKREALlartg_(const PetscReal*,const PetscReal*,PetscReal*,PetscReal*,PetscReal*);
#else
#define LAPACKlartg_(a,b,c,d,e) PetscMissingLapack("LARTG",a,b,c,d,e);
#define LAPACKREALlartg_(a,b,c,d,e) PetscMissingLapack("LARTG",a,b,c,d,e);
#endif
#if !defined(SLEPC_MISSING_LAPACK_LAED4)
BLAS_EXTERN void     LAPACKlaed4_(const PetscBLASInt*,const PetscBLASInt*,const PetscReal*,const PetscReal*,PetscReal*,const PetscReal*,PetscReal*,PetscBLASInt*);
#else
#define LAPACKlaed4_(a,b,c,d,e,f,g,h) PetscMissingLapack("LAED4",a,b,c,d,e,f,g,h);
#endif
#if !defined(SLEPC_MISSING_LAPACK_LAMRG)
BLAS_EXTERN void     LAPACKlamrg_(const PetscBLASInt*,const PetscBLASInt*,const PetscReal*,const PetscBLASInt*,const PetscBLASInt*,PetscBLASInt*);
#else
#define LAPACKlamrg_(a,b,c,d,e,f) PetscMissingLapack("LAMRG",a,b,c,d,e,f);
#endif
#if !defined(SLEPC_MISSING_LAPACK_ORGHR)
BLAS_EXTERN void     LAPACKorghr_(const PetscBLASInt*,const PetscBLASInt*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,const PetscScalar*,PetscScalar*,const PetscBLASInt*,PetscBLASInt*);
#else
#define LAPACKorghr_(a,b,c,d,e,f,g,h,i) PetscMissingLapack("ORGHR",a,b,c,d,e,f,g,h,i);
#endif
#if !defined(PETSC_USE_COMPLEX)
#if !defined(SLEPC_MISSING_LAPACK_TGEXC)
BLAS_EXTERN void     LAPACKtgexc_(const PetscBLASInt*,const PetscBLASInt*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscBLASInt*);
#else
#define LAPACKtgexc_(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p) PetscMissingLapack("TGEXC",a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p);
#endif
BLAS_EXTERN void     LAPACKgeqp3_(const PetscBLASInt*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscScalar*,const PetscBLASInt*,PetscBLASInt*);
#else
#if !defined(SLEPC_MISSING_LAPACK_TGEXC)
BLAS_EXTERN void     LAPACKtgexc_(const PetscBLASInt*,const PetscBLASInt*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,const PetscBLASInt*,PetscBLASInt*);
#else
#define LAPACKtgexc_(a,b,c,d,e,f,g,h,i,j,k,l,m,n) PetscMissingLapack("TGEXC",a,b,c,d,e,f,g,h,i,j,k,l,m,n);
#endif
BLAS_EXTERN void     LAPACKgeqp3_(const PetscBLASInt*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscScalar*,const PetscBLASInt*,PetscReal*,PetscBLASInt*);
#endif

/* LAPACK functions with string parameters */

/* same name for real and complex */
BLAS_EXTERN void     BLAStrmm_(const char*,const char*,const char*,const char*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*);
BLAS_EXTERN SlepcLRT LAPACKlange_(const char*,const PetscBLASInt*,const PetscBLASInt*,const PetscScalar*,const PetscBLASInt*,PetscReal*);
#if !defined(SLEPC_MISSING_LAPACK_LANHS)
BLAS_EXTERN SlepcLRT LAPACKlanhs_(const char*,const PetscBLASInt*,const PetscScalar*,const PetscBLASInt*,PetscReal*);
#else
static inline SlepcLRT LAPACKlanhs_(const char *norm,const PetscBLASInt *n,const PetscScalar *A,const PetscBLASInt *lda,PetscReal *work) {return LAPACKlange_(norm,n,n,A,lda,work);}
#endif
#if !defined(SLEPC_MISSING_LAPACK_LARF)
BLAS_EXTERN void     LAPACKlarf_(const char*,const PetscBLASInt*,const PetscBLASInt*,const PetscScalar*,const PetscBLASInt*,const PetscScalar*,PetscScalar*,const PetscBLASInt*,PetscScalar*);
#else
#define LAPACKlarf_(a,b,c,d,e,f,g,h,i) PetscMissingLapack("LARF",a,b,c,d,e,f,g,h,i);
#endif
BLAS_EXTERN SlepcLRT LAPACKlansy_(const char*,const char*,const PetscBLASInt*,const PetscScalar*,const PetscBLASInt*,PetscReal*);
#if !defined(SLEPC_MISSING_LAPACK_TRSYL)
BLAS_EXTERN void     LAPACKtrsyl_(const char*,const char*,const PetscBLASInt*,const PetscBLASInt*,const PetscBLASInt*,const PetscScalar*,const PetscBLASInt*,const PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscReal*,PetscBLASInt*);
#else
#define LAPACKtrsyl_(a,b,c,d,e,f,g,h,i,j,k,l,m) PetscMissingLapack("TRSYL",a,b,c,d,e,f,g,h,i,j,k,l,m);
#endif
#if !defined(SLEPC_MISSING_LAPACK_BDSQR)
BLAS_EXTERN void     LAPACKbdsqr_(const char*,const PetscBLASInt*,const PetscBLASInt*,const PetscBLASInt*,const PetscBLASInt*,PetscReal*,PetscReal*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscReal*,PetscBLASInt*);
#else
#define LAPACKbdsqr_(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o) PetscMissingLapack("BDSQR",a,b,c,d,e,f,g,h,i,j,k,l,m,n,o);
#endif

/* subroutines in which we use only the real version, do not care whether they have different name */
#if !defined(SLEPC_MISSING_LAPACK_STEVR)
BLAS_EXTERN void     LAPACKstevr_(const char*,const char*,PetscBLASInt*,PetscReal*,PetscReal*,PetscReal*,PetscReal*,PetscBLASInt*,PetscBLASInt*,PetscReal*,PetscBLASInt*,PetscReal*,PetscReal*,PetscBLASInt*,PetscBLASInt*,PetscReal*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*);
#else
#define LAPACKstevr_(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t) PetscMissingLapack("STEVR",a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t);
#endif
#if !defined(SLEPC_MISSING_LAPACK_BDSDC)
BLAS_EXTERN void     LAPACKbdsdc_(const char*,const char*,PetscBLASInt*,PetscReal*,PetscReal*,PetscReal*,PetscBLASInt*,PetscReal*,PetscBLASInt*,PetscReal*,PetscBLASInt*,PetscReal*,PetscBLASInt*,PetscBLASInt*);
#else
#define LAPACKbdsdc_(a,b,c,d,e,f,g,h,i,j,k,l,m,n) PetscMissingLapack("BDSDC",a,b,c,d,e,f,g,h,i,j,k,l,m,n);
#endif
BLAS_EXTERN SlepcLRT LAPACKlamch_(const char*);
BLAS_EXTERN SlepcLRT LAPACKlamc3_(PetscReal*,PetscReal*);

/* subroutines with different name in real/complex */
#if !defined(SLEPC_MISSING_LAPACK_ORGTR)
BLAS_EXTERN void     LAPACKorgtr_(const char*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
#else
#define LAPACKorgtr_(a,b,c,d,e,f,g,h) PetscMissingLapack("ORGTR",a,b,c,d,e,f,g,h);
#endif
#if !defined(SLEPC_MISSING_LAPACK_ORMBR)
BLAS_EXTERN void     LAPACKormbr_(const char*,const char*,const char*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
#else
#define LAPACKormbr_(a,b,c,d,e,f,g,h,i,j,k,l,m,n) PetscMissingLapack("ORMBR",a,b,c,d,e,f,g,h,i,j,k,l,m,n);
#endif
#if !defined(SLEPC_MISSING_LAPACK_SYTRD)
BLAS_EXTERN void     LAPACKsytrd_(const char*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscReal*,PetscReal*,PetscScalar*,PetscScalar*,const PetscBLASInt*,PetscBLASInt*);
#else
#define LAPACKsytrd_(a,b,c,d,e,f,g,h,i,j) PetscMissingLapack("SYTRD",a,b,c,d,e,f,g,h,i,j);
#endif
#if !defined(PETSC_USE_COMPLEX)
BLAS_EXTERN void     LAPACKsyevd_(const char*,const char*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*);
BLAS_EXTERN void     LAPACKsygvd_(PetscBLASInt*,const char*,const char*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*);
#else
BLAS_EXTERN void     LAPACKsyevd_(const char*,const char*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscReal*,PetscScalar*,PetscBLASInt*,PetscReal*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*);
BLAS_EXTERN void     LAPACKsygvd_(PetscBLASInt*,const char*,const char*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscReal*,PetscScalar*,PetscBLASInt*,PetscReal*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*);
#endif

/* subroutines with different signature in real/complex */
#if !defined(PETSC_USE_COMPLEX)
BLAS_EXTERN void     LAPACKggev_(const char*,const char*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,PetscScalar*,PetscScalar*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscBLASInt*);
#if !defined(SLEPC_MISSING_LAPACK_GGEV3)
BLAS_EXTERN void     LAPACKggev3_(const char*,const char*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,PetscScalar*,PetscScalar*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscBLASInt*);
#else
#define LAPACKggev3_(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q) PetscMissingLapack("GGEV3",a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q);
#endif
#if !defined(SLEPC_MISSING_LAPACK_GGES3)
BLAS_EXTERN void     LAPACKgges3_(const char*,const char*,const char*,PetscBLASInt(*)(void),const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscScalar*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscBLASInt*,PetscBLASInt*);
#else
#define LAPACKgges3_(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u) PetscMissingLapack("GGES3",a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u)
#endif
#if !defined(SLEPC_MISSING_LAPACK_GGSVD)
BLAS_EXTERN void     LAPACKggsvd_(const char*,const char*,const char*,const PetscBLASInt*,const PetscBLASInt*,const PetscBLASInt*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscReal*,PetscReal*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
#else
#define LAPACKggsvd_(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w) PetscMissingLapack("GGSVD",a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w);
#endif
#if !defined(SLEPC_MISSING_LAPACK_GGSVD3)
BLAS_EXTERN void     LAPACKggsvd3_(const char*,const char*,const char*,const PetscBLASInt*,const PetscBLASInt*,const PetscBLASInt*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscReal*,PetscReal*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscBLASInt*,PetscBLASInt*);
#else
#define LAPACKggsvd3_(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x) PetscMissingLapack("GGSVD3",a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x);
#endif
#if !defined(SLEPC_MISSING_LAPACK_TREVC)
BLAS_EXTERN void     LAPACKtrevc_(const char*,const char*,PetscBLASInt*,const PetscBLASInt*,const PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,const PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*);
#else
#define LAPACKtrevc_(a,b,c,d,e,f,g,h,i,j,k,l,m,n) PetscMissingLapack("TREVC",a,b,c,d,e,f,g,h,i,j,k,l,m,n);
#endif
BLAS_EXTERN void     LAPACKgeevx_(const char*,const char*,const char*,const char*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,PetscScalar*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscScalar*,PetscScalar*,PetscScalar*,const PetscBLASInt*,PetscBLASInt*,PetscBLASInt*);
BLAS_EXTERN void     LAPACKgees_(const char*,const char*,PetscBLASInt(*)(PetscReal,PetscReal),const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscBLASInt*,PetscBLASInt*);
#if !defined(SLEPC_MISSING_LAPACK_TREXC)
BLAS_EXTERN void     LAPACKtrexc_(const char*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*);
#else
#define LAPACKtrexc_(a,b,c,d,e,f,g,h,i,j) PetscMissingLapack("TREXC",a,b,c,d,e,f,g,h,i,j);
#endif
BLAS_EXTERN void     LAPACKgesdd_(const char*,const PetscBLASInt*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscReal*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscBLASInt*,PetscBLASInt*);
#if !defined(SLEPC_MISSING_LAPACK_TGEVC)
BLAS_EXTERN void     LAPACKtgevc_(const char*,const char*,const PetscBLASInt*,const PetscBLASInt*,const PetscScalar*,const PetscBLASInt*,const PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,const PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*);
#else
#define LAPACKtgevc_(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p) PetscMissingLapack("TGEVC",a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p);
#endif
#if !defined(SLEPC_MISSING_LAPACK_HSEIN)
BLAS_EXTERN void     LAPACKhsein_(const char*,const char*,const char*,PetscBLASInt*,const PetscBLASInt*,const PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscScalar*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,const PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*);
#else
#define LAPACKhsein_(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s) PetscMissingLapack("HSEIN",a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s);
#endif
#if !defined(SLEPC_MISSING_LAPACK_STEDC)
BLAS_EXTERN void     LAPACKstedc_(const char*,const PetscBLASInt*,PetscReal*,PetscReal*,PetscScalar*,const PetscBLASInt*,PetscReal*,const PetscBLASInt*,PetscBLASInt*,const PetscBLASInt*,PetscBLASInt*);
#else
#define LAPACKstedc_(a,b,c,d,e,f,g,h,i,j,k) PetscMissingLapack("STEDC",a,b,c,d,e,f,g,h,i,j,k);
#endif
#if !defined(SLEPC_MISSING_LAPACK_LASCL)
BLAS_EXTERN void     LAPACKlascl_(const char*,const PetscBLASInt*,const PetscBLASInt*,const PetscScalar*,const PetscScalar*,const PetscBLASInt*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscBLASInt*);
#else
#define LAPACKlascl_(a,b,c,d,e,f,g,h,i,j) PetscMissingLapack("LASCL",a,b,c,d,e,f,g,h,i,j);
#endif
#else
BLAS_EXTERN void     LAPACKggev_(const char*,const char*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,PetscScalar*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscReal*,PetscBLASInt*);
#if !defined(SLEPC_MISSING_LAPACK_GGEV3)
BLAS_EXTERN void     LAPACKggev3_(const char*,const char*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,PetscScalar*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscReal*,PetscBLASInt*);
#else
#define LAPACKggev3_(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q) PetscMissingLapack("GGEV3",a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q);
#endif
#if !defined(SLEPC_MISSING_LAPACK_GGES3)
BLAS_EXTERN void     LAPACKgges3_(const char*,const char*,const char*,PetscBLASInt(*)(void),const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscReal*,PetscBLASInt*,PetscBLASInt*);
#else
#define LAPACKgges3_(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u) PetscMissingLapack("GGES3",a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u)
#endif
#if !defined(SLEPC_MISSING_LAPACK_GGSVD)
BLAS_EXTERN void     LAPACKggsvd_(const char*,const char*,const char*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscReal*,PetscReal*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscReal*,PetscBLASInt*,PetscBLASInt*);
#else
#define LAPACKggsvd_(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x) PetscMissingLapack("GGSVD",a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x);
#endif
#if !defined(SLEPC_MISSING_LAPACK_GGSVD3)
BLAS_EXTERN void     LAPACKggsvd3_(const char*,const char*,const char*,const PetscBLASInt*,const PetscBLASInt*,const PetscBLASInt*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscReal*,PetscReal*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscReal*,PetscBLASInt*,PetscBLASInt*);
#else
#define LAPACKggsvd3_(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y) PetscMissingLapack("GGSVD3",a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y);
#endif
#if !defined(SLEPC_MISSING_LAPACK_TREVC)
BLAS_EXTERN void     LAPACKtrevc_(const char*,const char*,const PetscBLASInt*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,const PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscReal*,PetscBLASInt*);
#else
#define LAPACKtrevc_(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o) PetscMissingLapack("TREVC",a,b,c,d,e,f,g,h,i,j,k,l,m,n,o);
#endif
BLAS_EXTERN void     LAPACKgeevx_(const char*,const char*,const char*,const char*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscBLASInt*,PetscBLASInt*,PetscReal*,PetscReal*,PetscReal*,PetscReal*,PetscScalar*,const PetscBLASInt*,PetscReal*,PetscBLASInt*);
BLAS_EXTERN void     LAPACKgees_(const char*,const char*,PetscBLASInt(*)(PetscScalar),const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscReal*,PetscBLASInt*,PetscBLASInt*);
#if !defined(SLEPC_MISSING_LAPACK_TREXC)
BLAS_EXTERN void     LAPACKtrexc_(const char*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,const PetscBLASInt*,const PetscBLASInt*,PetscBLASInt*);
#else
#define LAPACKtrexc_(a,b,c,d,e,f,g,h,i) PetscMissingLapack("TREXC",a,b,c,d,e,f,g,h,i);
#endif
BLAS_EXTERN void     LAPACKgesdd_(const char*,const PetscBLASInt*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscReal*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscReal*,PetscBLASInt*,PetscBLASInt*);
#if !defined(SLEPC_MISSING_LAPACK_TGEVC)
BLAS_EXTERN void     LAPACKtgevc_(const char*,const char*,const PetscBLASInt*,const PetscBLASInt*,const PetscScalar*,const PetscBLASInt*,const PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,const PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscReal*,PetscBLASInt*);
#else
#define LAPACKtgevc_(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q) PetscMissingLapack("TGEVC",a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q);
#endif
#if !defined(SLEPC_MISSING_LAPACK_HSEIN)
BLAS_EXTERN void     LAPACKhsein_(const char*,const char*,const char*,const PetscBLASInt*,const PetscBLASInt*,const PetscScalar*,const PetscBLASInt*,PetscScalar*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,const PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscReal*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*);
#else
#define LAPACKhsein_(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s) PetscMissingLapack("HSEIN",a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s);
#endif
#if !defined(SLEPC_MISSING_LAPACK_STEDC)
BLAS_EXTERN void     LAPACKstedc_(const char*,const PetscBLASInt*,PetscReal*,PetscReal*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscReal*,const PetscBLASInt*,PetscBLASInt*,const PetscBLASInt*,PetscBLASInt*);
#else
#define LAPACKstedc_(a,b,c,d,e,f,g,h,i,j,k,l,m) PetscMissingLapack("STEDC",a,b,c,d,e,f,g,h,i,j,k,l,m);
#endif
#if !defined(SLEPC_MISSING_LAPACK_LASCL)
BLAS_EXTERN void     LAPACKlascl_(const char*,const PetscBLASInt*,const PetscBLASInt*,const PetscReal*,const PetscReal*,const PetscBLASInt*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscBLASInt*);
#else
#define LAPACKlascl_(a,b,c,d,e,f,g,h,i,j) PetscMissingLapack("LASCL",a,b,c,d,e,f,g,h,i,j);
#endif
#endif

#if defined(PETSC_HAVE_COMPLEX)
/* complex subroutines to be called with scalar-type=real */
BLAS_EXTERN void BLASCOMPLEXgemm_(const char*,const char*,const PetscBLASInt*,const PetscBLASInt*,const PetscBLASInt*,const PetscComplex*,const PetscComplex*,const PetscBLASInt*,const PetscComplex*,const PetscBLASInt*,const PetscComplex*,PetscComplex*,const PetscBLASInt*);
BLAS_EXTERN void BLASCOMPLEXscal_(const PetscBLASInt*,const PetscComplex*,PetscComplex*,const PetscBLASInt*);
BLAS_EXTERN void LAPACKCOMPLEXgesv_(const PetscBLASInt*,const PetscBLASInt*,PetscComplex*,const PetscBLASInt*,PetscBLASInt*,PetscComplex*,const PetscBLASInt*,PetscBLASInt*);
#endif
