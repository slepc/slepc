/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#pragma once

#include <slepclme.h>
#include <slepc/private/slepcimpl.h>

/* SUBMANSEC = LME */

SLEPC_EXTERN PetscBool LMERegisterAllCalled;
SLEPC_EXTERN PetscBool LMEMonitorRegisterAllCalled;
SLEPC_EXTERN PetscErrorCode LMERegisterAll(void);
SLEPC_EXTERN PetscErrorCode LMEMonitorRegisterAll(void);
SLEPC_EXTERN PetscLogEvent LME_SetUp,LME_Solve,LME_ComputeError;

typedef struct _LMEOps *LMEOps;

struct _LMEOps {
  PetscErrorCode (*solve[sizeof(LMEProblemType)])(LME);
  PetscErrorCode (*setup)(LME);
  PetscErrorCode (*setfromoptions)(LME,PetscOptionItems);
  PetscErrorCode (*publishoptions)(LME);
  PetscErrorCode (*destroy)(LME);
  PetscErrorCode (*reset)(LME);
  PetscErrorCode (*view)(LME,PetscViewer);
};

/*
     Maximum number of monitors you can run with a single LME
*/
#define MAXLMEMONITORS 5

/*
   Defines the LME data structure.
*/
struct _p_LME {
  PETSCHEADER(struct _LMEOps);
  /*------------------------- User parameters ---------------------------*/
  Mat            A,B,D,E;        /* the coefficient matrices */
  Mat            C;              /* the right-hand side */
  Mat            X;              /* the solution */
  LMEProblemType problem_type;   /* which kind of equation to be solved */
  PetscInt       max_it;         /* maximum number of iterations */
  PetscInt       ncv;            /* number of basis vectors */
  PetscReal      tol;            /* tolerance */
  PetscBool      errorifnotconverged;    /* error out if LMESolve() does not converge */

  /*-------------- User-provided functions and contexts -----------------*/
  PetscErrorCode    (*monitor[MAXLMEMONITORS])(LME,PetscInt,PetscReal,void*);
  PetscCtxDestroyFn *monitordestroy[MAXLMEMONITORS];
  void              *monitorcontext[MAXLMEMONITORS];
  PetscInt          numbermonitors;

  /*----------------- Child objects and working data -------------------*/
  BV             V;              /* set of basis vectors */
  PetscInt       nwork;          /* number of work vectors */
  Vec            *work;          /* work vectors */
  void           *data;          /* placeholder for solver-specific stuff */

  /* ----------------------- Status variables -------------------------- */
  PetscInt       its;            /* number of iterations so far computed */
  PetscReal      errest;         /* error estimate */
  PetscInt       setupcalled;
  LMEConvergedReason reason;
};

SLEPC_INTERN PetscErrorCode LMEDenseRankSVD(LME,PetscInt,PetscScalar*,PetscInt,PetscScalar*,PetscInt,PetscInt*);

/*
    Macros to test valid LME arguments
*/
#if !defined(PETSC_USE_DEBUG)

#define LMECheckCoeff(h,A,mat,eq) do {(void)(h);} while (0)

#else

#define LMECheckCoeff(h,A,mat,eq) \
  do { \
    PetscCheck(A,PetscObjectComm((PetscObject)(h)),PETSC_ERR_ARG_WRONGSTATE,"%s matrix equation requires coefficient matrix %s",eq,mat); \
  } while (0)

#endif

/* functions interfaced from Fortran library SLICOT */
#if defined(SLEPC_HAVE_SLICOT)

#if defined(SLEPC_SLICOT_HAVE_UNDERSCORE)
#define SLEPC_SLICOT(lcase,ucase) lcase##_
#elif defined(SLEPC_SLICOT_HAVE_CAPS)
#define SLEPC_SLICOT(lcase,ucase) ucase
#else
#define SLEPC_SLICOT(lcase,ucase) lcase
#endif

#define SLICOTsb03od_(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q) SLEPC_SLICOT(sb03od,SB03OD) ((a),(b),(c),(d),(e),(f),(g),(h),(i),(j),(k),(l),(m),(n),(o),(p),(q),1,1,1)
SLEPC_EXTERN void SLEPC_SLICOT(sb03od,SB03OD)(const char*,const char*,const char*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscReal*,PetscScalar*,PetscScalar*,PetscReal*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt,PetscBLASInt,PetscBLASInt);
#define SLICOTsb03md_(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t) SLEPC_SLICOT(sb03md,SB03MD) ((a),(b),(c),(d),(e),(f),(g),(h),(i),(j),(k),(l),(m),(n),(o),(p),(q),(r),(s),(t),1,1,1,1)
SLEPC_EXTERN void SLEPC_SLICOT(sb03md,SB03MD)(const char*,const char*,const char*,const char*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscReal*,PetscReal*,PetscReal*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscReal*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt,PetscBLASInt,PetscBLASInt,PetscBLASInt);

#endif
