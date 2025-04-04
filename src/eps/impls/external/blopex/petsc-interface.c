/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* @@@ BLOPEX (version 1.1) LGPL Version 2.1 or above.See www.gnu.org. */
/* @@@ Copyright 2010 BLOPEX team https://github.com/lobpcg/blopex     */
/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* This code was developed by Merico Argentati, Andrew Knyazev, Ilya Lashuk and Evgueni Ovtchinnikov */

#include <slepcsys.h>
#include <petscblaslapack.h>
#include <interpreter.h>
#include <temp_multivector.h>
#include <fortran_matrix.h>
#include "petsc-interface.h"

static PetscRandom LOBPCG_RandomContext = NULL;

#if !defined(PETSC_USE_COMPLEX)
BlopexInt PETSC_dpotrf_interface (char *uplo,BlopexInt *n,double *a,BlopexInt * lda,BlopexInt *info)
{
  PetscBLASInt n_,lda_,info_;

  /* type conversion */
  n_ = *n;
  lda_ = *lda;
  info_ = *info;

  LAPACKpotrf_(uplo,&n_,(PetscScalar*)a,&lda_,&info_);

  *info = info_;
  return 0;
}

BlopexInt PETSC_dsygv_interface (BlopexInt *itype,char *jobz,char *uplo,BlopexInt *n,double *a,BlopexInt *lda,double *b,BlopexInt *ldb,double *w,double *work,BlopexInt *lwork,BlopexInt *info)
{
  PetscBLASInt itype_,n_,lda_,ldb_,lwork_,info_;

  itype_ = *itype;
  n_ = *n;
  lda_ = *lda;
  ldb_ = *ldb;
  lwork_ = *lwork;
  info_ = *info;

  LAPACKsygv_(&itype_,jobz,uplo,&n_,(PetscScalar*)a,&lda_,(PetscScalar*)b,&ldb_,(PetscScalar*)w,(PetscScalar*)work,&lwork_,&info_);

  *info = info_;
  return 0;
}
#else
BlopexInt PETSC_zpotrf_interface (char *uplo,BlopexInt *n,komplex *a,BlopexInt* lda,BlopexInt *info)
{
  PetscBLASInt n_,lda_,info_;

  /* type conversion */
  n_ = *n;
  lda_ = (PetscBLASInt)*lda;

  LAPACKpotrf_(uplo,&n_,(PetscScalar*)a,&lda_,&info_);

  *info = info_;
  return 0;
}

BlopexInt PETSC_zsygv_interface (BlopexInt *itype,char *jobz,char *uplo,BlopexInt *n,komplex *a,BlopexInt *lda,komplex *b,BlopexInt *ldb,double *w,komplex *work,BlopexInt *lwork,double *rwork,BlopexInt *info)
{
  PetscBLASInt itype_,n_,lda_,ldb_,lwork_,info_;

  itype_ = *itype;
  n_ = *n;
  lda_ = *lda;
  ldb_ = *ldb;
  lwork_ = *lwork;
  info_ = *info;

  LAPACKsygv_(&itype_,jobz,uplo,&n_,(PetscScalar*)a,&lda_,(PetscScalar*)b,&ldb_,(PetscReal*)w,(PetscScalar*)work,&lwork_,(PetscReal*)rwork,&info_);

  *info = info_;
  return 0;
}
#endif

static void *PETSC_MimicVector(void *vvector)
{
  Vec temp;

  PetscCallAbort(PETSC_COMM_SELF,VecDuplicate((Vec)vvector,&temp));
  return (void*)temp;
}

static BlopexInt PETSC_DestroyVector(void *vvector)
{
  Vec v = (Vec)vvector;

  PetscCall(VecDestroy(&v));
  return 0;
}

static BlopexInt PETSC_InnerProd(void *x,void *y,void *result)
{

  PetscCall(VecDot((Vec)x,(Vec)y,(PetscScalar*)result));
  return 0;
}

static BlopexInt PETSC_CopyVector(void *x,void *y)
{

  PetscCall(VecCopy((Vec)x,(Vec)y));
  return 0;
}

static BlopexInt PETSC_ClearVector(void *x)
{

  PetscCall(VecSet((Vec)x,0.0));
  return 0;
}

static BlopexInt PETSC_SetRandomValues(void* v,BlopexInt seed)
{

  /* note: without previous call to LOBPCG_InitRandomContext LOBPCG_RandomContext will be null,
    and VecSetRandom will use internal PetscRandom context */

  PetscCall(VecSetRandom((Vec)v,LOBPCG_RandomContext));
  return 0;
}

static BlopexInt PETSC_ScaleVector(double alpha,void *x)
{

  PetscCall(VecScale((Vec)x,alpha));
  return 0;
}

static BlopexInt PETSC_Axpy(void *alpha,void *x,void *y)
{

  PetscCall(VecAXPY((Vec)y,*(PetscScalar*)alpha,(Vec)x));
  return 0;
}

static BlopexInt PETSC_VectorSize(void *x)
{
  PetscInt N;
  (void)VecGetSize((Vec)x,&N);
  return N;
}

int LOBPCG_InitRandomContext(MPI_Comm comm,PetscRandom rand)
{
  /* PetscScalar rnd_bound = 1.0; */

  if (rand) {
    PetscCall(PetscObjectReference((PetscObject)rand));
    PetscCall(PetscRandomDestroy(&LOBPCG_RandomContext));
    LOBPCG_RandomContext = rand;
  } else PetscCall(PetscRandomCreate(comm,&LOBPCG_RandomContext));
  return 0;
}

int LOBPCG_SetFromOptionsRandomContext(void)
{
  PetscCall(PetscRandomSetFromOptions(LOBPCG_RandomContext));

#if defined(PETSC_USE_COMPLEX)
  PetscCall(PetscRandomSetInterval(LOBPCG_RandomContext,(PetscScalar)PetscCMPLX(-1.0,-1.0),(PetscScalar)PetscCMPLX(1.0,1.0)));
#else
  PetscCall(PetscRandomSetInterval(LOBPCG_RandomContext,(PetscScalar)-1.0,(PetscScalar)1.0));
#endif
  return 0;
}

int LOBPCG_DestroyRandomContext(void)
{

  PetscCall(PetscRandomDestroy(&LOBPCG_RandomContext));
  return 0;
}

int PETSCSetupInterpreter(mv_InterfaceInterpreter *i)
{
  i->CreateVector = PETSC_MimicVector;
  i->DestroyVector = PETSC_DestroyVector;
  i->InnerProd = PETSC_InnerProd;
  i->CopyVector = PETSC_CopyVector;
  i->ClearVector = PETSC_ClearVector;
  i->SetRandomValues = PETSC_SetRandomValues;
  i->ScaleVector = PETSC_ScaleVector;
  i->Axpy = PETSC_Axpy;
  i->VectorSize = PETSC_VectorSize;

  /* Multivector part */

  i->CreateMultiVector = mv_TempMultiVectorCreateFromSampleVector;
  i->CopyCreateMultiVector = mv_TempMultiVectorCreateCopy;
  i->DestroyMultiVector = mv_TempMultiVectorDestroy;

  i->Width = mv_TempMultiVectorWidth;
  i->Height = mv_TempMultiVectorHeight;
  i->SetMask = mv_TempMultiVectorSetMask;
  i->CopyMultiVector = mv_TempMultiVectorCopy;
  i->ClearMultiVector = mv_TempMultiVectorClear;
  i->SetRandomVectors = mv_TempMultiVectorSetRandom;
  i->Eval = mv_TempMultiVectorEval;

#if defined(PETSC_USE_COMPLEX)
  i->MultiInnerProd = mv_TempMultiVectorByMultiVector_complex;
  i->MultiInnerProdDiag = mv_TempMultiVectorByMultiVectorDiag_complex;
  i->MultiVecMat = mv_TempMultiVectorByMatrix_complex;
  i->MultiVecMatDiag = mv_TempMultiVectorByDiagonal_complex;
  i->MultiAxpy = mv_TempMultiVectorAxpy_complex;
  i->MultiXapy = mv_TempMultiVectorXapy_complex;
#else
  i->MultiInnerProd = mv_TempMultiVectorByMultiVector;
  i->MultiInnerProdDiag = mv_TempMultiVectorByMultiVectorDiag;
  i->MultiVecMat = mv_TempMultiVectorByMatrix;
  i->MultiVecMatDiag = mv_TempMultiVectorByDiagonal;
  i->MultiAxpy = mv_TempMultiVectorAxpy;
  i->MultiXapy = mv_TempMultiVectorXapy;
#endif

  return 0;
}
