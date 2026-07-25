/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test DSSynchronize() on a NHEP.\n\n";

#include <slepcds.h>

PetscErrorCode CheckArray(PetscScalar *A,const char *label,PetscInt k)
{
  PetscInt       j;
  PetscMPIInt    p,size,rank,k_;
  PetscScalar    dif,*buf;
  PetscReal      error;

  PetscFunctionBeginUser;
  PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD,&size));
  PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD,&rank));
  PetscCall(PetscMPIIntCast(k,&k_));
  if (rank) PetscCallMPI(MPI_Send(A,k_,MPIU_SCALAR,0,111,PETSC_COMM_WORLD));
  else {
    PetscCall(PetscMalloc1(k,&buf));
    for (p=1;p<size;p++) {
      PetscCallMPI(MPI_Recv(buf,k_,MPIU_SCALAR,p,111,PETSC_COMM_WORLD,MPI_STATUS_IGNORE));
      dif = 0.0;
      for (j=0;j<k;j++) dif += A[j]-buf[j];
      error = PetscAbsScalar(dif);
      if (error>10*PETSC_MACHINE_EPSILON) PetscCall(PetscPrintf(PETSC_COMM_WORLD," Array %s differs in proc %d: %g\n",label,(int)p,(double)error));
    }
    PetscCall(PetscFree(buf));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

int main(int argc,char **argv)
{
  DS             ds;
  SlepcSC        sc;
  PetscScalar    *A,*Q,*wr,*wi;
  PetscReal      re,im;
  PetscInt       i,j,n=10;
  PetscMPIInt    size;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Solve a Dense System of type NHEP - dimension %" PetscInt_FMT ".\n",n));

  /* Create DS object */
  PetscCall(DSCreate(PETSC_COMM_WORLD,&ds));
  PetscCall(DSSetType(ds,DSNHEP));
  PetscCall(DSSetFromOptions(ds));
  PetscCall(DSAllocate(ds,n));
  PetscCall(DSSetDimensions(ds,n,0,0));

  /* Fill with Grcar matrix */
  PetscCall(DSGetArray(ds,DS_MAT_A,&A));
  for (i=1;i<n;i++) A[i+(i-1)*n]=-1.0;
  for (j=0;j<4;j++) {
    for (i=0;i<n-j;i++) A[i+(i+j)*n]=1.0;
  }
  PetscCall(DSRestoreArray(ds,DS_MAT_A,&A));
  PetscCall(DSSetState(ds,DS_STATE_INTERMEDIATE));

  /* Solve */
  PetscCall(PetscMalloc2(n,&wr,n,&wi));
  PetscCall(DSGetSlepcSC(ds,&sc));
  sc->comparison    = SlepcCompareLargestMagnitude;
  sc->comparisonctx = NULL;
  sc->map           = NULL;
  sc->mapobj        = NULL;
  PetscCall(DSSolve(ds,wr,wi));
  PetscCall(DSSort(ds,wr,wi,NULL,NULL,NULL));

  /* Print eigenvalues */
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Computed eigenvalues =\n"));
  for (i=0;i<n;i++) {
#if PetscDefined(USE_COMPLEX)
    re = PetscRealPart(wr[i]);
    im = PetscImaginaryPart(wr[i]);
#else
    re = wr[i];
    im = wi[i];
#endif
    if (PetscAbs(im)<1e-10) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"  %.5f\n",(double)re));
    else PetscCall(PetscPrintf(PETSC_COMM_WORLD,"  %.5f%+.5fi\n",(double)re,(double)im));
  }

  /* Synchronize data and check */
  PetscCall(DSSynchronize(ds,wr,wi));
  PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD,&size));
  if (size>1) {
    PetscCall(CheckArray(wr,"wr",n));
#if !PetscDefined(USE_COMPLEX)
    PetscCall(CheckArray(wi,"wi",n));
#endif
    PetscCall(DSGetArray(ds,DS_MAT_A,&A));
    PetscCall(CheckArray(A,"A",n*n));
    PetscCall(DSRestoreArray(ds,DS_MAT_A,&A));
    PetscCall(DSGetArray(ds,DS_MAT_Q,&Q));
    PetscCall(CheckArray(Q,"Q",n*n));
    PetscCall(DSRestoreArray(ds,DS_MAT_Q,&Q));
  }

  PetscCall(PetscFree2(wr,wi));
  PetscCall(DSDestroy(&ds));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   testset:
      nsize: {{1 2 3}}
      filter: sed -e "s/[+-]\([0-9]\.[0-9]*i\)/+-\\1/" | sed -e "s/1.58254/1.58255/" | sed -e "s/1.75989/1.75988/" | sed -e "s/1.01945/1.01946/"
      output_file: output/test18_1.out
      test:
         suffix: 1
         args: -ds_parallel redundant
      test:
         suffix: 2
         args: -ds_parallel synchronized

TEST*/
