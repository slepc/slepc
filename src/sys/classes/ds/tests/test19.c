/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test DSSortWithPermutation() on a NHEP.\n\n";

#include <slepcds.h>

int main(int argc,char **argv)
{
  DS             ds;
  SlepcSC        sc;
  PetscScalar    *A,*wr,*wi;
  PetscReal      re,im;
  PetscInt       i,j,n=12,*perm;

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
  PetscCall(PetscMalloc3(n,&wr,n,&wi,n,&perm));
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

  /* Reorder eigenvalues */
  for (i=0;i<n/2;i++) perm[i] = n/2+i;
  for (i=0;i<n/2;i++) perm[i+n/2] = i;
  PetscCall(DSSortWithPermutation(ds,perm,wr,wi));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Reordered eigenvalues =\n"));
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

  PetscCall(PetscFree3(wr,wi,perm));
  PetscCall(DSDestroy(&ds));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   test:
      suffix: 1
      filter: sed -e "s/[+-]\([0-9]\.[0-9]*i\)/+-\\1/" | sed -e "s/2.16612/2.16613/"

TEST*/
