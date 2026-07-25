/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test DSGHIEP with compact storage.\n\n";

#include <slepcds.h>

int main(int argc,char **argv)
{
  DS             ds;
  SlepcSC        sc;
  PetscReal      *T,*s,re,im;
  PetscScalar    *eigr,*eigi;
  PetscInt       i,n=10,l=2,k=5,ld;
  PetscViewer    viewer;
  PetscBool      verbose;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Solve a Dense System of type GHIEP with compact storage - dimension %" PetscInt_FMT ".\n",n));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-l",&l,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-k",&k,NULL));
  PetscCheck(l<=n && k<=n && l<=k,PETSC_COMM_WORLD,PETSC_ERR_USER_INPUT,"Wrong value of dimensions");
  PetscCall(PetscOptionsHasName(NULL,NULL,"-verbose",&verbose));

  /* Create DS object */
  PetscCall(DSCreate(PETSC_COMM_WORLD,&ds));
  PetscCall(DSSetType(ds,DSGHIEP));
  PetscCall(DSSetFromOptions(ds));
  ld = n+2;  /* test leading dimension larger than n */
  PetscCall(DSAllocate(ds,ld));
  PetscCall(DSSetDimensions(ds,n,l,k));
  PetscCall(DSSetCompact(ds,PETSC_TRUE));

  /* Set up viewer */
  PetscCall(PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&viewer));
  PetscCall(PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_INFO_DETAIL));
  PetscCall(DSView(ds,viewer));
  PetscCall(PetscViewerPopFormat(viewer));
  PetscCall(PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB));

  /* Fill arrow-tridiagonal matrix */
  PetscCall(DSGetArrayReal(ds,DS_MAT_T,&T));
  PetscCall(DSGetArrayReal(ds,DS_MAT_D,&s));
  for (i=0;i<n;i++) T[i] = (PetscReal)(i+1);
  for (i=k;i<n-1;i++) T[i+ld] = 1.0;
  for (i=l;i<k;i++) T[i+2*ld] = 1.0;
  T[2*ld+l+1] = -7; T[ld+k+1] = -7;
  /* Signature matrix */
  for (i=0;i<n;i++) s[i] = 1.0;
  s[l+1] = -1.0;
  s[k+1] = -1.0;
  PetscCall(DSRestoreArrayReal(ds,DS_MAT_T,&T));
  PetscCall(DSRestoreArrayReal(ds,DS_MAT_D,&s));
  if (l==0 && k==0) PetscCall(DSSetState(ds,DS_STATE_INTERMEDIATE));
  else PetscCall(DSSetState(ds,DS_STATE_RAW));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Initial - - - - - - - - -\n"));
  PetscCall(DSView(ds,viewer));

  /* Solve */
  PetscCall(PetscCalloc2(n,&eigr,n,&eigi));
  PetscCall(DSGetSlepcSC(ds,&sc));
  sc->comparison    = SlepcCompareLargestMagnitude;
  sc->comparisonctx = NULL;
  sc->map           = NULL;
  sc->mapobj        = NULL;
  PetscCall(DSSolve(ds,eigr,eigi));
  PetscCall(DSSort(ds,eigr,eigi,NULL,NULL,NULL));
  if (verbose) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"After solve - - - - - - - - -\n"));
    PetscCall(DSView(ds,viewer));
  }

  /* Print eigenvalues */
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Computed eigenvalues =\n"));
  for (i=0;i<n;i++) {
#if PetscDefined(USE_COMPLEX)
    re = PetscRealPart(eigr[i]);
    im = PetscImaginaryPart(eigr[i]);
#else
    re = eigr[i];
    im = eigi[i];
#endif
    if (PetscAbs(im)<1e-10) PetscCall(PetscViewerASCIIPrintf(viewer,"  %.5f\n",(double)re));
    else PetscCall(PetscViewerASCIIPrintf(viewer,"  %.5f%+.5fi\n",(double)re,(double)im));
  }
  PetscCall(PetscFree2(eigr,eigi));
  PetscCall(DSDestroy(&ds));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   test:
      suffix: 1
      requires: !single
      args: -ds_method {{0 1 2}}
      filter: grep -v "solving the problem"

   test:
      suffix: 2
      args: -n 5 -k 4 -l 4 -ds_method {{0 1 2}}
      filter: grep -v "solving the problem"

TEST*/
