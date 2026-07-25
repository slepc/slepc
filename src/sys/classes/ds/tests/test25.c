/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test for DSPEP and DSNEP.\n\n";

#include <slepcds.h>

#define NMAT 5

int main(int argc,char **argv)
{
  DS             ds;
  FN             f[NMAT],qfun;
  SlepcSC        sc;
  PetscScalar    *A,*wr,*wi,*X,*y,*r,numer[NMAT],alpha;
  PetscReal      c[10] = { 0.6, 1.3, 1.3, 0.1, 0.1, 1.2, 1.0, 1.0, 1.2, 1.0 };
  PetscReal      tol,radius=1.5,re,im,nrm;
  PetscInt       i,j,ii,jj,II,k,m=3,n,ld,nev,nfun,d,*inside;
  PetscViewer    viewer;
  PetscBool      verbose,isnep=PETSC_FALSE;
  RG             rg;
  DSMatType      mat[5]={DS_MAT_E0,DS_MAT_E1,DS_MAT_E2,DS_MAT_E3,DS_MAT_E4};
#if !PetscDefined(USE_COMPLEX)
  PetscScalar    *yi,*ri,alphai=0.0,t;
#endif

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-m",&m,NULL));
  PetscCall(PetscOptionsGetBool(NULL,NULL,"-isnep",&isnep,NULL));
  n = m*m;
  k = 10;
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nButterfly problem, n=%" PetscInt_FMT " (m=%" PetscInt_FMT ")\n\n",n,m));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-verbose",&verbose));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-radius",&radius,NULL));

  /* Create DS object */
  PetscCall(DSCreate(PETSC_COMM_WORLD,&ds));
  tol  = 1000*n*PETSC_MACHINE_EPSILON;
  if (isnep) {
    PetscCall(DSSetType(ds,DSNEP));
    PetscCall(DSSetMethod(ds,1));
    PetscCall(DSNEPSetRefine(ds,tol,PETSC_DECIDE));
  } else PetscCall(DSSetType(ds,DSPEP));
  PetscCall(DSSetFromOptions(ds));

  /* Set functions (prior to DSAllocate) f_i=x^i */
  if (isnep) {
    numer[0] = 1.0;
    for (j=1;j<NMAT;j++) numer[j] = 0.0;
    for (i=0;i<NMAT;i++) {
      PetscCall(FNCreate(PETSC_COMM_WORLD,&f[i]));
      PetscCall(FNSetType(f[i],FNRATIONAL));
      PetscCall(FNRationalSetNumerator(f[i],i+1,numer));
    }
    PetscCall(DSNEPSetFN(ds,NMAT,f));
  } else PetscCall(DSPEPSetDegree(ds,NMAT-1));

  /* Set dimensions */
  ld = n+2;  /* test leading dimension larger than n */
  PetscCall(DSAllocate(ds,ld));
  PetscCall(DSSetDimensions(ds,n,0,0));

  /* Set region (used only in method=1) */
  PetscCall(RGCreate(PETSC_COMM_WORLD,&rg));
  PetscCall(RGSetType(rg,RGELLIPSE));
  PetscCall(RGEllipseSetParameters(rg,1.5,radius,.5));
  PetscCall(RGSetFromOptions(rg));
  if (isnep) PetscCall(DSNEPSetRG(ds,rg));

  /* Set up viewer */
  PetscCall(PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&viewer));
  PetscCall(DSViewFromOptions(ds,NULL,"-ds_view"));
  if (verbose) {
    PetscCall(PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB));
    /* Show info about functions */
    if (isnep) {
      PetscCall(DSNEPGetNumFN(ds,&nfun));
      for (i=0;i<nfun;i++) {
        PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Function %" PetscInt_FMT ":\n",i));
        PetscCall(DSNEPGetFN(ds,i,&qfun));
        PetscCall(FNView(qfun,NULL));
      }
    }
  }

  /* Fill matrices */
  /* A0 */
  PetscCall(DSGetArray(ds,DS_MAT_E0,&A));
  for (II=0;II<n;II++) {
    i = II/m; j = II-i*m;
    A[II+II*ld] = 4.0*c[0]/6.0+4.0*c[1]/6.0;
    if (j>0) A[II+(II-1)*ld] = c[0]/6.0;
    if (j<m-1) A[II+ld*(II+1)] = c[0]/6.0;
    if (i>0) A[II+ld*(II-m)] = c[1]/6.0;
    if (i<m-1) A[II+ld*(II+m)] = c[1]/6.0;
  }
  PetscCall(DSRestoreArray(ds,DS_MAT_E0,&A));

  /* A1 */
  PetscCall(DSGetArray(ds,DS_MAT_E1,&A));
  for (II=0;II<n;II++) {
    i = II/m; j = II-i*m;
    if (j>0) A[II+ld*(II-1)] = c[2];
    if (j<m-1) A[II+ld*(II+1)] = -c[2];
    if (i>0) A[II+ld*(II-m)] = c[3];
    if (i<m-1) A[II+ld*(II+m)] = -c[3];
  }
  PetscCall(DSRestoreArray(ds,DS_MAT_E1,&A));

  /* A2 */
  PetscCall(DSGetArray(ds,DS_MAT_E2,&A));
  for (II=0;II<n;II++) {
    i = II/m; j = II-i*m;
    A[II+ld*II] = -2.0*c[4]-2.0*c[5];
    if (j>0) A[II+ld*(II-1)] = c[4];
    if (j<m-1) A[II+ld*(II+1)] = c[4];
    if (i>0) A[II+ld*(II-m)] = c[5];
    if (i<m-1) A[II+ld*(II+m)] = c[5];
  }
  PetscCall(DSRestoreArray(ds,DS_MAT_E2,&A));

  /* A3 */
  PetscCall(DSGetArray(ds,DS_MAT_E3,&A));
  for (II=0;II<n;II++) {
    i = II/m; j = II-i*m;
    if (j>0) A[II+ld*(II-1)] = c[6];
    if (j<m-1) A[II+ld*(II+1)] = -c[6];
    if (i>0) A[II+ld*(II-m)] = c[7];
    if (i<m-1) A[II+ld*(II+m)] = -c[7];
  }
  PetscCall(DSRestoreArray(ds,DS_MAT_E3,&A));

  /* A4 */
  PetscCall(DSGetArray(ds,DS_MAT_E4,&A));
  for (II=0;II<n;II++) {
    i = II/m; j = II-i*m;
    A[II+ld*II] = 2.0*c[8]+2.0*c[9];
    if (j>0) A[II+ld*(II-1)] = -c[8];
    if (j<m-1) A[II+ld*(II+1)] = -c[8];
    if (i>0) A[II+ld*(II-m)] = -c[9];
    if (i<m-1) A[II+ld*(II+m)] = -c[9];
  }
  PetscCall(DSRestoreArray(ds,DS_MAT_E4,&A));

  if (verbose) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Initial - - - - - - - - -\n"));
    PetscCall(DSView(ds,viewer));
  }

  /* Solve */
  if (isnep) PetscCall(DSNEPGetMinimality(ds,&d));
  else PetscCall(DSPEPGetDegree(ds,&d));
  PetscCall(PetscCalloc3(n*d,&wr,n*d,&wi,n*d,&inside));
  PetscCall(DSGetSlepcSC(ds,&sc));
  sc->comparison    = SlepcCompareLargestMagnitude;
  sc->comparisonctx = NULL;
  sc->map           = NULL;
  sc->mapobj        = NULL;
  PetscCall(DSSolve(ds,wr,wi));
  PetscCall(DSSort(ds,wr,wi,NULL,NULL,NULL));

  if (verbose) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"After solve - - - - - - - - -\n"));
    PetscCall(DSView(ds,viewer));
  }
  if (isnep) {
    PetscCall(DSGetDimensions(ds,NULL,NULL,NULL,&nev));
    for (i=0;i<nev;i++) inside[i] = i;
  } else {
    PetscCall(RGCheckInside(rg,d*n,wr,wi,inside));
    nev = 0;
    for (i=0;i<d*n;i++) if (inside[i]>0) inside[nev++] = i;
  }

  /* Print computed eigenvalues */
  PetscCall(PetscMalloc2(ld,&y,ld,&r));
#if !PetscDefined(USE_COMPLEX)
  PetscCall(PetscMalloc2(ld,&yi,ld,&ri));
#endif
  PetscCall(DSVectors(ds,DS_MAT_X,NULL,NULL));
  PetscCall(DSGetArray(ds,DS_MAT_X,&X));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Computed eigenvalues in the region: %" PetscInt_FMT "\n",nev));
  for (i=0;i<nev;i++) {
#if PetscDefined(USE_COMPLEX)
    re = PetscRealPart(wr[inside[i]]);
    im = PetscImaginaryPart(wr[inside[i]]);
#else
    re = wr[inside[i]];
    im = wi[inside[i]];
#endif
    PetscCall(PetscArrayzero(r,n));
#if !PetscDefined(USE_COMPLEX)
    PetscCall(PetscArrayzero(ri,n));
#endif
    /* Residual */
    alpha = 1.0;
    for (k=0;k<NMAT;k++) {
      PetscCall(DSGetArray(ds,mat[k],&A));
      for (ii=0;ii<n;ii++) {
        y[ii] = 0.0;
        for (jj=0;jj<n;jj++) y[ii] += A[jj*ld+ii]*X[inside[i]*ld+jj];
      }
#if !PetscDefined(USE_COMPLEX)
      for (ii=0;ii<n;ii++) {
        yi[ii] = 0.0;
        for (jj=0;jj<n;jj++) yi[ii] += A[jj*ld+ii]*X[inside[i+1]*ld+jj];
      }
#endif
      PetscCall(DSRestoreArray(ds,mat[k],&A));
      if (isnep) PetscCall(FNEvaluateFunction(f[k],wr[inside[i]],&alpha));
      for (ii=0;ii<n;ii++) r[ii] += alpha*y[ii];
#if !PetscDefined(USE_COMPLEX)
      for (ii=0;ii<n;ii++) r[ii]  -= alphai*yi[ii];
      for (ii=0;ii<n;ii++) ri[ii] += alpha*yi[ii]+alphai*y[ii];
#endif
      if (!isnep) {
#if PetscDefined(USE_COMPLEX)
        alpha *= wr[inside[i]];
#else
        t      = alpha;
        alpha  = alpha*re-alphai*im;
        alphai = alphai*re+t*im;
#endif
      }
    }
    nrm = 0.0;
    for (k=0;k<n;k++) {
#if !PetscDefined(USE_COMPLEX)
      nrm += r[k]*r[k]+ri[k]*ri[k];
#else
      nrm += PetscRealPart(r[k]*PetscConj(r[k]));
#endif
    }
    nrm = PetscSqrtReal(nrm);
    if (nrm/SlepcAbsEigenvalue(wr[inside[i]],wi[inside[i]])>tol) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Warning: the residual norm of the %" PetscInt_FMT "-th computed eigenpair %g\n",i,(double)nrm));
    if (PetscAbs(im)<1e-10) PetscCall(PetscViewerASCIIPrintf(viewer,"  %.5f\n",(double)re));
    else PetscCall(PetscViewerASCIIPrintf(viewer,"  %.5f%+.5fi\n",(double)re,(double)im));
#if !PetscDefined(USE_COMPLEX)
    if (im!=0.0) i++;
    if (PetscAbs(im)<1e-10) PetscCall(PetscViewerASCIIPrintf(viewer,"  %.5f\n",(double)re));
    else PetscCall(PetscViewerASCIIPrintf(viewer,"  %.5f%+.5fi\n",(double)re,(double)-im));
#endif
  }
  PetscCall(DSRestoreArray(ds,DS_MAT_X,&X));
  PetscCall(PetscFree3(wr,wi,inside));
  PetscCall(PetscFree2(y,r));
#if !PetscDefined(USE_COMPLEX)
  PetscCall(PetscFree2(yi,ri));
#endif
  if (isnep) {
    for (i=0;i<NMAT;i++) PetscCall(FNDestroy(&f[i]));
  }
  PetscCall(DSDestroy(&ds));
  PetscCall(RGDestroy(&rg));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   testset:
      filter: sed -e "s/[+-]\([0-9]\.[0-9]*i\)/+-\\1/" | sed -e "s/56808/56807/" | sed -e "s/34719/34720/"
      output_file: output/test25_1.out
      test:
         suffix: 1
      test:
         suffix: 2
         args: -isnep
         requires: complex !single

TEST*/
