/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test some NLEIGS interface functions.\n\n"
  "Based on ex27.c. The command line options are:\n"
  "  -n <n>, where <n> = matrix dimension.\n";

/*
   Solve T(lambda)x=0 using NLEIGS solver
      with T(lambda) = -D+sqrt(lambda)*I
      where D is the Laplacian operator in 1 dimension
      and with the interpolation interval [.01,16]
*/

#include <slepcnep.h>

/*
   User-defined routines
*/
PetscErrorCode ComputeSingularities(NEP,PetscInt*,PetscScalar*,void*);

int main(int argc,char **argv)
{
  NEP            nep;             /* nonlinear eigensolver context */
  Mat            A[2];
  PetscInt       n=100,Istart,Iend,i,ns,nsin;
  PetscBool      terse,fb;
  RG             rg;
  FN             f[2];
  PetscScalar    coeffs,shifts[]={1.06,1.1,1.12,1.15},*rkshifts,val;
  PetscErrorCode (*fsing)(NEP,PetscInt*,PetscScalar*,void*);

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nSquare root eigenproblem, n=%" PetscInt_FMT "\n\n",n));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create nonlinear eigensolver and set some options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(NEPCreate(PETSC_COMM_WORLD,&nep));
  PetscCall(NEPSetType(nep,NEPNLEIGS));
  PetscCall(NEPNLEIGSSetSingularitiesFunction(nep,ComputeSingularities,NULL));
  PetscCall(NEPGetRG(nep,&rg));
  PetscCall(RGSetType(rg,RGINTERVAL));
#if defined(PETSC_USE_COMPLEX)
  PetscCall(RGIntervalSetEndpoints(rg,0.01,16.0,-0.001,0.001));
#else
  PetscCall(RGIntervalSetEndpoints(rg,0.01,16.0,0,0));
#endif
  PetscCall(NEPSetTarget(nep,1.1));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Define the nonlinear problem in split form
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* Create matrices */
  PetscCall(MatCreate(PETSC_COMM_WORLD,&A[0]));
  PetscCall(MatSetSizes(A[0],PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetFromOptions(A[0]));
  PetscCall(MatGetOwnershipRange(A[0],&Istart,&Iend));
  for (i=Istart;i<Iend;i++) {
    if (i>0) PetscCall(MatSetValue(A[0],i,i-1,1.0,INSERT_VALUES));
    if (i<n-1) PetscCall(MatSetValue(A[0],i,i+1,1.0,INSERT_VALUES));
    PetscCall(MatSetValue(A[0],i,i,-2.0,INSERT_VALUES));
  }
  PetscCall(MatAssemblyBegin(A[0],MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A[0],MAT_FINAL_ASSEMBLY));

  PetscCall(MatCreateConstantDiagonal(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,n,n,1.0,&A[1]));

  /* Define functions */
  PetscCall(FNCreate(PETSC_COMM_WORLD,&f[0]));
  PetscCall(FNSetType(f[0],FNRATIONAL));
  coeffs = 1.0;
  PetscCall(FNRationalSetNumerator(f[0],1,&coeffs));
  PetscCall(FNCreate(PETSC_COMM_WORLD,&f[1]));
  PetscCall(FNSetType(f[1],FNSQRT));
  PetscCall(NEPSetSplitOperator(nep,2,A,f,SUBSET_NONZERO_PATTERN));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                        Set some options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(NEPNLEIGSSetFullBasis(nep,PETSC_FALSE));
  PetscCall(NEPNLEIGSSetRKShifts(nep,4,shifts));
  PetscCall(NEPSetFromOptions(nep));

  PetscCall(NEPNLEIGSGetFullBasis(nep,&fb));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Using full basis = %s\n",fb?"true":"false"));
  PetscCall(NEPNLEIGSGetRKShifts(nep,&ns,&rkshifts));
  if (ns) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD," Using %" PetscInt_FMT " RK shifts =",ns));
    for (i=0;i<ns;i++) PetscCall(PetscPrintf(PETSC_COMM_WORLD," %g",(double)PetscRealPart(rkshifts[i])));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n"));
    PetscCall(PetscFree(rkshifts));
  }
  PetscCall(NEPNLEIGSGetSingularitiesFunction(nep,&fsing,NULL));
  nsin = 1;
  PetscCall((*fsing)(nep,&nsin,&val,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," First returned singularity = %g\n",(double)PetscRealPart(val)));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(NEPSolve(nep));

  /* show detailed info unless -terse option is given by user */
  PetscCall(PetscOptionsHasName(NULL,NULL,"-terse",&terse));
  if (terse) PetscCall(NEPErrorView(nep,NEP_ERROR_BACKWARD,NULL));
  else {
    PetscCall(PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL));
    PetscCall(NEPConvergedReasonView(nep,PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(NEPErrorView(nep,NEP_ERROR_BACKWARD,PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD));
  }
  PetscCall(NEPDestroy(&nep));
  PetscCall(MatDestroy(&A[0]));
  PetscCall(MatDestroy(&A[1]));
  PetscCall(FNDestroy(&f[0]));
  PetscCall(FNDestroy(&f[1]));
  PetscCall(SlepcFinalize());
  return 0;
}

/* ------------------------------------------------------------------- */
/*
   ComputeSingularities - Computes maxnp points (at most) in the complex plane where
   the function T(.) is not analytic.

   In this case, we discretize the singularity region (-inf,0)~(-1e+6,-1e-5)
*/
PetscErrorCode ComputeSingularities(NEP nep,PetscInt *maxnp,PetscScalar *xi,void *pt)
{
  PetscReal h;
  PetscInt  i;

  PetscFunctionBeginUser;
  h = 11.0/(*maxnp-1);
  xi[0] = -1e-5; xi[*maxnp-1] = -1e+6;
  for (i=1;i<*maxnp-1;i++) xi[i] = -PetscPowReal(10,-5+h*i);
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*TEST

   test:
      suffix: 1
      args: -nep_nev 3 -nep_nleigs_interpolation_degree 20 -terse -nep_view
      requires: double
      filter: grep -v tolerance | sed -e "s/[+-]0\.0*i//g"

TEST*/
