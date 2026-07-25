/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   This example implements one of the problems found at
       NLEVP: A Collection of Nonlinear Eigenvalue Problems,
       The University of Manchester.
   The details of the collection can be found at:
       [1] T. Betcke et al., "NLEVP: A Collection of Nonlinear Eigenvalue
           Problems", ACM Trans. Math. Software 39(2), Article 7, 2013.

   The loaded_string problem is a rational eigenvalue problem for the
   finite element model of a loaded vibrating string.
*/

static char help[] = "Illustrates computation of left eigenvectors and resolvent.\n\n"
  "This is based on loaded_string from the NLEVP collection.\n"
  "The command line options are:\n"
  "  -n <n>, dimension of the matrices.\n"
  "  -kappa <kappa>, stiffness of elastic spring.\n"
  "  -mass <m>, mass of the attached load.\n\n";

#include <slepcnep.h>

#define NMAT 3

int main(int argc,char **argv)
{
  Mat            A[NMAT];         /* problem matrices */
  FN             f[NMAT];         /* functions to define the nonlinear operator */
  NEP            nep;             /* nonlinear eigensolver context */
  RG             rg;
  Vec            v,r,z,w;
  PetscInt       n=100,Istart,Iend,i,nconv;
  PetscReal      kappa=1.0,m=1.0,nrm,tol;
  PetscScalar    lambda,sigma,numer[2],denom[2],omega1,omega2;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));

  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-kappa",&kappa,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-mass",&m,NULL));
  sigma = kappa/m;
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Loaded vibrating string, n=%" PetscInt_FMT " kappa=%g m=%g\n\n",n,(double)kappa,(double)m));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                       Build the problem matrices
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* initialize matrices */
  for (i=0;i<NMAT;i++) {
    PetscCall(MatCreate(PETSC_COMM_WORLD,&A[i]));
    PetscCall(MatSetSizes(A[i],PETSC_DECIDE,PETSC_DECIDE,n,n));
    PetscCall(MatSetFromOptions(A[i]));
  }
  PetscCall(MatGetOwnershipRange(A[0],&Istart,&Iend));

  /* A0 */
  for (i=Istart;i<Iend;i++) {
    PetscCall(MatSetValue(A[0],i,i,(i==n-1)?1.0*n:2.0*n,INSERT_VALUES));
    if (i>0) PetscCall(MatSetValue(A[0],i,i-1,-1.0*n,INSERT_VALUES));
    if (i<n-1) PetscCall(MatSetValue(A[0],i,i+1,-1.0*n,INSERT_VALUES));
  }

  /* A1 */
  for (i=Istart;i<Iend;i++) {
    PetscCall(MatSetValue(A[1],i,i,(i==n-1)?2.0/(6.0*n):4.0/(6.0*n),INSERT_VALUES));
    if (i>0) PetscCall(MatSetValue(A[1],i,i-1,1.0/(6.0*n),INSERT_VALUES));
    if (i<n-1) PetscCall(MatSetValue(A[1],i,i+1,1.0/(6.0*n),INSERT_VALUES));
  }

  /* A2 */
  if (Istart<=n-1 && n-1<Iend) PetscCall(MatSetValue(A[2],n-1,n-1,kappa,INSERT_VALUES));

  /* assemble matrices */
  for (i=0;i<NMAT;i++) PetscCall(MatAssemblyBegin(A[i],MAT_FINAL_ASSEMBLY));
  for (i=0;i<NMAT;i++) PetscCall(MatAssemblyEnd(A[i],MAT_FINAL_ASSEMBLY));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                       Create the problem functions
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* f1=1 */
  PetscCall(FNCreate(PETSC_COMM_WORLD,&f[0]));
  PetscCall(FNSetType(f[0],FNRATIONAL));
  numer[0] = 1.0;
  PetscCall(FNRationalSetNumerator(f[0],1,numer));

  /* f2=-lambda */
  PetscCall(FNCreate(PETSC_COMM_WORLD,&f[1]));
  PetscCall(FNSetType(f[1],FNRATIONAL));
  numer[0] = -1.0; numer[1] = 0.0;
  PetscCall(FNRationalSetNumerator(f[1],2,numer));

  /* f3=lambda/(lambda-sigma) */
  PetscCall(FNCreate(PETSC_COMM_WORLD,&f[2]));
  PetscCall(FNSetType(f[2],FNRATIONAL));
  numer[0] = 1.0; numer[1] = 0.0;
  denom[0] = 1.0; denom[1] = -sigma;
  PetscCall(FNRationalSetNumerator(f[2],2,numer));
  PetscCall(FNRationalSetDenominator(f[2],2,denom));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and solve the problem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(NEPCreate(PETSC_COMM_WORLD,&nep));
  PetscCall(NEPSetSplitOperator(nep,3,A,f,SUBSET_NONZERO_PATTERN));
  PetscCall(NEPSetProblemType(nep,NEP_RATIONAL));
  PetscCall(NEPSetDimensions(nep,8,PETSC_DETERMINE,PETSC_DETERMINE));

  /* set two-sided NLEIGS solver */
  PetscCall(NEPSetType(nep,NEPNLEIGS));
  PetscCall(NEPNLEIGSSetFullBasis(nep,PETSC_TRUE));
  PetscCall(NEPSetTwoSided(nep,PETSC_TRUE));
  PetscCall(NEPGetRG(nep,&rg));
  PetscCall(RGSetType(rg,RGINTERVAL));
#if PetscDefined(USE_COMPLEX)
  PetscCall(RGIntervalSetEndpoints(rg,4.0,700.0,-0.001,0.001));
#else
  PetscCall(RGIntervalSetEndpoints(rg,4.0,700.0,0,0));
#endif
  PetscCall(NEPSetTarget(nep,5.0));

  PetscCall(NEPSetFromOptions(nep));
  PetscCall(NEPSolve(nep));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                       Check left residual
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(MatCreateVecs(A[0],&v,&r));
  PetscCall(VecDuplicate(v,&w));
  PetscCall(VecDuplicate(v,&z));
  PetscCall(NEPGetConverged(nep,&nconv));
  PetscCall(NEPGetTolerances(nep,&tol,NULL));
  for (i=0;i<nconv;i++) {
    PetscCall(NEPGetEigenpair(nep,i,&lambda,NULL,NULL,NULL));
    PetscCall(NEPGetLeftEigenvector(nep,i,v,NULL));
    PetscCall(NEPApplyAdjoint(nep,lambda,v,w,r,NULL,NULL));
    PetscCall(VecNorm(r,NORM_2,&nrm));
    if (nrm>tol*PetscAbsScalar(lambda)) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Left residual i=%" PetscInt_FMT " is above tolerance --> %g\n",i,(double)(nrm/PetscAbsScalar(lambda))));
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Operate with resolvent
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  omega1 = 20.0;
  omega2 = 150.0;
  PetscCall(VecSet(v,0.0));
  PetscCall(VecSetValue(v,0,-1.0,INSERT_VALUES));
  PetscCall(VecSetValue(v,1,3.0,INSERT_VALUES));
  PetscCall(VecAssemblyBegin(v));
  PetscCall(VecAssemblyEnd(v));
  PetscCall(NEPApplyResolvent(nep,NULL,omega1,v,r));
  PetscCall(VecNorm(r,NORM_2,&nrm));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"resolvent, omega=%g: norm of computed vector=%g\n",(double)PetscRealPart(omega1),(double)nrm));
  PetscCall(NEPApplyResolvent(nep,NULL,omega2,v,r));
  PetscCall(VecNorm(r,NORM_2,&nrm));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"resolvent, omega=%g: norm of computed vector=%g\n",(double)PetscRealPart(omega2),(double)nrm));
  PetscCall(VecSet(v,1.0));
  PetscCall(NEPApplyResolvent(nep,NULL,omega1,v,r));
  PetscCall(VecNorm(r,NORM_2,&nrm));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"resolvent, omega=%g: norm of computed vector=%g\n",(double)PetscRealPart(omega1),(double)nrm));
  PetscCall(NEPApplyResolvent(nep,NULL,omega2,v,r));
  PetscCall(VecNorm(r,NORM_2,&nrm));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"resolvent, omega=%g: norm of computed vector=%g\n",(double)PetscRealPart(omega2),(double)nrm));

  /* clean up */
  PetscCall(NEPDestroy(&nep));
  for (i=0;i<NMAT;i++) {
    PetscCall(MatDestroy(&A[i]));
    PetscCall(FNDestroy(&f[i]));
  }
  PetscCall(VecDestroy(&v));
  PetscCall(VecDestroy(&r));
  PetscCall(VecDestroy(&w));
  PetscCall(VecDestroy(&z));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   test:
      suffix: 1
      requires: !single

TEST*/
