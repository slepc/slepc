/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test NEP interface functions.\n\n";

#include <slepcnep.h>

int main(int argc,char **argv)
{
  Mat                  A[3],B;      /* problem matrices */
  FN                   f[3],g;      /* problem functions */
  NEP                  nep;         /* eigenproblem solver context */
  DS                   ds;
  RG                   rg;
  PetscReal            tol;
  PetscScalar          coeffs[2],target;
  PetscInt             n=20,i,its,nev,ncv,mpd,Istart,Iend,nterm;
  PetscBool            twoside;
  NEPWhich             which;
  NEPConvergedReason   reason;
  NEPType              type;
  NEPRefine            refine;
  NEPRefineScheme      rscheme;
  NEPConv              conv;
  NEPStop              stop;
  NEPProblemType       ptype;
  MatStructure         mstr;
  PetscViewerAndFormat *vf;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nDiagonal Nonlinear Eigenproblem, n=%" PetscInt_FMT "\n\n",n));

  /*
     Matrices
  */
  PetscCall(MatCreate(PETSC_COMM_WORLD,&A[0]));
  PetscCall(MatSetSizes(A[0],PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetFromOptions(A[0]));
  PetscCall(MatGetOwnershipRange(A[0],&Istart,&Iend));
  for (i=Istart;i<Iend;i++) PetscCall(MatSetValue(A[0],i,i,i+1,INSERT_VALUES));
  PetscCall(MatAssemblyBegin(A[0],MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A[0],MAT_FINAL_ASSEMBLY));

  PetscCall(MatCreate(PETSC_COMM_WORLD,&A[1]));
  PetscCall(MatSetSizes(A[1],PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetFromOptions(A[1]));
  PetscCall(MatGetOwnershipRange(A[1],&Istart,&Iend));
  for (i=Istart;i<Iend;i++) PetscCall(MatSetValue(A[1],i,i,1.0,INSERT_VALUES));
  PetscCall(MatAssemblyBegin(A[1],MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A[1],MAT_FINAL_ASSEMBLY));

  PetscCall(MatCreate(PETSC_COMM_WORLD,&A[2]));
  PetscCall(MatSetSizes(A[2],PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetFromOptions(A[2]));
  PetscCall(MatGetOwnershipRange(A[1],&Istart,&Iend));
  for (i=Istart;i<Iend;i++) PetscCall(MatSetValue(A[2],i,i,n/(PetscReal)(i+1),INSERT_VALUES));
  PetscCall(MatAssemblyBegin(A[2],MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A[2],MAT_FINAL_ASSEMBLY));

  /*
     Functions: f0=-lambda, f1=1.0, f2=sqrt(lambda)
  */
  PetscCall(FNCreate(PETSC_COMM_WORLD,&f[0]));
  PetscCall(FNSetType(f[0],FNRATIONAL));
  coeffs[0] = -1.0; coeffs[1] = 0.0;
  PetscCall(FNRationalSetNumerator(f[0],2,coeffs));

  PetscCall(FNCreate(PETSC_COMM_WORLD,&f[1]));
  PetscCall(FNSetType(f[1],FNRATIONAL));
  coeffs[0] = 1.0;
  PetscCall(FNRationalSetNumerator(f[1],1,coeffs));

  PetscCall(FNCreate(PETSC_COMM_WORLD,&f[2]));
  PetscCall(FNSetType(f[2],FNSQRT));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
             Create eigensolver and test interface functions
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(NEPCreate(PETSC_COMM_WORLD,&nep));
  PetscCall(NEPSetSplitOperator(nep,3,A,f,SAME_NONZERO_PATTERN));
  PetscCall(NEPGetSplitOperatorInfo(nep,&nterm,&mstr));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Nonlinear function with %" PetscInt_FMT " terms, with %s nonzero pattern\n",nterm,MatStructures[mstr]));
  PetscCall(NEPGetSplitOperatorTerm(nep,0,&B,&g));
  PetscCall(MatView(B,NULL));
  PetscCall(FNView(g,NULL));

  PetscCall(NEPSetType(nep,NEPRII));
  PetscCall(NEPGetType(nep,&type));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Type set to %s\n",type));
  PetscCall(NEPGetTwoSided(nep,&twoside));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Two-sided flag = %s\n",twoside?"true":"false"));

  PetscCall(NEPGetProblemType(nep,&ptype));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Problem type before changing = %d",(int)ptype));
  PetscCall(NEPSetProblemType(nep,NEP_RATIONAL));
  PetscCall(NEPGetProblemType(nep,&ptype));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," ... changed to %d.\n",(int)ptype));

  PetscCall(NEPSetRefine(nep,NEP_REFINE_SIMPLE,1,1e-9,2,NEP_REFINE_SCHEME_EXPLICIT));
  PetscCall(NEPGetRefine(nep,&refine,NULL,&tol,&its,&rscheme));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Refinement: %s, tol=%g, its=%" PetscInt_FMT ", scheme=%s\n",NEPRefineTypes[refine],(double)tol,its,NEPRefineSchemes[rscheme]));

  PetscCall(NEPSetTarget(nep,1.1));
  PetscCall(NEPGetTarget(nep,&target));
  PetscCall(NEPSetWhichEigenpairs(nep,NEP_TARGET_MAGNITUDE));
  PetscCall(NEPGetWhichEigenpairs(nep,&which));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Which = %d, target = %g\n",(int)which,(double)PetscRealPart(target)));

  PetscCall(NEPSetDimensions(nep,1,12,PETSC_CURRENT));
  PetscCall(NEPGetDimensions(nep,&nev,&ncv,&mpd));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Dimensions: nev=%" PetscInt_FMT ", ncv=%" PetscInt_FMT ", mpd=%" PetscInt_FMT "\n",nev,ncv,mpd));

  PetscCall(NEPSetTolerances(nep,1.0e-6,200));
  PetscCall(NEPGetTolerances(nep,&tol,&its));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Tolerance = %.6f, max_its = %" PetscInt_FMT "\n",(double)tol,its));

  PetscCall(NEPSetConvergenceTest(nep,NEP_CONV_ABS));
  PetscCall(NEPGetConvergenceTest(nep,&conv));
  PetscCall(NEPSetStoppingTest(nep,NEP_STOP_BASIC));
  PetscCall(NEPGetStoppingTest(nep,&stop));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Convergence test = %d, stopping test = %d\n",(int)conv,(int)stop));

  PetscCall(PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_DEFAULT,&vf));
  PetscCall(NEPMonitorSet(nep,(PetscErrorCode (*)(NEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,void*))NEPMonitorFirst,vf,(PetscCtxDestroyFn*)PetscViewerAndFormatDestroy));
  PetscCall(NEPMonitorCancel(nep));

  PetscCall(NEPGetDS(nep,&ds));
  PetscCall(DSView(ds,NULL));
  PetscCall(NEPSetFromOptions(nep));

  PetscCall(NEPGetRG(nep,&rg));
  PetscCall(RGView(rg,NULL));

  PetscCall(NEPSolve(nep));
  PetscCall(NEPGetConvergedReason(nep,&reason));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Finished - converged reason = %d\n",(int)reason));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(NEPErrorView(nep,NEP_ERROR_RELATIVE,NULL));
  PetscCall(NEPDestroy(&nep));
  PetscCall(MatDestroy(&A[0]));
  PetscCall(MatDestroy(&A[1]));
  PetscCall(MatDestroy(&A[2]));
  PetscCall(FNDestroy(&f[0]));
  PetscCall(FNDestroy(&f[1]));
  PetscCall(FNDestroy(&f[2]));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   test:
      suffix: 1
      args: -nep_view
      filter: grep -v tolerance

TEST*/
