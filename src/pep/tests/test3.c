/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test PEP interface functions.\n\n";

#include <slepcpep.h>

int main(int argc,char **argv)
{
  Mat                A[3],B;      /* problem matrices */
  PEP                pep;         /* eigenproblem solver context */
  ST                 st;
  KSP                ksp;
  Vec                Dl,Dr;
  DS                 ds;
  PetscReal          tol,alpha;
  PetscScalar        target;
  PetscInt           n=20,i,its,nev,ncv,mpd,Istart,Iend,nmat;
  PEPWhich           which;
  PEPConvergedReason reason;
  PEPType            type;
  PEPExtract         extr;
  PEPBasis           basis;
  PEPScale           scale;
  PEPRefine          refine;
  PEPRefineScheme    rscheme;
  PEPConv            conv;
  PEPStop            stop;
  PEPProblemType     ptype;
  PetscViewerAndFormat *vf;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nDiagonal Quadratic Eigenproblem, n=%" PetscInt_FMT "\n\n",n));

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

  PetscCall(MatCreateVecs(A[0],&Dr,&Dl));
  PetscCall(VecSet(Dl,1.0));
  PetscCall(VecSet(Dr,0.95));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
             Create eigensolver and test interface functions
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(PEPCreate(PETSC_COMM_WORLD,&pep));
  PetscCall(PEPSetOperators(pep,3,A));
  PetscCall(PEPGetNumMatrices(pep,&nmat));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Polynomial of degree %" PetscInt_FMT "\n",nmat-1));
  PetscCall(PEPGetOperators(pep,0,&B));
  PetscCall(MatView(B,NULL));

  PetscCall(PEPSetType(pep,PEPTOAR));
  PetscCall(PEPGetType(pep,&type));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Type set to %s\n",type));

  PetscCall(PEPGetProblemType(pep,&ptype));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Problem type before changing = %d",(int)ptype));
  PetscCall(PEPSetProblemType(pep,PEP_HERMITIAN));
  PetscCall(PEPGetProblemType(pep,&ptype));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," ... changed to %d.\n",(int)ptype));

  PetscCall(PEPGetExtract(pep,&extr));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Extraction before changing = %d",(int)extr));
  PetscCall(PEPSetExtract(pep,PEP_EXTRACT_STRUCTURED));
  PetscCall(PEPGetExtract(pep,&extr));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," ... changed to %d\n",(int)extr));

  PetscCall(PEPSetScale(pep,PEP_SCALE_BOTH,.1,Dl,Dr,5,1.0));
  PetscCall(PEPGetScale(pep,&scale,&alpha,NULL,NULL,&its,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Scaling: %s, alpha=%g, its=%" PetscInt_FMT "\n",PEPScaleTypes[scale],(double)alpha,its));

  PetscCall(PEPSetBasis(pep,PEP_BASIS_CHEBYSHEV1));
  PetscCall(PEPGetBasis(pep,&basis));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Polynomial basis: %s\n",PEPBasisTypes[basis]));

  PetscCall(PEPSetRefine(pep,PEP_REFINE_SIMPLE,1,1e-9,2,PEP_REFINE_SCHEME_SCHUR));
  PetscCall(PEPGetRefine(pep,&refine,NULL,&tol,&its,&rscheme));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Refinement: %s, tol=%g, its=%" PetscInt_FMT ", scheme=%s\n",PEPRefineTypes[refine],(double)tol,its,PEPRefineSchemes[rscheme]));

  PetscCall(PEPSetTarget(pep,4.8));
  PetscCall(PEPGetTarget(pep,&target));
  PetscCall(PEPSetWhichEigenpairs(pep,PEP_TARGET_MAGNITUDE));
  PetscCall(PEPGetWhichEigenpairs(pep,&which));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Which = %d, target = %g\n",(int)which,(double)PetscRealPart(target)));

  PetscCall(PEPSetDimensions(pep,4,PETSC_DETERMINE,PETSC_DETERMINE));
  PetscCall(PEPGetDimensions(pep,&nev,&ncv,&mpd));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Dimensions: nev=%" PetscInt_FMT ", ncv=%" PetscInt_FMT ", mpd=%" PetscInt_FMT "\n",nev,ncv,mpd));

  PetscCall(PEPSetTolerances(pep,2.2e-4,200));
  PetscCall(PEPGetTolerances(pep,&tol,&its));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Tolerance = %.5f, max_its = %" PetscInt_FMT "\n",(double)tol,its));

  PetscCall(PEPSetConvergenceTest(pep,PEP_CONV_ABS));
  PetscCall(PEPGetConvergenceTest(pep,&conv));
  PetscCall(PEPSetStoppingTest(pep,PEP_STOP_BASIC));
  PetscCall(PEPGetStoppingTest(pep,&stop));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Convergence test = %d, stopping test = %d\n",(int)conv,(int)stop));

  PetscCall(PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_DEFAULT,&vf));
  PetscCall(PEPMonitorSet(pep,(PetscErrorCode (*)(PEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,void*))PEPMonitorFirst,vf,(PetscCtxDestroyFn*)PetscViewerAndFormatDestroy));
  PetscCall(PEPMonitorCancel(pep));

  PetscCall(PEPGetST(pep,&st));
  PetscCall(STGetKSP(st,&ksp));
  PetscCall(KSPSetTolerances(ksp,1e-8,1e-35,PETSC_CURRENT,PETSC_CURRENT));
  PetscCall(STView(st,NULL));
  PetscCall(PEPGetDS(pep,&ds));
  PetscCall(DSView(ds,NULL));

  PetscCall(PEPSetFromOptions(pep));
  PetscCall(PEPSolve(pep));
  PetscCall(PEPGetConvergedReason(pep,&reason));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Finished - converged reason = %d\n",(int)reason));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(PEPErrorView(pep,PEP_ERROR_RELATIVE,NULL));
  PetscCall(PEPDestroy(&pep));
  PetscCall(MatDestroy(&A[0]));
  PetscCall(MatDestroy(&A[1]));
  PetscCall(MatDestroy(&A[2]));
  PetscCall(VecDestroy(&Dl));
  PetscCall(VecDestroy(&Dr));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   test:
      suffix: 1
      args: -pep_tol 1e-6 -pep_ncv 22
      filter: sed -e "s/[+-]0\.0*i//g"

TEST*/
