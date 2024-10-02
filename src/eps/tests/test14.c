/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test EPS interface functions.\n\n";

#include <slepceps.h>

int main(int argc,char **argv)
{
  Mat                A,B;         /* problem matrix */
  EPS                eps;         /* eigenproblem solver context */
  ST                 st;
  KSP                ksp;
  DS                 ds;
  PetscReal          cut,tol;
  PetscScalar        target;
  PetscInt           n=20,i,its,nev,ncv,mpd,Istart,Iend;
  PetscBool          flg,pur,track;
  EPSConvergedReason reason;
  EPSType            type;
  EPSExtraction      extr;
  EPSBalance         bal;
  EPSWhich           which;
  EPSConv            conv;
  EPSStop            stop;
  EPSProblemType     ptype;
  PetscViewerAndFormat *vf;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nDiagonal Eigenproblem, n=%" PetscInt_FMT "\n\n",n));

  PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
  PetscCall(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetFromOptions(A));
  PetscCall(MatGetOwnershipRange(A,&Istart,&Iend));
  for (i=Istart;i<Iend;i++) PetscCall(MatSetValue(A,i,i,i+1,INSERT_VALUES));
  PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
             Create eigensolver and test interface functions
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(EPSCreate(PETSC_COMM_WORLD,&eps));
  PetscCall(EPSSetOperators(eps,A,NULL));
  PetscCall(EPSGetOperators(eps,&B,NULL));
  PetscCall(MatView(B,NULL));

  PetscCall(EPSSetType(eps,EPSKRYLOVSCHUR));
  PetscCall(EPSGetType(eps,&type));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Type set to %s\n",type));

  PetscCall(EPSGetProblemType(eps,&ptype));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Problem type before changing = %d",(int)ptype));
  PetscCall(EPSSetProblemType(eps,EPS_HEP));
  PetscCall(EPSGetProblemType(eps,&ptype));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," ... changed to %d.",(int)ptype));
  PetscCall(EPSIsGeneralized(eps,&flg));
  if (flg) PetscCall(PetscPrintf(PETSC_COMM_WORLD," generalized"));
  PetscCall(EPSIsHermitian(eps,&flg));
  if (flg) PetscCall(PetscPrintf(PETSC_COMM_WORLD," hermitian"));
  PetscCall(EPSIsPositive(eps,&flg));
  if (flg) PetscCall(PetscPrintf(PETSC_COMM_WORLD," positive"));

  PetscCall(EPSGetExtraction(eps,&extr));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n Extraction before changing = %d",(int)extr));
  PetscCall(EPSSetExtraction(eps,EPS_HARMONIC));
  PetscCall(EPSGetExtraction(eps,&extr));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," ... changed to %d\n",(int)extr));

  PetscCall(EPSSetBalance(eps,EPS_BALANCE_ONESIDE,8,1e-6));
  PetscCall(EPSGetBalance(eps,&bal,&its,&cut));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Balance: %s, its=%" PetscInt_FMT ", cutoff=%g\n",EPSBalanceTypes[bal],its,(double)cut));

  PetscCall(EPSSetPurify(eps,PETSC_FALSE));
  PetscCall(EPSGetPurify(eps,&pur));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Eigenvector purification: %s\n",pur?"on":"off"));
  PetscCall(EPSGetTrackAll(eps,&track));

  PetscCall(EPSSetTarget(eps,4.8));
  PetscCall(EPSGetTarget(eps,&target));
  PetscCall(EPSSetWhichEigenpairs(eps,EPS_TARGET_MAGNITUDE));
  PetscCall(EPSGetWhichEigenpairs(eps,&which));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Which = %d, target = %g\n",(int)which,(double)PetscRealPart(target)));

  PetscCall(EPSSetDimensions(eps,4,PETSC_DETERMINE,PETSC_DETERMINE));
  PetscCall(EPSGetDimensions(eps,&nev,&ncv,&mpd));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Dimensions: nev=%" PetscInt_FMT ", ncv=%" PetscInt_FMT ", mpd=%" PetscInt_FMT "\n",nev,ncv,mpd));

  PetscCall(EPSSetTolerances(eps,2.2e-4,200));
  PetscCall(EPSGetTolerances(eps,&tol,&its));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Tolerance = %.5f, max_its = %" PetscInt_FMT "\n",(double)tol,its));

  PetscCall(EPSSetConvergenceTest(eps,EPS_CONV_ABS));
  PetscCall(EPSGetConvergenceTest(eps,&conv));
  PetscCall(EPSSetStoppingTest(eps,EPS_STOP_BASIC));
  PetscCall(EPSGetStoppingTest(eps,&stop));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Convergence test = %d, stopping test = %d\n",(int)conv,(int)stop));

  PetscCall(PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_DEFAULT,&vf));
  PetscCall(EPSMonitorSet(eps,(PetscErrorCode (*)(EPS,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,void*))EPSMonitorFirst,vf,(PetscCtxDestroyFn*)PetscViewerAndFormatDestroy));
  PetscCall(EPSMonitorCancel(eps));

  PetscCall(EPSGetST(eps,&st));
  PetscCall(STGetKSP(st,&ksp));
  PetscCall(KSPSetTolerances(ksp,1e-8,1e-35,PETSC_CURRENT,PETSC_CURRENT));
  PetscCall(STView(st,NULL));
  PetscCall(EPSGetDS(eps,&ds));
  PetscCall(DSView(ds,NULL));

  PetscCall(EPSSetFromOptions(eps));
  PetscCall(EPSSolve(eps));
  PetscCall(EPSGetConvergedReason(eps,&reason));
  PetscCall(EPSGetIterationNumber(eps,&its));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Finished - converged reason = %d, its=%" PetscInt_FMT "\n",(int)reason,its));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(EPSErrorView(eps,EPS_ERROR_RELATIVE,NULL));
  PetscCall(EPSDestroy(&eps));
  PetscCall(MatDestroy(&A));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   test:
      suffix: 1
      args: -eps_ncv 14
      filter: sed -e "s/00001/00000/" | sed -e "s/4.99999/5.00000/" | sed -e "s/5.99999/6.00000/"

TEST*/
