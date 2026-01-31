/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test LME interface functions, based on ex32.c.\n\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = number of grid subdivisions in x dimension.\n"
  "  -m <m>, where <m> = number of grid subdivisions in y dimension.\n\n";

#include <slepclme.h>

int main(int argc,char **argv)
{
  Mat                  A,B,C,C1,D;
  LME                  lme;
  PetscReal            tol,errest,error;
  PetscScalar          *u;
  PetscInt             N,n=10,m,Istart,Iend,II,maxit,ncv,i,j;
  PetscBool            flg,testprefix=PETSC_FALSE,viewmatrices=PETSC_FALSE;
  const char           *prefix;
  LMEType              type;
  LMEProblemType       ptype;
  PetscViewerAndFormat *vf;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));

  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-m",&m,&flg));
  if (!flg) m=n;
  N = n*m;
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nLyapunov equation, N=%" PetscInt_FMT " (%" PetscInt_FMT "x%" PetscInt_FMT " grid)\n\n",N,n,m));
  PetscCall(PetscOptionsGetBool(NULL,NULL,"-test_prefix",&testprefix,NULL));
  PetscCall(PetscOptionsGetBool(NULL,NULL,"-view_matrices",&viewmatrices,NULL));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                       Create the 2-D Laplacian, A
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
  PetscCall(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,N,N));
  PetscCall(MatSetFromOptions(A));
  PetscCall(MatGetOwnershipRange(A,&Istart,&Iend));
  for (II=Istart;II<Iend;II++) {
    i = II/n; j = II-i*n;
    if (i>0) PetscCall(MatSetValue(A,II,II-n,1.0,INSERT_VALUES));
    if (i<m-1) PetscCall(MatSetValue(A,II,II+n,1.0,INSERT_VALUES));
    if (j>0) PetscCall(MatSetValue(A,II,II-1,1.0,INSERT_VALUES));
    if (j<n-1) PetscCall(MatSetValue(A,II,II+1,1.0,INSERT_VALUES));
    PetscCall(MatSetValue(A,II,II,-4.0,INSERT_VALUES));
  }
  PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Create a low-rank Mat to store the right-hand side C = C1*C1'
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(MatCreate(PETSC_COMM_WORLD,&C1));
  PetscCall(MatSetSizes(C1,PETSC_DECIDE,PETSC_DECIDE,N,2));
  PetscCall(MatSetType(C1,MATDENSE));
  PetscCall(MatGetOwnershipRange(C1,&Istart,&Iend));
  PetscCall(MatDenseGetArray(C1,&u));
  for (i=Istart;i<Iend;i++) {
    if (i<N/2) u[i-Istart] = 1.0;
    if (i==0) u[i+Iend-2*Istart] = -2.0;
    if (i==1) u[i+Iend-2*Istart] = -1.0;
    if (i==2) u[i+Iend-2*Istart] = -1.0;
  }
  PetscCall(MatDenseRestoreArray(C1,&u));
  PetscCall(MatAssemblyBegin(C1,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(C1,MAT_FINAL_ASSEMBLY));
  PetscCall(MatCreateLRC(NULL,C1,NULL,NULL,&C));
  PetscCall(MatDestroy(&C1));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the solver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(LMECreate(PETSC_COMM_WORLD,&lme));
  PetscCall(LMESetProblemType(lme,LME_SYLVESTER));
  PetscCall(LMEGetProblemType(lme,&ptype));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Equation type set to %d\n",ptype));
  PetscCall(LMESetProblemType(lme,LME_LYAPUNOV));
  PetscCall(LMEGetProblemType(lme,&ptype));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Equation type changed to %d\n",ptype));
  PetscCall(LMESetCoefficients(lme,A,NULL,NULL,NULL));
  PetscCall(LMESetRHS(lme,C));

  /* test prefix usage */
  if (testprefix) {
    PetscCall(LMESetOptionsPrefix(lme,"check_"));
    PetscCall(LMEAppendOptionsPrefix(lme,"myprefix_"));
    PetscCall(LMEGetOptionsPrefix(lme,&prefix));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD," LME prefix is currently: %s\n",prefix));
  }

  /* test some interface functions */
  PetscCall(LMEGetCoefficients(lme,&B,NULL,NULL,NULL));
  if (viewmatrices) PetscCall(MatView(B,PETSC_VIEWER_STDOUT_WORLD));
  PetscCall(LMEGetRHS(lme,&D));
  if (viewmatrices) PetscCall(MatView(D,PETSC_VIEWER_STDOUT_WORLD));
  PetscCall(LMESetTolerances(lme,PETSC_CURRENT,100));
  PetscCall(LMESetDimensions(lme,21));
  PetscCall(LMESetErrorIfNotConverged(lme,PETSC_TRUE));
  /* test monitors */
  PetscCall(PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_DEFAULT,&vf));
  PetscCall(LMEMonitorSet(lme,(PetscErrorCode (*)(LME,PetscInt,PetscReal,void*))LMEMonitorDefault,vf,(PetscCtxDestroyFn*)PetscViewerAndFormatDestroy));
  /* PetscCall(LMEMonitorCancel(lme)); */
  PetscCall(LMESetFromOptions(lme));

  PetscCall(LMEGetType(lme,&type));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Solver being used: %s\n",type));

  /* query properties and print them */
  PetscCall(LMEGetTolerances(lme,&tol,&maxit));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Tolerance: %g, max iterations: %" PetscInt_FMT "\n",(double)tol,maxit));
  PetscCall(LMEGetDimensions(lme,&ncv));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Subspace dimension: %" PetscInt_FMT "\n",ncv));
  PetscCall(LMEGetErrorIfNotConverged(lme,&flg));
  if (flg) PetscCall(PetscPrintf(PETSC_COMM_WORLD," Erroring out if convergence fails\n"));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Solve the matrix equation and compute residual error
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(LMESolve(lme));
  PetscCall(LMEGetErrorEstimate(lme,&errest));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Error estimate reported by the solver: %.4g\n",(double)errest));
  PetscCall(LMEComputeError(lme,&error));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Computed residual norm: %.4g\n\n",(double)error));

  /*
     Free work space
  */
  PetscCall(LMEDestroy(&lme));
  PetscCall(MatDestroy(&A));
  PetscCall(MatDestroy(&C));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   test:
      suffix: 1
      args: -lme_monitor_cancel -lme_converged_reason -lme_view -view_matrices -log_exclude lme,bv
      requires: double
      filter: sed -e "s/4.0[0-9]*e-10/4.03e-10/"

   test:
      suffix: 2
      args: -test_prefix -check_myprefix_lme_monitor
      requires: double
      filter: sed -e "s/estimate [0-9]\.[0-9]*e[+-]\([0-9]*\)/estimate (removed)/g" | sed -e "s/4.0[0-9]*e-10/4.03e-10/"

   test:
      suffix: 3
      args: -lme_monitor_cancel -info -lme_monitor draw::draw_lg -draw_virtual
      requires: x double
      filter: sed -e "s/equation = [0-9]\.[0-9]*e[+-]\([0-9]*\)/equation = (removed)/g" | sed -e "s/4.0[0-9]*e-10/4.03e-10/" | grep -v Comm | grep -v machine | grep -v PetscGetHostName | grep -v OpenMP | grep -v Colormap | grep -v "Rank of the Cholesky factor" | grep -v "potrf failed" | grep -v "querying" | grep -v FPTrap | grep -v Device | grep -v SetNumThread

TEST*/
