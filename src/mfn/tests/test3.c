/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test MFN interface functions.\n\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = number of grid subdivisions in x dimension.\n\n";

#include <slepcmfn.h>

int main(int argc,char **argv)
{
  Mat                  A,B;
  MFN                  mfn;
  FN                   f;
  MFNConvergedReason   reason;
  MFNType              type;
  PetscReal            norm,tol;
  Vec                  v,y;
  PetscInt             N,n=4,Istart,Iend,i,j,II,ncv,its,maxit;
  PetscBool            flg,testprefix=PETSC_FALSE;
  const char           *prefix;
  PetscViewerAndFormat *vf;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  N = n*n;
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nSquare root of Laplacian y=sqrt(A)*e_1, N=%" PetscInt_FMT " (%" PetscInt_FMT "x%" PetscInt_FMT " grid)\n\n",N,n,n));
  PetscCall(PetscOptionsGetBool(NULL,NULL,"-test_prefix",&testprefix,NULL));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                 Compute the discrete 2-D Laplacian, A
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
  PetscCall(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,N,N));
  PetscCall(MatSetFromOptions(A));

  PetscCall(MatGetOwnershipRange(A,&Istart,&Iend));
  for (II=Istart;II<Iend;II++) {
    i = II/n; j = II-i*n;
    if (i>0) PetscCall(MatSetValue(A,II,II-n,-1.0,INSERT_VALUES));
    if (i<n-1) PetscCall(MatSetValue(A,II,II+n,-1.0,INSERT_VALUES));
    if (j>0) PetscCall(MatSetValue(A,II,II-1,-1.0,INSERT_VALUES));
    if (j<n-1) PetscCall(MatSetValue(A,II,II+1,-1.0,INSERT_VALUES));
    PetscCall(MatSetValue(A,II,II,4.0,INSERT_VALUES));
  }

  PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));
  PetscCall(MatSetOption(A,MAT_HERMITIAN,PETSC_TRUE));

  PetscCall(MatCreateVecs(A,&v,&y));
  PetscCall(VecSetValue(v,0,1.0,INSERT_VALUES));
  PetscCall(VecAssemblyBegin(v));
  PetscCall(VecAssemblyEnd(v));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
             Create the solver, set the matrix and the function
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(MFNCreate(PETSC_COMM_WORLD,&mfn));
  PetscCall(MFNSetOperator(mfn,A));
  PetscCall(MFNGetFN(mfn,&f));
  PetscCall(FNSetType(f,FNSQRT));

  PetscCall(MFNSetType(mfn,MFNKRYLOV));
  PetscCall(MFNGetType(mfn,&type));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Type set to %s\n",type));

  /* test prefix usage */
  if (testprefix) {
    PetscCall(MFNSetOptionsPrefix(mfn,"check_"));
    PetscCall(MFNAppendOptionsPrefix(mfn,"myprefix_"));
    PetscCall(MFNGetOptionsPrefix(mfn,&prefix));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD," MFN prefix is currently: %s\n",prefix));
  }

  /* test some interface functions */
  PetscCall(MFNGetOperator(mfn,&B));
  PetscCall(MatView(B,PETSC_VIEWER_STDOUT_WORLD));
  PetscCall(MFNSetTolerances(mfn,1e-4,500));
  PetscCall(MFNSetDimensions(mfn,6));
  PetscCall(MFNSetErrorIfNotConverged(mfn,PETSC_TRUE));
  /* test monitors */
  PetscCall(PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_DEFAULT,&vf));
  PetscCall(MFNMonitorSet(mfn,(PetscErrorCode (*)(MFN,PetscInt,PetscReal,void*))MFNMonitorDefault,vf,(PetscCtxDestroyFn*)PetscViewerAndFormatDestroy));
  /* PetscCall(MFNMonitorCancel(mfn)); */
  PetscCall(MFNSetFromOptions(mfn));

  /* query properties and print them */
  PetscCall(MFNGetTolerances(mfn,&tol,&maxit));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Tolerance: %g, max iterations: %" PetscInt_FMT "\n",(double)tol,maxit));
  PetscCall(MFNGetDimensions(mfn,&ncv));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Subspace dimension: %" PetscInt_FMT "\n",ncv));
  PetscCall(MFNGetErrorIfNotConverged(mfn,&flg));
  if (flg) PetscCall(PetscPrintf(PETSC_COMM_WORLD," Erroring out if convergence fails\n"));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                           Solve  y=sqrt(A)*v
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(MFNSolve(mfn,v,y));
  PetscCall(MFNGetConvergedReason(mfn,&reason));
  PetscCall(MFNGetIterationNumber(mfn,&its));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Finished - converged reason = %d\n",(int)reason));
  /* PetscCall(PetscPrintf(PETSC_COMM_WORLD," its = %" PetscInt_FMT "\n",its)); */
  PetscCall(VecNorm(y,NORM_2,&norm));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," sqrt(A)*v has norm %g\n",(double)norm));

  /*
     Free work space
  */
  PetscCall(MFNDestroy(&mfn));
  PetscCall(MatDestroy(&A));
  PetscCall(VecDestroy(&v));
  PetscCall(VecDestroy(&y));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   testset:
      args: -mfn_monitor_cancel -mfn_converged_reason -mfn_view -log_exclude mfn,bv,fn
      output_file: output/test3_1.out
      test:
         suffix: 1
      test:
         suffix: 1_x
         args: -mfn_monitor draw::draw_lg -draw_virtual
         requires: x

   test:
      suffix: 2
      args: -test_prefix -check_myprefix_mfn_monitor
      filter: sed -e "s/estimate [0-9]\.[0-9]*e[+-]\([0-9]*\)/estimate (removed)/g" | sed -e "s/4.0[0-9]*e-10/4.03e-10/"

TEST*/
