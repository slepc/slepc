/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test NEP view and monitor functionality.\n\n";

#include <slepcnep.h>

int main(int argc,char **argv)
{
  Mat                A[3];
  FN                 f[3];
  NEP                nep;
  Vec                xr,xi;
  PetscScalar        kr,ki,coeffs[3];
  PetscInt           n=6,i,Istart,Iend,nconv,its;
  PetscReal          errest;
  PetscBool          checkfile;
  char               filename[PETSC_MAX_PATH_LEN];
  PetscViewer        viewer;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nDiagonal Nonlinear Eigenproblem, n=%" PetscInt_FMT "\n\n",n));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Generate the matrices
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

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
  for (i=Istart;i<Iend;i++) PetscCall(MatSetValue(A[1],i,i,-1.5,INSERT_VALUES));
  PetscCall(MatAssemblyBegin(A[1],MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A[1],MAT_FINAL_ASSEMBLY));

  PetscCall(MatCreate(PETSC_COMM_WORLD,&A[2]));
  PetscCall(MatSetSizes(A[2],PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetFromOptions(A[2]));
  for (i=Istart;i<Iend;i++) PetscCall(MatSetValue(A[2],i,i,-1.0/(i+1),INSERT_VALUES));
  PetscCall(MatAssemblyBegin(A[2],MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A[2],MAT_FINAL_ASSEMBLY));

  /*
     Functions: f0=1.0, f1=lambda, f2=lambda^2
  */
  PetscCall(FNCreate(PETSC_COMM_WORLD,&f[0]));
  PetscCall(FNSetType(f[0],FNRATIONAL));
  coeffs[0] = 1.0;
  PetscCall(FNRationalSetNumerator(f[0],1,coeffs));

  PetscCall(FNCreate(PETSC_COMM_WORLD,&f[1]));
  PetscCall(FNSetType(f[1],FNRATIONAL));
  coeffs[0] = 1.0; coeffs[1] = 0.0;
  PetscCall(FNRationalSetNumerator(f[1],2,coeffs));

  PetscCall(FNCreate(PETSC_COMM_WORLD,&f[2]));
  PetscCall(FNSetType(f[2],FNRATIONAL));
  coeffs[0] = 1.0; coeffs[1] = 0.0; coeffs[2] = 0.0;
  PetscCall(FNRationalSetNumerator(f[2],3,coeffs));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Create the NEP solver
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(NEPCreate(PETSC_COMM_WORLD,&nep));
  PetscCall(PetscObjectSetName((PetscObject)nep,"nep"));
  PetscCall(NEPSetSplitOperator(nep,3,A,f,SAME_NONZERO_PATTERN));
  PetscCall(NEPSetTarget(nep,1.1));
  PetscCall(NEPSetWhichEigenpairs(nep,NEP_TARGET_MAGNITUDE));
  PetscCall(NEPSetFromOptions(nep));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Solve the eigensystem and display solution
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(NEPSolve(nep));
  PetscCall(NEPGetConverged(nep,&nconv));
  PetscCall(NEPGetIterationNumber(nep,&its));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," %" PetscInt_FMT " converged eigenpairs after %" PetscInt_FMT " iterations\n",nconv,its));
  if (nconv>0) {
    PetscCall(MatCreateVecs(A[0],&xr,&xi));
    PetscCall(NEPGetEigenpair(nep,0,&kr,&ki,xr,xi));
    PetscCall(VecDestroy(&xr));
    PetscCall(VecDestroy(&xi));
    PetscCall(NEPGetErrorEstimate(nep,0,&errest));
  }
  PetscCall(NEPErrorView(nep,NEP_ERROR_BACKWARD,NULL));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                   Check file containing the eigenvalues
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(PetscOptionsGetString(NULL,NULL,"-checkfile",filename,sizeof(filename),&checkfile));
  if (checkfile) {
#if PetscDefined(HAVE_COMPLEX)
    PetscComplex *eigs,eval;
    PetscCall(PetscMalloc1(nconv,&eigs));
    PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_READ,&viewer));
    PetscCall(PetscViewerBinaryRead(viewer,eigs,nconv,NULL,PETSC_COMPLEX));
    PetscCall(PetscViewerDestroy(&viewer));
    for (i=0;i<nconv;i++) {
      PetscCall(NEPGetEigenpair(nep,i,&kr,&ki,NULL,NULL));
#if PetscDefined(USE_COMPLEX)
      eval = kr;
#else
      eval = PetscCMPLX(kr,ki);
#endif
      PetscCheck(eval==eigs[i],PETSC_COMM_WORLD,PETSC_ERR_FILE_UNEXPECTED,"Eigenvalues in the file do not match");
    }
    PetscCall(PetscFree(eigs));
#else
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"The -checkfile option requires C99 complex numbers");
#endif
  }

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
      args: -nep_type slp -nep_target -.5 -nep_error_backward ::ascii_info_detail -nep_view_values -nep_error_absolute ::ascii_matlab -nep_monitor_all -nep_converged_reason -nep_view
      filter: grep -v "tolerance" | grep -v "problem type" | sed -e "s/[+-]0\.0*i//g" -e "s/+0i//" -e "s/[+-][0-9]\.[0-9]*e-[0-9]*i//g" -e "s/[0-9]\.[0-9]*e[+-]\([0-9]*\)/removed/g"
      requires: double

   test:
      suffix: 2
      args: -nep_type rii -nep_target -.5 -nep_rii_hermitian -nep_monitor -nep_view_values ::ascii_matlab
      filter: sed -e "s/[+-][0-9]\.[0-9]*e-[0-9]*i//" -e "s/([0-9]\.[0-9]*e[+-]\([0-9]*\))/(removed)/g"
      requires: double

   test:
      suffix: 3
      args: -nep_type slp -nep_nev 4 -nep_view_values binary:myvalues.bin -checkfile myvalues.bin -nep_error_relative ::ascii_matlab
      filter: sed -e "s/[0-9]\.[0-9]*e[+-]\([0-9]*\)/removed/g"
      requires: double c99_complex

   test:
      suffix: 4
      args: -nep_type slp -nep_nev 4 -nep_monitor draw::draw_lg -nep_monitor_all draw::draw_lg -nep_view_values draw -draw_save myeigen.ppm -draw_virtual
      requires: x double

TEST*/
