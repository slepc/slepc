/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test PEP view and monitor functionality.\n\n";

#include <slepcpep.h>

int main(int argc,char **argv)
{
  Mat            A[3];
  PEP            pep;
  Vec            xr,xi;
  PetscScalar    kr,ki;
  PetscInt       n=6,Istart,Iend,i,nconv,its;
  PetscReal      errest;
  PetscBool      checkfile;
  char           filename[PETSC_MAX_PATH_LEN];
  PetscViewer    viewer;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nPEP of diagonal problem, n=%" PetscInt_FMT "\n\n",n));

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

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                     Create the PEP solver
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(PEPCreate(PETSC_COMM_WORLD,&pep));
  PetscCall(PetscObjectSetName((PetscObject)pep,"pep"));
  PetscCall(PEPSetOperators(pep,3,A));
  PetscCall(PEPSetFromOptions(pep));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Solve the eigensystem and display solution
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(PEPSolve(pep));
  PetscCall(PEPGetConverged(pep,&nconv));
  PetscCall(PEPGetIterationNumber(pep,&its));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," %" PetscInt_FMT " converged eigenpairs after %" PetscInt_FMT " iterations\n",nconv,its));
  if (nconv>0) {
    PetscCall(MatCreateVecs(A[0],&xr,&xi));
    PetscCall(PEPGetEigenpair(pep,0,&kr,&ki,xr,xi));
    PetscCall(VecDestroy(&xr));
    PetscCall(VecDestroy(&xi));
    PetscCall(PEPGetErrorEstimate(pep,0,&errest));
  }
  PetscCall(PEPErrorView(pep,PEP_ERROR_RELATIVE,NULL));

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
      PetscCall(PEPGetEigenpair(pep,i,&kr,&ki,NULL,NULL));
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

  PetscCall(PEPDestroy(&pep));
  PetscCall(MatDestroy(&A[0]));
  PetscCall(MatDestroy(&A[1]));
  PetscCall(MatDestroy(&A[2]));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   test:
      suffix: 1
      args: -pep_error_backward ::ascii_info_detail -pep_largest_real -pep_view_values -pep_monitor_conv -pep_error_absolute ::ascii_matlab -pep_monitor_all -pep_converged_reason -pep_view
      requires: !single
      filter: grep -v "tolerance" | grep -v "problem type" | sed -e "s/[+-]0\.0*i//g" -e "s/\([0-9]\.[5]*\)[+-][0-9]\.[0-9]*e-[0-9]*i/\\1/g" -e "s/[0-9]\.[0-9]*e-\([0-9]*\)/removed/g"

   test:
      suffix: 2
      args: -n 12 -pep_largest_real -pep_monitor -pep_view_values ::ascii_matlab
      requires: double
      filter: sed -e "s/[+-][0-9]\.[0-9]*e-[0-9]*i//" -e "s/[0-9]\.[0-9]*e-\([0-9]*\)/removed/g" -e "s/5\.\([49]\)999999[0-9]*e+00/5.\\1999999999999999e+00/"

   test:
      suffix: 3
      args: -pep_nev 4 -pep_view_values binary:myvalues.bin -checkfile myvalues.bin
      requires: double c99_complex

   test:
      suffix: 4
      args: -pep_nev 4 -pep_ncv 10 -pep_refine -pep_conv_norm -pep_extract none -pep_scale scalar -pep_view -pep_monitor -pep_error_relative ::ascii_info_detail
      requires: double !complex
      filter: grep -v "tolerance" | sed -e "s/[0-9]\.[0-9]*e-\([0-9]*\)/removed/g"

   test:
      suffix: 5
      args: -n 12 -pep_largest_real -pep_monitor draw::draw_lg -pep_monitor_all draw::draw_lg -pep_view_values draw -draw_save myeigen.ppm -draw_virtual
      requires: x double

TEST*/
