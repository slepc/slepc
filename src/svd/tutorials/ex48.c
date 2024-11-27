/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Solves a GSVD problem with matrices loaded from a file.\n"
  "The command line options are:\n"
  "  -f1 <filename>, PETSc binary file containing matrix A.\n"
  "  -f2 <filename>, PETSc binary file containing matrix B (optional). Instead of"
  "     a file it is possible to specify one of 'identity', 'bidiagonal' or 'tridiagonal'"
  "  -p <p>, in case B is not taken from a file.\n\n";

#include <slepcsvd.h>

int main(int argc,char **argv)
{
  Mat            A,B;             /* matrices */
  SVD            svd;             /* singular value problem solver context */
  PetscInt       i,m,n,p,Istart,Iend,col[3];
  PetscScalar    vals[3];
  char           filename[PETSC_MAX_PATH_LEN];
  PetscViewer    viewer;
  PetscBool      flg,terse;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Load matrices that define the generalized singular value problem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nGeneralized singular value problem stored in file.\n\n"));
  PetscCall(PetscOptionsGetString(NULL,NULL,"-f1",filename,sizeof(filename),&flg));
  PetscCheck(flg,PETSC_COMM_WORLD,PETSC_ERR_USER_INPUT,"Must indicate a file name for matrix A with the -f1 option");

#if defined(PETSC_USE_COMPLEX)
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Reading COMPLEX matrix from a binary file...\n"));
#else
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Reading REAL matrix from a binary file...\n"));
#endif
  PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_READ,&viewer));
  PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
  PetscCall(MatSetFromOptions(A));
  PetscCall(MatLoad(A,viewer));
  PetscCall(PetscViewerDestroy(&viewer));

  PetscCall(MatGetSize(A,&m,&n));

  PetscCall(PetscOptionsGetString(NULL,NULL,"-f2",filename,sizeof(filename),&flg));
  PetscCheck(flg,PETSC_COMM_WORLD,PETSC_ERR_USER_INPUT,"Must indicate a file name for matrix B with the -f2 option, or alternatively the strings 'identity', 'bidiagonal', or 'tridiagonal'");
  PetscCall(PetscStrcmp(filename,"identity",&flg));
  if (flg) {
    p = n;
    PetscCall(PetscOptionsGetInt(NULL,NULL,"-p",&p,&flg));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD," Using B=I with p=%" PetscInt_FMT "\n\n",p));
    PetscCall(MatCreate(PETSC_COMM_WORLD,&B));
    PetscCall(MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,p,n));
    PetscCall(MatSetFromOptions(B));
    PetscCall(MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY));
    PetscCall(MatShift(B,1.0));
  } else {
    PetscCall(PetscStrcmp(filename,"bidiagonal",&flg));
    if (flg) {
      p = n+1;
      PetscCall(PetscOptionsGetInt(NULL,NULL,"-p",&p,&flg));
      vals[0]=-1; vals[1]=1;
      PetscCall(PetscPrintf(PETSC_COMM_WORLD," Using B=bidiag(1,-1) with p=%" PetscInt_FMT "\n\n",p));
      PetscCall(MatCreate(PETSC_COMM_WORLD,&B));
      PetscCall(MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,p,n));
      PetscCall(MatSetFromOptions(B));
      PetscCall(MatGetOwnershipRange(B,&Istart,&Iend));
      for (i=Istart;i<Iend;i++) {
        col[0]=i-1; col[1]=i;
        if (i==0) PetscCall(MatSetValue(B,i,col[1],vals[1],INSERT_VALUES));
        else if (i<n) PetscCall(MatSetValues(B,1,&i,2,col,vals,INSERT_VALUES));
        else if (i==n) PetscCall(MatSetValue(B,i,col[0],vals[0],INSERT_VALUES));
      }
      PetscCall(MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY));
      PetscCall(MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY));
    } else {
      PetscCall(PetscStrcmp(filename,"tridiagonal",&flg));
      if (flg) {
        p = n-2;
        PetscCall(PetscOptionsGetInt(NULL,NULL,"-p",&p,&flg));
        vals[0]=-1; vals[1]=2; vals[2]=-1;
        PetscCall(PetscPrintf(PETSC_COMM_WORLD," Using B=tridiag(-1,2,-1) with p=%" PetscInt_FMT "\n\n",p));
        PetscCall(MatCreate(PETSC_COMM_WORLD,&B));
        PetscCall(MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,p,n));
        PetscCall(MatSetFromOptions(B));
        PetscCall(MatGetOwnershipRange(B,&Istart,&Iend));
        for (i=Istart;i<Iend;i++) {
          col[0]=i; col[1]=i+1; col[2]=i+2;
          PetscCall(MatSetValues(B,1,&i,3,col,vals,INSERT_VALUES));
        }
        PetscCall(MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY));
        PetscCall(MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY));
      } else {  /* load file */
        PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_READ,&viewer));
        PetscCall(MatCreate(PETSC_COMM_WORLD,&B));
        PetscCall(MatSetFromOptions(B));
        PetscCall(MatLoad(B,viewer));
        PetscCall(PetscViewerDestroy(&viewer));
      }
    }
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
           Create the singular value solver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Create singular value solver context
  */
  PetscCall(SVDCreate(PETSC_COMM_WORLD,&svd));

  /*
     Set operators of GSVD problem
  */
  PetscCall(SVDSetOperators(svd,A,B));
  PetscCall(SVDSetProblemType(svd,SVD_GENERALIZED));

  /*
     Set solver parameters at runtime
  */
  PetscCall(SVDSetFromOptions(svd));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the problem and print solution
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(SVDSolve(svd));

  /* show detailed info unless -terse option is given by user */
  PetscCall(PetscOptionsHasName(NULL,NULL,"-terse",&terse));
  if (terse) PetscCall(SVDErrorView(svd,SVD_ERROR_NORM,NULL));
  else {
    PetscCall(PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL));
    PetscCall(SVDConvergedReasonView(svd,PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(SVDErrorView(svd,SVD_ERROR_NORM,PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD));
  }
  PetscCall(SVDDestroy(&svd));
  PetscCall(MatDestroy(&A));
  PetscCall(MatDestroy(&B));
  PetscCall(SlepcFinalize());
  return 0;
}
/*TEST

   testset:
      requires: double !complex !defined(PETSC_USE_64BIT_INDICES)
      args: -f1 ${SLEPC_DIR}/share/slepc/datafiles/matrices/bfw62a.petsc -f2 ${SLEPC_DIR}/share/slepc/datafiles/matrices/bfw62b.petsc -svd_nsv 3 -terse
      output_file: output/ex48_1.out
      test:
         suffix: 1
         args: -svd_type trlanczos -svd_trlanczos_explicitmatrix {{0 1}} -svd_trlanczos_scale 1e5 -svd_trlanczos_ksp_rtol 1e-15
      test:
         suffix: 1_spqr
         args: -svd_type trlanczos -svd_trlanczos_explicitmatrix -svd_trlanczos_pc_type qr -svd_trlanczos_scale 1e5 -svd_trlanczos_oneside {{0 1}}
         requires: suitesparse
      test:
         suffix: 1_autoscale
         args: -svd_type trlanczos -svd_trlanczos_gbidiag {{lower upper}} -svd_trlanczos_scale -5 -svd_trlanczos_ksp_rtol 1e-16 -svd_trlanczos_oneside {{0 1}}
      test:
         suffix: 1_cross
         args: -svd_type cross -svd_cross_explicitmatrix
      test:
         suffix: 1_cyclic
         args: -svd_type cyclic -svd_cyclic_explicitmatrix

   testset:
      requires: double complex datafilespath !defined(PETSC_USE_64BIT_INDICES)
      args: -f1 ${DATAFILESPATH}/matrices/complex/qc324.petsc -f2 bidiagonal -svd_nsv 3 -terse
      output_file: output/ex48_2.out
      filter: sed -e "s/30749/30748/"
      timeoutfactor: 2
      test:
         suffix: 2
         args: -svd_type trlanczos -svd_trlanczos_explicitmatrix {{0 1}} -svd_trlanczos_ksp_rtol 1e-10 -svd_trlanczos_scale 100
         requires: !defined(PETSCTEST_VALGRIND)
      test:
         suffix: 2_spqr
         args: -svd_type trlanczos -svd_trlanczos_explicitmatrix -svd_trlanczos_pc_type qr -svd_trlanczos_ksp_rtol 1e-10
         requires: suitesparse
      test:
         suffix: 2_cross
         args: -svd_type cross -svd_cross_explicitmatrix
      test:
         suffix: 2_cyclic
         args: -svd_type cyclic -svd_cyclic_explicitmatrix

   test:
      requires: double complex datafilespath !defined(PETSC_USE_64BIT_INDICES) !defined(PETSCTEST_VALGRIND)
      args: -f1 ${DATAFILESPATH}/matrices/complex/qc324.petsc -f2 bidiagonal -p 320 -svd_nsv 3 -svd_type trlanczos -svd_trlanczos_ksp_rtol 1e-14 -svd_trlanczos_scale 100 -terse
      timeoutfactor: 2
      suffix: 3

   testset:
      requires: double !complex !defined(PETSC_USE_64BIT_INDICES)
      args: -f1 ${SLEPC_DIR}/share/slepc/datafiles/matrices/rdb200.petsc -f2 identity -svd_nsv 3 -svd_ncv 24 -svd_smallest -terse
      output_file: output/ex48_4.out
      test:
         suffix: 4
         args: -svd_type trlanczos
      test:
         suffix: 4_spqr
         args: -svd_type trlanczos -svd_trlanczos_explicitmatrix -svd_trlanczos_pc_type qr
         requires: suitesparse
      test:
         suffix: 4_cross
         args: -svd_type cross -svd_cross_explicitmatrix
      test:
         suffix: 4_cross_implicit
         args: -svd_type cross -svd_cross_eps_type lobpcg -svd_cross_st_ksp_type cg -svd_cross_st_pc_type jacobi -svd_max_it 1000
      test:
         suffix: 4_cyclic
         args: -svd_type cyclic -svd_cyclic_explicitmatrix
      test:
         suffix: 4_hpddm
         nsize: 4
         args: -svd_type trlanczos -svd_trlanczos_explicitmatrix -svd_trlanczos_pc_type hpddm
         args: -prefix_push svd_trlanczos_pc_hpddm_ -levels_1_st_share_sub_ksp -levels_1_eps_nev 10 -levels_1_eps_threshold 0.005 -levels_1_pc_asm_type basic -define_subdomains -levels_1_pc_asm_sub_mat_type sbaij -levels_1_sub_pc_type cholesky -prefix_pop
         requires: hpddm

   testset:
      args: -f1 ${SLEPC_DIR}/share/slepc/datafiles/matrices/rdb200.petsc -f2 identity -svd_threshold_relative 0.87 -terse
      output_file: output/ex48_5.out
      requires: double !complex !defined(PETSC_USE_64BIT_INDICES)
      test:
         suffix: 5
         args: -svd_type trlanczos -svd_trlanczos_explicitmatrix -svd_trlanczos_scale 100 -svd_trlanczos_gbidiag {{lower upper single}}
      test:
         suffix: 5_cross
         args: -svd_type cross -svd_cross_explicitmatrix

   testset:
      args: -f1 ${DATAFILESPATH}/matrices/complex/qc324.petsc -f2 bidiagonal -svd_threshold_relative 0.3 -svd_ncv 8 -svd_tol 1e-7 -terse
      output_file: output/ex48_5_complex.out
      requires: double complex datafilespath !defined(PETSC_USE_64BIT_INDICES)
      test:
         suffix: 5_complex
         args: -svd_type trlanczos -svd_trlanczos_explicitmatrix -svd_trlanczos_scale 100 -svd_trlanczos_ksp_rtol 1e-6 -svd_trlanczos_gbidiag {{lower upper single}}
      test:
         suffix: 5_complex_cross
         args: -svd_type cross -svd_cross_explicitmatrix

TEST*/
