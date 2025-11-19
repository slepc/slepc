/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Example based on spring problem in NLEVP collection [1]. See the parameters
   meaning at Example 2 in [2].

   [1] T. Betcke, N. J. Higham, V. Mehrmann, C. Schroder, and F. Tisseur,
       NLEVP: A Collection of Nonlinear Eigenvalue Problems, MIMS EPrint
       2010.98, November 2010.
   [2] F. Tisseur, Backward error and condition of polynomial eigenvalue
       problems, Linear Algebra and its Applications, 309 (2000), pp. 339--361,
       April 2000.
*/

static char help[] = "Test the solution of a PEP from a finite element model of damped mass-spring system.\n\n"
  "Problem from NLEVP collection.\n"
  "The command line options are:\n"
  "  -n <n> ... number of grid subdivisions.\n"
  "  -mu <value> ... mass (default 1).\n"
  "  -tau <value> ... damping constant of the dampers (default 10).\n"
  "  -kappa <value> ... damping constant of the springs (default 5).\n"
  "  -initv ... set an initial vector.\n\n";

#include <slepcpep.h>

/*
   Check if computed eigenvectors have unit norm
*/
PetscErrorCode CheckNormalizedVectors(PEP pep)
{
  PetscInt       i,nconv;
  Mat            A;
  Vec            xr,xi;
  PetscReal      error=0.0,normr;
#if !defined(PETSC_USE_COMPLEX)
  PetscReal      normi;
#endif

  PetscFunctionBeginUser;
  PetscCall(PEPGetConverged(pep,&nconv));
  if (nconv>0) {
    PetscCall(PEPGetOperators(pep,0,&A));
    PetscCall(MatCreateVecs(A,&xr,&xi));
    for (i=0;i<nconv;i++) {
      PetscCall(PEPGetEigenpair(pep,i,NULL,NULL,xr,xi));
#if defined(PETSC_USE_COMPLEX)
      PetscCall(VecNorm(xr,NORM_2,&normr));
      error = PetscMax(error,PetscAbsReal(normr-PetscRealConstant(1.0)));
#else
      PetscCall(VecNormBegin(xr,NORM_2,&normr));
      PetscCall(VecNormBegin(xi,NORM_2,&normi));
      PetscCall(VecNormEnd(xr,NORM_2,&normr));
      PetscCall(VecNormEnd(xi,NORM_2,&normi));
      error = PetscMax(error,PetscAbsReal(SlepcAbsEigenvalue(normr,normi)-PetscRealConstant(1.0)));
#endif
    }
    PetscCall(VecDestroy(&xr));
    PetscCall(VecDestroy(&xi));
    if (error>100*PETSC_MACHINE_EPSILON) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Vectors are not normalized. Error=%g\n",(double)error));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

int main(int argc,char **argv)
{
  Mat            M,C,K,A[3];      /* problem matrices */
  PEP            pep;             /* polynomial eigenproblem solver context */
  PetscInt       n=30,Istart,Iend,i,nev;
  PetscReal      mu=1.0,tau=10.0,kappa=5.0;
  PetscBool      initv=PETSC_FALSE,skipnorm=PETSC_FALSE;
  Vec            IV[2];

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));

  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-mu",&mu,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-tau",&tau,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-kappa",&kappa,NULL));
  PetscCall(PetscOptionsGetBool(NULL,NULL,"-initv",&initv,NULL));
  PetscCall(PetscOptionsGetBool(NULL,NULL,"-skipnorm",&skipnorm,NULL));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the matrices that define the eigensystem, (k^2*M+k*C+K)x=0
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* K is a tridiagonal */
  PetscCall(MatCreate(PETSC_COMM_WORLD,&K));
  PetscCall(MatSetSizes(K,PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetFromOptions(K));

  PetscCall(MatGetOwnershipRange(K,&Istart,&Iend));
  for (i=Istart;i<Iend;i++) {
    if (i>0) PetscCall(MatSetValue(K,i,i-1,-kappa,INSERT_VALUES));
    PetscCall(MatSetValue(K,i,i,kappa*3.0,INSERT_VALUES));
    if (i<n-1) PetscCall(MatSetValue(K,i,i+1,-kappa,INSERT_VALUES));
  }

  PetscCall(MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY));

  /* C is a tridiagonal */
  PetscCall(MatCreate(PETSC_COMM_WORLD,&C));
  PetscCall(MatSetSizes(C,PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetFromOptions(C));

  PetscCall(MatGetOwnershipRange(C,&Istart,&Iend));
  for (i=Istart;i<Iend;i++) {
    if (i>0) PetscCall(MatSetValue(C,i,i-1,-tau,INSERT_VALUES));
    PetscCall(MatSetValue(C,i,i,tau*3.0,INSERT_VALUES));
    if (i<n-1) PetscCall(MatSetValue(C,i,i+1,-tau,INSERT_VALUES));
  }

  PetscCall(MatAssemblyBegin(C,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(C,MAT_FINAL_ASSEMBLY));

  /* M is a diagonal matrix */
  PetscCall(MatCreate(PETSC_COMM_WORLD,&M));
  PetscCall(MatSetSizes(M,PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetFromOptions(M));
  PetscCall(MatGetOwnershipRange(M,&Istart,&Iend));
  for (i=Istart;i<Iend;i++) PetscCall(MatSetValue(M,i,i,mu,INSERT_VALUES));
  PetscCall(MatAssemblyBegin(M,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(M,MAT_FINAL_ASSEMBLY));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(PEPCreate(PETSC_COMM_WORLD,&pep));
  A[0] = K; A[1] = C; A[2] = M;
  PetscCall(PEPSetOperators(pep,3,A));
  PetscCall(PEPSetProblemType(pep,PEP_GENERAL));
  PetscCall(PEPSetTolerances(pep,PETSC_SMALL,PETSC_CURRENT));
  if (initv) { /* initial vector */
    PetscCall(MatCreateVecs(K,&IV[0],NULL));
    PetscCall(VecSetValue(IV[0],0,-1.0,INSERT_VALUES));
    PetscCall(VecSetValue(IV[0],1,0.5,INSERT_VALUES));
    PetscCall(VecAssemblyBegin(IV[0]));
    PetscCall(VecAssemblyEnd(IV[0]));
    PetscCall(MatCreateVecs(K,&IV[1],NULL));
    PetscCall(VecSetValue(IV[1],0,4.0,INSERT_VALUES));
    PetscCall(VecSetValue(IV[1],2,1.5,INSERT_VALUES));
    PetscCall(VecAssemblyBegin(IV[1]));
    PetscCall(VecAssemblyEnd(IV[1]));
    PetscCall(PEPSetInitialSpace(pep,2,IV));
    PetscCall(VecDestroy(&IV[0]));
    PetscCall(VecDestroy(&IV[1]));
  }
  PetscCall(PEPSetFromOptions(pep));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(PEPSolve(pep));
  PetscCall(PEPGetDimensions(pep,&nev,NULL,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %" PetscInt_FMT "\n",nev));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(PEPErrorView(pep,PEP_ERROR_BACKWARD,NULL));
  if (!skipnorm) PetscCall(CheckNormalizedVectors(pep));
  PetscCall(PEPDestroy(&pep));
  PetscCall(MatDestroy(&M));
  PetscCall(MatDestroy(&C));
  PetscCall(MatDestroy(&K));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   testset:
      args: -pep_nev 4 -initv
      requires: !single
      output_file: output/test2_1.out
      test:
         suffix: 1
         args: -pep_type {{toar linear}}
      test:
         suffix: 1_toar_mgs
         args: -pep_type toar -bv_orthog_type mgs
      test:
         suffix: 1_qarnoldi
         args: -pep_type qarnoldi -bv_orthog_refine never
      test:
         suffix: 1_linear_gd
         args: -pep_type linear -pep_linear_eps_type gd -pep_linear_explicitmatrix

   testset:
      args: -pep_target -0.43 -pep_nev 4 -pep_ncv 20 -st_type sinvert
      output_file: output/test2_2.out
      test:
         suffix: 2
         args: -pep_type {{toar linear}}
      test:
         suffix: 2_toar_scaleboth
         args: -pep_type toar -pep_scale both
      test:
         suffix: 2_toar_transform
         args: -pep_type toar -st_transform
      test:
         suffix: 2_qarnoldi
         args: -pep_type qarnoldi -bv_orthog_refine always
      test:
         suffix: 2_linear_explicit
         args: -pep_type linear -pep_linear_explicitmatrix -pep_linear_linearization 0,1
      test:
         suffix: 2_linear_explicit_her
         args: -pep_type linear -pep_linear_explicitmatrix -pep_hermitian -pep_linear_linearization {{0,1 1,0 .3,.7}}
      test:
         suffix: 2_stoar
         args: -pep_type stoar -pep_hermitian
         requires: !single
      test:
         suffix: 2_jd
         args: -pep_type jd -st_type precond -pep_max_it 200 -pep_ncv 24
         requires: !single

   test:
      suffix: 3
      args: -pep_nev 12 -pep_extract {{none norm residual structured}} -pep_monitor_cancel
      requires: !single

   testset:
      args: -st_type sinvert -pep_target -0.43 -pep_nev 4
      output_file: output/test2_2.out
      test:
         suffix: 4_schur
         args: -pep_refine simple -pep_refine_scheme schur
      test:
         suffix: 4_mbe
         args: -pep_refine simple -pep_refine_scheme mbe -pep_refine_ksp_type preonly -pep_refine_pc_type lu
      test:
         suffix: 4_explicit
         args: -pep_refine simple -pep_refine_scheme explicit
      test:
         suffix: 4_multiple_schur
         args: -pep_refine multiple -pep_refine_scheme schur
         requires: !single
      test:
         suffix: 4_multiple_mbe
         args: -pep_refine multiple -pep_refine_scheme mbe -pep_refine_ksp_type preonly -pep_refine_pc_type lu -pep_refine_pc_factor_shift_type nonzero
      test:
         suffix: 4_multiple_explicit
         args: -pep_refine multiple -pep_refine_scheme explicit
         requires: !single

   test:
      suffix: 5
      nsize: 2
      args: -pep_type linear -pep_linear_explicitmatrix -pep_general -pep_target -0.43 -pep_nev 4 -pep_ncv 20 -st_type sinvert -pep_linear_st_ksp_type bcgs -pep_linear_st_pc_type bjacobi
      output_file: output/test2_2.out

   test:
      suffix: 6
      args: -pep_type linear -pep_nev 12 -pep_extract {{none norm residual}}
      requires: !single
      output_file: output/test2_3.out

   test:
      suffix: 7
      args: -pep_nev 12 -pep_extract {{none norm residual structured}} -pep_refine multiple
      requires: !single
      output_file: output/test2_3.out

   testset:
      args: -st_type sinvert -pep_target -0.43 -pep_nev 4 -st_transform
      output_file: output/test2_2.out
      test:
         suffix: 8_schur
         args: -pep_refine simple -pep_refine_scheme schur
      test:
         suffix: 8_mbe
         args: -pep_refine simple -pep_refine_scheme mbe -pep_refine_ksp_type preonly -pep_refine_pc_type lu
      test:
         suffix: 8_explicit
         args: -pep_refine simple -pep_refine_scheme explicit
      test:
         suffix: 8_multiple_schur
         args: -pep_refine multiple -pep_refine_scheme schur
      test:
         suffix: 8_multiple_mbe
         args: -pep_refine multiple -pep_refine_scheme mbe -pep_refine_ksp_type preonly -pep_refine_pc_type lu
      test:
         suffix: 8_multiple_explicit
         args: -pep_refine multiple -pep_refine_scheme explicit

   testset:
      nsize: 2
      args: -st_type sinvert -pep_target -0.49 -pep_nev 4 -pep_refine_partitions 2 -st_ksp_type bcgs -st_pc_type bjacobi -pep_scale diagonal -pep_scale_its 4
      output_file: output/test2_2.out
      test:
         suffix: 9_mbe
         args: -pep_refine simple -pep_refine_scheme mbe -pep_refine_ksp_type preonly -pep_refine_pc_type lu
      test:
         suffix: 9_explicit
         args: -pep_refine simple -pep_refine_scheme explicit
      test:
         suffix: 9_multiple_mbe
         args: -pep_refine multiple -pep_refine_scheme mbe -pep_refine_ksp_type preonly -pep_refine_pc_type lu
         requires: !single
      test:
         suffix: 9_multiple_explicit
         args: -pep_refine multiple -pep_refine_scheme explicit
         requires: !single

   test:
      suffix: 10
      nsize: 4
      args: -st_type sinvert -pep_target -0.43 -pep_nev 4 -pep_refine simple -pep_refine_scheme explicit -pep_refine_partitions 2 -st_ksp_type bcgs -st_pc_type bjacobi -pep_scale diagonal -pep_scale_its 4
      output_file: output/test2_2.out

   testset:
      args: -pep_nev 4 -initv -mat_type aijcusparse
      output_file: output/test2_1.out
      requires: cuda !single
      test:
         suffix: 11_cuda
         args: -pep_type {{toar linear}}
      test:
         suffix: 11_cuda_qarnoldi
         args: -pep_type qarnoldi -bv_orthog_refine never
      test:
         suffix: 11_cuda_linear_gd
         args: -pep_type linear -pep_linear_eps_type gd -pep_linear_explicitmatrix

   test:
      suffix: 12
      nsize: 2
      args: -pep_type jd -ds_parallel synchronized -pep_target -0.43 -pep_nev 4 -pep_ncv 20
      requires: !single

   test:
      suffix: 13
      args: -pep_nev 12 -pep_view_values draw -pep_monitor draw::draw_lg
      requires: x !single
      output_file: output/test2_3.out

   test:
      suffix: 14
      requires: complex double
      args: -pep_type ciss -rg_type ellipse -rg_ellipse_center -48.5 -rg_ellipse_radius 1.5 -pep_ciss_delta 1e-10

   testset:
      args: -pep_nev 4 -initv -mat_type aijhipsparse
      output_file: output/test2_1.out
      requires: hip !single
      test:
         suffix: 15_hip
         args: -pep_type {{toar linear}}
      test:
         suffix: 15_hip_qarnoldi
         args: -pep_type qarnoldi -bv_orthog_refine never
      test:
         suffix: 15_hip_linear_gd
         args: -pep_type linear -pep_linear_eps_type gd -pep_linear_explicitmatrix

TEST*/
