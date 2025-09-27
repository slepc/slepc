/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Eigenvalue problem associated with a Markov model of a random walk on a triangular grid.\n\n"
  "It is a standard nonsymmetric eigenproblem with real eigenvalues and the rightmost eigenvalue is known to be 1.\n"
  "This example illustrates how the user can set the initial vector.\n\n"
  "The command line options are:\n"
  "  -m <m>, where <m> = number of grid subdivisions in each dimension.\n\n";

#include <slepceps.h>

/*
   User-defined routines
*/
PetscErrorCode MatMarkovModel(PetscInt m,Mat A);
PetscErrorCode MyEigenSort(PetscScalar ar,PetscScalar ai,PetscScalar br,PetscScalar bi,PetscInt *r,void *ctx);

/*
   Check if computed eigenvectors have unit norm
*/
PetscErrorCode CheckNormalizedVectors(EPS eps)
{
  PetscInt       i,nconv;
  Mat            A;
  Vec            xr,xi;
  PetscReal      error=0.0,normr;
#if !defined(PETSC_USE_COMPLEX)
  PetscReal      normi;
#endif

  PetscFunctionBeginUser;
  PetscCall(EPSGetConverged(eps,&nconv));
  if (nconv>0) {
    PetscCall(EPSGetOperators(eps,&A,NULL));
    PetscCall(MatCreateVecs(A,&xr,&xi));
    for (i=0;i<nconv;i++) {
      PetscCall(EPSGetEigenvector(eps,i,xr,xi));
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
  Vec            v0;              /* initial vector */
  Mat            A;               /* operator matrix */
  EPS            eps;             /* eigenproblem solver context */
  PetscReal      tol=0.5*PETSC_SMALL;
  PetscInt       N,m=15,nev;
  PetscScalar    origin=0.0;
  PetscBool      flg,delay,skipnorm=PETSC_FALSE;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));

  PetscCall(PetscOptionsGetInt(NULL,NULL,"-m",&m,NULL));
  N = m*(m+1)/2;
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nMarkov Model, N=%" PetscInt_FMT " (m=%" PetscInt_FMT ")\n\n",N,m));
  PetscCall(PetscOptionsGetBool(NULL,NULL,"-skipnorm",&skipnorm,NULL));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the operator matrix that defines the eigensystem, Ax=kx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
  PetscCall(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,N,N));
  PetscCall(MatSetFromOptions(A));
  PetscCall(MatMarkovModel(m,A));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Create eigensolver context
  */
  PetscCall(EPSCreate(PETSC_COMM_WORLD,&eps));

  /*
     Set operators. In this case, it is a standard eigenvalue problem
  */
  PetscCall(EPSSetOperators(eps,A,NULL));
  PetscCall(EPSSetProblemType(eps,EPS_NHEP));
  PetscCall(EPSSetTolerances(eps,tol,PETSC_CURRENT));

  /*
     Set the custom comparing routine in order to obtain the eigenvalues
     closest to the target on the right only
  */
  PetscCall(EPSSetEigenvalueComparison(eps,MyEigenSort,&origin));

  /*
     Set solver parameters at runtime
  */
  PetscCall(EPSSetFromOptions(eps));
  PetscCall(PetscObjectTypeCompare((PetscObject)eps,EPSARNOLDI,&flg));
  if (flg) {
    PetscCall(EPSArnoldiGetDelayed(eps,&delay));
    if (delay) PetscCall(PetscPrintf(PETSC_COMM_WORLD," Warning: delayed reorthogonalization may be unstable\n"));
  }

  /*
     Set the initial vector. This is optional, if not done the initial
     vector is set to random values
  */
  PetscCall(MatCreateVecs(A,&v0,NULL));
  PetscCall(VecSetValue(v0,0,-1.5,INSERT_VALUES));
  PetscCall(VecSetValue(v0,1,2.1,INSERT_VALUES));
  PetscCall(VecAssemblyBegin(v0));
  PetscCall(VecAssemblyEnd(v0));
  PetscCall(EPSSetInitialSpace(eps,1,&v0));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSSolve(eps));
  PetscCall(EPSGetDimensions(eps,&nev,NULL,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %" PetscInt_FMT "\n",nev));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSErrorView(eps,EPS_ERROR_RELATIVE,NULL));
  if (!skipnorm) PetscCall(CheckNormalizedVectors(eps));
  PetscCall(EPSDestroy(&eps));
  PetscCall(MatDestroy(&A));
  PetscCall(VecDestroy(&v0));
  PetscCall(SlepcFinalize());
  return 0;
}

PetscErrorCode MatMarkovModel(PetscInt m,Mat A)
{
  const PetscReal cst = 0.5/(PetscReal)(m-1);
  PetscReal       pd,pu;
  PetscInt        Istart,Iend,i,j,jmax,ix=0;

  PetscFunctionBeginUser;
  PetscCall(MatGetOwnershipRange(A,&Istart,&Iend));
  for (i=1;i<=m;i++) {
    jmax = m-i+1;
    for (j=1;j<=jmax;j++) {
      ix = ix + 1;
      if (ix-1<Istart || ix>Iend) continue;  /* compute only owned rows */
      if (j!=jmax) {
        pd = cst*(PetscReal)(i+j-1);
        /* north */
        if (i==1) PetscCall(MatSetValue(A,ix-1,ix,2*pd,INSERT_VALUES));
        else PetscCall(MatSetValue(A,ix-1,ix,pd,INSERT_VALUES));
        /* east */
        if (j==1) PetscCall(MatSetValue(A,ix-1,ix+jmax-1,2*pd,INSERT_VALUES));
        else PetscCall(MatSetValue(A,ix-1,ix+jmax-1,pd,INSERT_VALUES));
      }
      /* south */
      pu = 0.5 - cst*(PetscReal)(i+j-3);
      if (j>1) PetscCall(MatSetValue(A,ix-1,ix-2,pu,INSERT_VALUES));
      /* west */
      if (i>1) PetscCall(MatSetValue(A,ix-1,ix-jmax-2,pu,INSERT_VALUES));
    }
  }
  PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
    Function for user-defined eigenvalue ordering criterion.

    Given two eigenvalues ar+i*ai and br+i*bi, the subroutine must choose
    one of them as the preferred one according to the criterion.
    In this example, the preferred value is the one furthest away from the origin.
*/
PetscErrorCode MyEigenSort(PetscScalar ar,PetscScalar ai,PetscScalar br,PetscScalar bi,PetscInt *r,void *ctx)
{
  PetscScalar origin = *(PetscScalar*)ctx;
  PetscReal   d;

  PetscFunctionBeginUser;
  d = (SlepcAbsEigenvalue(br-origin,bi) - SlepcAbsEigenvalue(ar-origin,ai))/PetscMax(SlepcAbsEigenvalue(ar-origin,ai),SlepcAbsEigenvalue(br-origin,bi));
  *r = d > PETSC_SQRT_MACHINE_EPSILON ? 1 : (d < -PETSC_SQRT_MACHINE_EPSILON ? -1 : PetscSign(PetscRealPart(br)));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*TEST

   testset:
      args: -eps_nev 4
      output_file: output/test9_1.out
      filter: sed -e "s/97136/97137/g"
      test:
         suffix: 1
         args: -eps_type {{krylovschur arnoldi lapack}} -eps_ncv 8 -eps_max_it 300
      test:
         suffix: 1_gd
         args: -eps_type gd -st_pc_type none
      test:
         suffix: 1_gd2
         args: -eps_type gd -eps_gd_double_expansion -st_pc_type none

   test:
      suffix: 2
      args: -eps_balance {{none oneside twoside}} -eps_krylovschur_locking {{0 1}} -eps_nev 4 -eps_max_it 1500
      requires: double
      output_file: output/test9_1.out

   test:
      suffix: 3
      nsize: 2
      args: -eps_type arnoldi -eps_arnoldi_delayed -eps_largest_real -eps_nev 3 -eps_tol 1e-7 -bv_orthog_refine {{never ifneeded}} -skipnorm
      requires: !single
      output_file: output/test9_3.out

   test:
      suffix: 4
      args: -eps_nev 4 -eps_true_residual
      requires: !single
      output_file: output/test9_1.out

   test:
      suffix: 5
      args: -eps_type jd -eps_nev 3 -eps_target .5 -eps_harmonic -st_ksp_type bicg -st_pc_type lu -eps_jd_minv 2
      filter: sed -e "s/[+-]0\.0*i//g"
      requires: !single

   test:
      suffix: 5_arpack
      args: -eps_nev 3 -st_type sinvert -eps_target .5 -eps_type arpack -eps_ncv 10
      requires: arpack !single
      output_file: output/test9_5.out

   testset:
      args: -eps_type ciss -eps_tol 1e-9 -rg_type ellipse -rg_ellipse_center 0.55 -rg_ellipse_radius 0.05 -rg_ellipse_vscale 0.1 -eps_ciss_usest 0 -eps_all
      requires: !single
      output_file: output/test9_6.out
      test:
         suffix: 6
      test:
         suffix: 6_hankel
         args: -eps_ciss_extraction hankel -eps_ciss_spurious_threshold 1e-6 -eps_ncv 64
      test:
         suffix: 6_cheby
         args: -eps_ciss_quadrule chebyshev
      test:
         suffix: 6_hankel_cheby
         args: -eps_ciss_extraction hankel -eps_ciss_quadrule chebyshev -eps_ncv 64
      test:
         suffix: 6_refine
         args: -eps_ciss_moments 4 -eps_ciss_blocksize 5 -eps_ciss_refine_inner 1 -eps_ciss_refine_blocksize 2
      test:
         suffix: 6_bcgs
         args: -eps_ciss_realmats -eps_ciss_ksp_type bcgs -eps_ciss_pc_type ilu -eps_ciss_integration_points 8

   test:
      suffix: 6_cheby_interval
      args: -eps_type ciss -eps_tol 1e-9 -rg_type interval -rg_interval_endpoints 0.5,0.6 -eps_ciss_quadrule chebyshev -eps_ciss_usest 0 -eps_all
      requires: !single
      output_file: output/test9_6.out

   testset:
      args: -eps_nev 4 -eps_two_sided -eps_view_vectors ::ascii_info -eps_view_values
      filter: sed -e "s/\(0x[0-9a-fA-F]*\)/objectid/"
      test:
         suffix: 7_real
         requires: !single !complex
      test:
         suffix: 7
         requires: !single complex

   test:
      suffix: 8
      args: -eps_nev 4 -eps_ncv 7 -eps_view_values draw -eps_monitor draw::draw_lg
      requires: x
      output_file: output/test9_1.out

   test:
      suffix: 5_ksphpddm
      args: -eps_nev 3 -st_type sinvert -eps_target .5 -st_ksp_type hpddm -st_ksp_hpddm_type gcrodr -eps_ncv 10
      requires: hpddm
      output_file: output/test9_5.out

   test:
      suffix: 5_pchpddm
      args: -eps_nev 3 -st_type sinvert -eps_target .5 -st_pc_type hpddm -st_pc_hpddm_coarse_pc_type lu -eps_ncv 10
      requires: hpddm
      output_file: output/test9_5.out

TEST*/
