/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Computes a partial generalized singular value decomposition (GSVD).\n\n"
  "The command line options are:\n"
  "  -m <m>, where <m> = number of rows of A.\n"
  "  -n <n>, where <n> = number of columns of A.\n"
  "  -p <p>, where <p> = number of rows of B.\n\n";

#include <slepcsvd.h>

int main(int argc,char **argv)
{
  Mat            A,B;             /* operator matrices */
  Vec            u,v,x;           /* singular vectors */
  SVD            svd;             /* singular value problem solver context */
  SVDType        type;
  Vec            uv,aux[2],*U,*V;
  PetscReal      error,tol,sigma,lev1=0.0,lev2=0.0;
  PetscInt       m=100,n,p=14,i,j,d,Istart,Iend,nsv,maxit,its,nconv;
  PetscBool      flg,skiporth=PETSC_FALSE;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));

  PetscCall(PetscOptionsGetInt(NULL,NULL,"-m",&m,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,&flg));
  if (!flg) n = m;
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-p",&p,&flg));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nGeneralized singular value decomposition, (%" PetscInt_FMT "+%" PetscInt_FMT ")x%" PetscInt_FMT "\n\n",m,p,n));
  PetscCall(PetscOptionsGetBool(NULL,NULL,"-skiporth",&skiporth,NULL));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                          Build the matrices
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
  PetscCall(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,m,n));
  PetscCall(MatSetFromOptions(A));

  PetscCall(MatGetOwnershipRange(A,&Istart,&Iend));
  for (i=Istart;i<Iend;i++) {
    if (i>0 && i-1<n) PetscCall(MatSetValue(A,i,i-1,-1.0,INSERT_VALUES));
    if (i+1<n) PetscCall(MatSetValue(A,i,i+1,-1.0,INSERT_VALUES));
    if (i<n) PetscCall(MatSetValue(A,i,i,2.0,INSERT_VALUES));
    if (i>n) PetscCall(MatSetValue(A,i,n-1,1.0,INSERT_VALUES));
  }
  PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));

  PetscCall(MatCreate(PETSC_COMM_WORLD,&B));
  PetscCall(MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,p,n));
  PetscCall(MatSetFromOptions(B));

  PetscCall(MatGetOwnershipRange(B,&Istart,&Iend));
  d = PetscMax(0,n-p);
  for (i=Istart;i<Iend;i++) {
    for (j=0;j<=PetscMin(i,n-1);j++) PetscCall(MatSetValue(B,i,j+d,1.0,INSERT_VALUES));
  }
  PetscCall(MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          Create the singular value solver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Create singular value solver context
  */
  PetscCall(SVDCreate(PETSC_COMM_WORLD,&svd));

  /*
     Set operators and problem type
  */
  PetscCall(SVDSetOperators(svd,A,B));
  PetscCall(SVDSetProblemType(svd,SVD_GENERALIZED));

  /*
     Set solver parameters at runtime
  */
  PetscCall(SVDSetFromOptions(svd));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the singular value system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(SVDSolve(svd));
  PetscCall(SVDGetIterationNumber(svd,&its));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %" PetscInt_FMT "\n",its));

  /*
     Optional: Get some information from the solver and display it
  */
  PetscCall(SVDGetType(svd,&type));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type));
  PetscCall(SVDGetDimensions(svd,&nsv,NULL,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Number of requested generalized singular values: %" PetscInt_FMT "\n",nsv));
  PetscCall(SVDGetTolerances(svd,&tol,&maxit));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g, maxit=%" PetscInt_FMT "\n",(double)tol,maxit));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Get number of converged singular triplets
  */
  PetscCall(SVDGetConverged(svd,&nconv));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Number of converged approximate singular triplets: %" PetscInt_FMT "\n\n",nconv));

  if (nconv>0) {
    /*
       Create vectors. The interface returns u and v as stacked on top of each other
       [u;v] so need to create a special vector (VecNest) to extract them
    */
    PetscCall(MatCreateVecs(A,&x,&u));
    PetscCall(MatCreateVecs(B,NULL,&v));
    aux[0] = u;
    aux[1] = v;
    PetscCall(VecCreateNest(PETSC_COMM_WORLD,2,NULL,aux,&uv));

    PetscCall(VecDuplicateVecs(u,nconv,&U));
    PetscCall(VecDuplicateVecs(v,nconv,&V));

    /*
       Display singular values and errors relative to the norms
    */
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
         "          sigma           ||r||/||[A;B]||\n"
         "  --------------------- ------------------\n"));
    for (i=0;i<nconv;i++) {
      /*
         Get converged singular triplets: i-th singular value is stored in sigma
      */
      PetscCall(SVDGetSingularTriplet(svd,i,&sigma,uv,x));

      /* at this point, u and v can be used normally as individual vectors */
      PetscCall(VecCopy(u,U[i]));
      PetscCall(VecCopy(v,V[i]));

      /*
         Compute the error associated to each singular triplet
      */
      PetscCall(SVDComputeError(svd,i,SVD_ERROR_NORM,&error));

      PetscCall(PetscPrintf(PETSC_COMM_WORLD,"       % 6f      ",(double)sigma));
      PetscCall(PetscPrintf(PETSC_COMM_WORLD,"   % 12g\n",(double)error));
    }
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n"));

    if (!skiporth) {
      PetscCall(VecCheckOrthonormality(U,nconv,NULL,nconv,NULL,NULL,&lev1));
      PetscCall(VecCheckOrthonormality(V,nconv,NULL,nconv,NULL,NULL,&lev2));
    }
    if (lev1+lev2<20*tol) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Level of orthogonality below the tolerance\n"));
    else PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Level of orthogonality: %g (U) %g (V)\n",(double)lev1,(double)lev2));
    PetscCall(VecDestroyVecs(nconv,&U));
    PetscCall(VecDestroyVecs(nconv,&V));
    PetscCall(VecDestroy(&x));
    PetscCall(VecDestroy(&u));
    PetscCall(VecDestroy(&v));
    PetscCall(VecDestroy(&uv));
  }

  /*
     Free work space
  */
  PetscCall(SVDDestroy(&svd));
  PetscCall(MatDestroy(&A));
  PetscCall(MatDestroy(&B));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   testset:
      filter: grep -v "Solution method" | grep -v "Number of iterations" | sed -e "s/, maxit=1[0]*$//" | sed -e "s/[0-9]\.[0-9]*e[+-]\([0-9]*\)/removed/g"
      requires: double
      test:
         args: -svd_type lapack -m 20 -n 10 -p 6
         suffix: 1
      test:
         args: -svd_type lapack -m 15 -n 20 -p 10 -svd_smallest
         suffix: 2
      test:
         args: -svd_type lapack -m 15 -n 20 -p 21
         suffix: 3
      test:
         args: -svd_type lapack -m 20 -n 15 -p 21
         suffix: 4

   testset:
      args: -m 25 -n 20 -p 21 -svd_smallest -svd_nsv 2
      filter: grep -v "Solution method" | grep -v "Number of iterations" | sed -e "s/, maxit=1[0]*$//" | sed -e "s/[0-9]\.[0-9]*e[+-]\([0-9]*\)/removed/g"
      requires: double
      output_file: output/ex45_5.out
      test:
         args: -svd_type trlanczos -svd_ncv 8 -svd_trlanczos_gbidiag {{upper lower}} -svd_trlanczos_oneside {{0 1}}
         suffix: 5
      test:
         args: -svd_type cross -svd_ncv 10 -svd_cross_explicitmatrix
         suffix: 5_cross
      test:
         args: -svd_type cross -svd_ncv 10 -svd_cross_eps_type krylovschur -svd_cross_st_type sinvert -svd_cross_eps_target 0 -svd_cross_st_ksp_type gmres -svd_cross_st_pc_type jacobi
         suffix: 5_cross_implicit
      test:
         args: -svd_type cyclic -svd_ncv 12 -svd_cyclic_explicitmatrix {{0 1}}
         suffix: 5_cyclic
         requires: !complex

   testset:
      args: -m 15 -n 20 -p 21 -svd_nsv 4 -svd_ncv 9
      filter: grep -v "Solution method" | grep -v "Number of iterations" | sed -e "s/7.884967/7.884968/" | sed -e "s/, maxit=1[0]*$//" | sed -e "s/[0-9]\.[0-9]*e[+-]\([0-9]*\)/removed/g"
      requires: double
      output_file: output/ex45_6.out
      test:
         args: -svd_type trlanczos -svd_trlanczos_gbidiag {{single upper lower}} -svd_trlanczos_locking {{0 1}} -svd_trlanczos_oneside {{0 1}}
         suffix: 6
      test:
         args: -svd_type cross -svd_cross_explicitmatrix {{0 1}}
         suffix: 6_cross

   test:
      args: -m 15 -n 20 -p 21 -svd_nsv 4 -svd_ncv 9 -svd_type cyclic -svd_cyclic_explicitmatrix {{0 1}}
      filter: grep -v "Number of iterations" | sed -e "s/7.884967/7.884968/" | sed -e "s/[0-9]\.[0-9]*e[+-]\([0-9]*\)/removed/g"
      requires: double
      suffix: 6_cyclic
      output_file: output/ex45_6_cyclic.out

   testset:
      args: -m 20 -n 15 -p 21 -svd_nsv 4 -svd_ncv 9
      filter: grep -v "Solution method" | grep -v "Number of iterations" | sed -e "s/, maxit=1[0]*$//" | sed -e "s/[0-9]\.[0-9]*e[+-]\([0-9]*\)/removed/g"
      requires: double
      output_file: output/ex45_7.out
      test:
         args: -svd_type trlanczos -svd_trlanczos_gbidiag {{single upper lower}} -svd_trlanczos_restart 0.4 -svd_trlanczos_oneside {{0 1}}
         suffix: 7
      test:
         args: -svd_type cross -svd_cross_explicitmatrix {{0 1}}
         suffix: 7_cross

   test:
      args: -m 20 -n 15 -p 21 -svd_nsv 4 -svd_ncv 16 -svd_type cyclic -svd_cyclic_explicitmatrix {{0 1}}
      filter: grep -v "Number of iterations" | sed -e "s/[0-9]\.[0-9]*e[+-]\([0-9]*\)/removed/g"
      requires: double
      suffix: 7_cyclic
      output_file: output/ex45_7_cyclic.out

   test:
       args: -m 25 -n 20 -p 21 -svd_smallest -svd_nsv 2 -svd_ncv 5 -svd_type trlanczos -svd_trlanczos_gbidiag {{upper lower}} -svd_trlanczos_scale {{0.1 -20}}
       filter: grep -v "Solution method" | grep -v "Number of iterations" | grep -v "Stopping condition" | sed -e "s/, maxit=1[0]*$//" | sed -e "s/[0-9]\.[0-9]*e[+-]\([0-9]*\)/removed/g"
       suffix: 8

TEST*/
