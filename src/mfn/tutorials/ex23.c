/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Computes exp(t*A)*v for a matrix associated with a Markov model.\n\n"
  "The command line options are:\n"
  "  -t <t>, where <t> = time parameter (multiplies the matrix).\n"
  "  -m <m>, where <m> = number of grid subdivisions in each dimension.\n\n"
  "  -shift <shift>, where <shift> = diagonal shift applied to the matrix.\n\n"
  "To draw the solution run with -mfn_view_solution draw -draw_pause -1\n\n";

#include <slepcmfn.h>

/*
   User-defined routines
*/
PetscErrorCode MatMarkovModel(PetscInt m,Mat A);

int main(int argc,char **argv)
{
  Mat                A;           /* problem matrix */
  MFN                mfn;
  FN                 f;
  PetscReal          tol,norm;
  PetscScalar        t=2.0,shift=0.0;
  Vec                v,y;
  PetscInt           N,m=15,ncv,maxit,its;
  MFNConvergedReason reason;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));

  PetscCall(PetscOptionsGetInt(NULL,NULL,"-m",&m,NULL));
  PetscCall(PetscOptionsGetScalar(NULL,NULL,"-t",&t,NULL));
  PetscCall(PetscOptionsGetScalar(NULL,NULL,"-shift",&shift,NULL));
  N = m*(m+1)/2;
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nMarkov y=exp(t*A)*e_1, N=%" PetscInt_FMT " (m=%" PetscInt_FMT ")",N,m));
  if (shift!=0.0) {
    if (PetscImaginaryPart(shift)==0.0) PetscCall(PetscPrintf(PETSC_COMM_WORLD,", shift=%g",(double)PetscRealPart(shift)));
    else PetscCall(PetscPrintf(PETSC_COMM_WORLD,", shift=%g%+gi",(double)PetscRealPart(shift),(double)PetscImaginaryPart(shift)));
  }
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n\n"));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            Compute the transition probability matrix, A
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
  PetscCall(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,N,N));
  PetscCall(MatSetFromOptions(A));
  PetscCall(MatMarkovModel(m,A));
  PetscCall(MatShift(A,shift));

  /* set v = e_1 */
  PetscCall(MatCreateVecs(A,NULL,&y));
  PetscCall(MatCreateVecs(A,NULL,&v));
  PetscCall(VecSetValue(v,0,1.0,INSERT_VALUES));
  PetscCall(VecAssemblyBegin(v));
  PetscCall(VecAssemblyEnd(v));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the solver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Create matrix function solver context
  */
  PetscCall(MFNCreate(PETSC_COMM_WORLD,&mfn));

  /*
     Set operator matrix, the function to compute, and other options
  */
  PetscCall(MFNSetOperator(mfn,A));
  PetscCall(MFNGetFN(mfn,&f));
  PetscCall(FNSetType(f,FNEXP));
  PetscCall(FNSetScale(f,t,1.0));
  PetscCall(MFNSetTolerances(mfn,1e-07,PETSC_CURRENT));

  /*
     Set solver parameters at runtime
  */
  PetscCall(MFNSetFromOptions(mfn));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the problem, y=exp(t*A)*v
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(MFNSolve(mfn,v,y));
  PetscCall(MFNGetConvergedReason(mfn,&reason));
  PetscCheck(reason>=0,PETSC_COMM_WORLD,PETSC_ERR_CONV_FAILED,"Solver did not converge");
  PetscCall(VecNorm(y,NORM_2,&norm));

  /*
     Optional: Get some information from the solver and display it
  */
  PetscCall(MFNGetIterationNumber(mfn,&its));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %" PetscInt_FMT "\n",its));
  PetscCall(MFNGetDimensions(mfn,&ncv));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Subspace dimension: %" PetscInt_FMT "\n",ncv));
  PetscCall(MFNGetTolerances(mfn,&tol,&maxit));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g, maxit=%" PetscInt_FMT "\n",(double)tol,maxit));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Computed vector at time t=%.4g has norm %g\n\n",(double)PetscRealPart(t),(double)norm));

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

/*
    Matrix generator for a Markov model of a random walk on a triangular grid.
    See ex5.c for additional details.
*/
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

/*TEST

   test:
      suffix: 1
      args: -mfn_ncv 6

   test:
      suffix: 2
      requires: complex
      args: -mfn_type cram -shift -2
      filter: sed -e "s/0.0452997/0.0452998/"

TEST*/
