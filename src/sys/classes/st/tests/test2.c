/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test ST with one matrix.\n\n";

#include <slepcst.h>

int main(int argc,char **argv)
{
  Mat            A,B,mat[1];
  ST             st;
  Vec            v,w;
  STType         type;
  PetscScalar    sigma,tau,val;
  PetscInt       n=10,i,Istart,Iend;
  PetscBool      test_compl=PETSC_FALSE;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscOptionsGetBool(NULL,NULL,"-complex",&test_compl,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n1-D Laplacian, n=%" PetscInt_FMT "\n\n",n));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the operator matrix for the 1-D Laplacian
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
  PetscCall(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetFromOptions(A));

  PetscCall(MatGetOwnershipRange(A,&Istart,&Iend));
#if PetscDefined(USE_COMPLEX)
  val = test_compl? PetscCMPLX(-1.0,0.4): -1.0;
#else
  val = -1.0;
#endif
  for (i=Istart;i<Iend;i++) {
    if (i>0) PetscCall(MatSetValue(A,i,i-1,val,INSERT_VALUES));
    if (i<n-1) PetscCall(MatSetValue(A,i,i+1,PetscConj(val),INSERT_VALUES));
    PetscCall(MatSetValue(A,i,i,2.0,INSERT_VALUES));
  }
  PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));
  PetscCall(MatCreateVecs(A,&v,&w));
#if PetscDefined(USE_COMPLEX)
  val = test_compl? PetscCMPLX(1.0,-0.01): 1.0;
#else
  val = 1.0;
#endif
  PetscCall(VecSet(v,val));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the spectral transformation object
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(STCreate(PETSC_COMM_WORLD,&st));
  mat[0] = A;
  PetscCall(STSetMatrices(st,1,mat));
  PetscCall(STSetTransform(st,PETSC_TRUE));
  PetscCall(STSetFromOptions(st));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
              Apply the transformed operator for several ST's
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* shift, sigma=0.0 */
  PetscCall(STSetUp(st));
  PetscCall(STGetType(st,&type));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"ST type %s\n",type));
  PetscCall(STApply(st,v,w));
  PetscCall(VecView(w,NULL));

  /* shift, sigma=0.1 */
  sigma = 0.1;
  PetscCall(STSetShift(st,sigma));
  PetscCall(STGetShift(st,&sigma));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"With shift=%g\n",(double)PetscRealPart(sigma)));
  PetscCall(STApply(st,v,w));
  PetscCall(VecView(w,NULL));
  PetscCall(STApplyTranspose(st,v,w));
  PetscCall(VecView(w,NULL));
  if (test_compl) {
    PetscCall(STApplyHermitianTranspose(st,v,w));
    PetscCall(VecView(w,NULL));
  }

  /* sinvert, sigma=0.1 */
  PetscCall(STPostSolve(st));   /* undo changes if inplace */
  PetscCall(STSetType(st,STSINVERT));
  PetscCall(STGetType(st,&type));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"ST type %s\n",type));
  PetscCall(STGetShift(st,&sigma));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"With shift=%g\n",(double)PetscRealPart(sigma)));
  PetscCall(STApply(st,v,w));
  PetscCall(VecView(w,NULL));

  /* sinvert, sigma=-0.5 */
  sigma = -0.5;
  PetscCall(STSetShift(st,sigma));
  PetscCall(STGetShift(st,&sigma));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"With shift=%g\n",(double)PetscRealPart(sigma)));
  PetscCall(STApply(st,v,w));
  PetscCall(VecView(w,NULL));
  if (test_compl) {
    PetscCall(STApplyHermitianTranspose(st,v,w));
    PetscCall(VecView(w,NULL));
  }

  /* cayley, sigma=-0.5, tau=-0.5 (equal to sigma by default) */
  PetscCall(STPostSolve(st));   /* undo changes if inplace */
  PetscCall(STSetType(st,STCAYLEY));
  PetscCall(STSetUp(st));
  PetscCall(STGetType(st,&type));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"ST type %s\n",type));
  PetscCall(STGetShift(st,&sigma));
  PetscCall(STCayleyGetAntishift(st,&tau));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"With shift=%g, antishift=%g\n",(double)PetscRealPart(sigma),(double)PetscRealPart(tau)));
  PetscCall(STApply(st,v,w));
  PetscCall(VecView(w,NULL));

  /* cayley, sigma=1.1, tau=1.1 (still equal to sigma) */
  sigma = 1.1;
  PetscCall(STSetShift(st,sigma));
  PetscCall(STGetShift(st,&sigma));
  PetscCall(STCayleyGetAntishift(st,&tau));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"With shift=%g, antishift=%g\n",(double)PetscRealPart(sigma),(double)PetscRealPart(tau)));
  PetscCall(STApply(st,v,w));
  PetscCall(VecView(w,NULL));

  /* cayley, sigma=1.1, tau=-1.0 */
  tau = -1.0;
  PetscCall(STCayleySetAntishift(st,tau));
  PetscCall(STGetShift(st,&sigma));
  PetscCall(STCayleyGetAntishift(st,&tau));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"With shift=%g, antishift=%g\n",(double)PetscRealPart(sigma),(double)PetscRealPart(tau)));
  PetscCall(STApply(st,v,w));
  PetscCall(VecView(w,NULL));

  /* check inner product matrix in Cayley */
  PetscCall(STGetBilinearForm(st,&B));
  PetscCall(MatMult(B,v,w));
  PetscCall(VecView(w,NULL));

  /* check again sinvert, sigma=0.1 */
  PetscCall(STPostSolve(st));   /* undo changes if inplace */
  PetscCall(STSetType(st,STSINVERT));
  PetscCall(STGetType(st,&type));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"ST type %s\n",type));
  sigma = 0.1;
  PetscCall(STSetShift(st,sigma));
  PetscCall(STGetShift(st,&sigma));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"With shift=%g\n",(double)PetscRealPart(sigma)));
  PetscCall(STApply(st,v,w));
  PetscCall(VecView(w,NULL));

  PetscCall(STDestroy(&st));
  PetscCall(MatDestroy(&A));
  PetscCall(MatDestroy(&B));
  PetscCall(VecDestroy(&v));
  PetscCall(VecDestroy(&w));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   test:
      suffix: 1
      args: -st_matmode {{copy inplace shell}}
      requires: !single

   test:
      suffix: 2
      args: -complex -st_matmode {{copy inplace shell}}
      requires: complex !single
      filter: sed -e "s/3.46945e-18 - 8.67362e-18 i/0./"

TEST*/
