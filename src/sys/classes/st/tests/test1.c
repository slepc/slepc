/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test ST with shell matrices.\n\n";

#include <slepcst.h>

static PetscErrorCode MatMult_Shell(Mat S,Vec x,Vec y);
static PetscErrorCode MatMultTranspose_Shell(Mat S,Vec x,Vec y);
#if defined(PETSC_USE_COMPLEX)
static PetscErrorCode MatMultHermitianTranspose_Shell(Mat S,Vec x,Vec y);
#endif
static PetscErrorCode MatGetDiagonal_Shell(Mat S,Vec diag);
static PetscErrorCode MatDuplicate_Shell(Mat S,MatDuplicateOption op,Mat *M);

static PetscErrorCode MyShellMatCreate(Mat *A,Mat *M)
{
  MPI_Comm       comm;
  PetscInt       n;

  PetscFunctionBeginUser;
  PetscCall(MatGetSize(*A,&n,NULL));
  PetscCall(PetscObjectGetComm((PetscObject)*A,&comm));
  PetscCall(MatCreateShell(comm,PETSC_DECIDE,PETSC_DECIDE,n,n,A,M));
  PetscCall(MatShellSetOperation(*M,MATOP_MULT,(PetscErrorCodeFn*)MatMult_Shell));
  PetscCall(MatShellSetOperation(*M,MATOP_MULT_TRANSPOSE,(PetscErrorCodeFn*)MatMultTranspose_Shell));
#if defined(PETSC_USE_COMPLEX)
  PetscCall(MatShellSetOperation(*M,MATOP_MULT_HERMITIAN_TRANSPOSE,(PetscErrorCodeFn*)MatMultHermitianTranspose_Shell));
#endif
  PetscCall(MatShellSetOperation(*M,MATOP_GET_DIAGONAL,(PetscErrorCodeFn*)MatGetDiagonal_Shell));
  PetscCall(MatShellSetOperation(*M,MATOP_DUPLICATE,(PetscErrorCodeFn*)MatDuplicate_Shell));
  PetscFunctionReturn(PETSC_SUCCESS);
}

int main(int argc,char **argv)
{
  Mat            A,S,mat[1];
  ST             st;
  Vec            v,w;
  STType         type;
  KSP            ksp;
  PC             pc;
  PetscScalar    sigma;
  PetscInt       n=10,i,Istart,Iend;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n1-D Laplacian with shell matrices, n=%" PetscInt_FMT "\n\n",n));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the operator matrix for the 1-D Laplacian
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
  PetscCall(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetFromOptions(A));

  PetscCall(MatGetOwnershipRange(A,&Istart,&Iend));
  for (i=Istart;i<Iend;i++) {
    if (i>0) PetscCall(MatSetValue(A,i,i-1,-1.0,INSERT_VALUES));
    if (i<n-1) PetscCall(MatSetValue(A,i,i+1,-1.0,INSERT_VALUES));
    PetscCall(MatSetValue(A,i,i,2.0,INSERT_VALUES));
  }
  PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));

  /* create the shell version of A */
  PetscCall(MyShellMatCreate(&A,&S));

  /* work vectors */
  PetscCall(MatCreateVecs(A,&v,&w));
  PetscCall(VecSet(v,1.0));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the spectral transformation object
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(STCreate(PETSC_COMM_WORLD,&st));
  mat[0] = S;
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
  PetscCall(STApplyHermitianTranspose(st,v,w));
  PetscCall(VecView(w,NULL));

  /* shift, sigma=0.1 */
  sigma = 0.1;
  PetscCall(STSetShift(st,sigma));
  PetscCall(STGetShift(st,&sigma));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"With shift=%g\n",(double)PetscRealPart(sigma)));
  PetscCall(STApply(st,v,w));
  PetscCall(VecView(w,NULL));

  /* sinvert, sigma=0.1 */
  PetscCall(STPostSolve(st));   /* undo changes if inplace */
  PetscCall(STSetType(st,STSINVERT));
  PetscCall(STGetKSP(st,&ksp));
  PetscCall(KSPSetType(ksp,KSPGMRES));
  PetscCall(KSPGetPC(ksp,&pc));
  PetscCall(PCSetType(pc,PCJACOBI));
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

  PetscCall(STDestroy(&st));
  PetscCall(MatDestroy(&A));
  PetscCall(MatDestroy(&S));
  PetscCall(VecDestroy(&v));
  PetscCall(VecDestroy(&w));
  PetscCall(SlepcFinalize());
  return 0;
}

static PetscErrorCode MatMult_Shell(Mat S,Vec x,Vec y)
{
  Mat               *A;

  PetscFunctionBeginUser;
  PetscCall(MatShellGetContext(S,&A));
  PetscCall(MatMult(*A,x,y));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatMultTranspose_Shell(Mat S,Vec x,Vec y)
{
  Mat               *A;

  PetscFunctionBeginUser;
  PetscCall(MatShellGetContext(S,&A));
  PetscCall(MatMultTranspose(*A,x,y));
  PetscFunctionReturn(PETSC_SUCCESS);
}

#if defined(PETSC_USE_COMPLEX)
static PetscErrorCode MatMultHermitianTranspose_Shell(Mat S,Vec x,Vec y)
{
  Mat               *A;

  PetscFunctionBeginUser;
  PetscCall(MatShellGetContext(S,&A));
  PetscCall(MatMultHermitianTranspose(*A,x,y));
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif

static PetscErrorCode MatGetDiagonal_Shell(Mat S,Vec diag)
{
  Mat               *A;

  PetscFunctionBeginUser;
  PetscCall(MatShellGetContext(S,&A));
  PetscCall(MatGetDiagonal(*A,diag));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatDuplicate_Shell(Mat S,MatDuplicateOption op,Mat *M)
{
  Mat            *A;

  PetscFunctionBeginUser;
  PetscCall(MatShellGetContext(S,&A));
  PetscCall(MyShellMatCreate(A,M));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*TEST

   test:
      suffix: 1
      args: -st_matmode {{inplace shell}}
      requires: !single

TEST*/
