/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Solves a linear system using PCHPDDM.\n\n"
  "Modification of ${PETSC_DIR}/src/ksp/ksp/tutorials/ex76.c where concurrent EPS are instantiated explicitly by the user.\n\n";

#include <slepceps.h>

int main(int argc,char **args)
{
  Mat            A,aux,a,P,B,X;
  Vec            b;
  KSP            ksp;
  PC             pc;
  EPS            eps;
  ST             st;
  IS             is,sizes;
  const PetscInt *idx;
  PetscInt       m,rstart,rend,location,nev,nconv;
  PetscMPIInt    rank,size;
  PetscViewer    viewer;
  char           dir[PETSC_MAX_PATH_LEN],name[PETSC_MAX_PATH_LEN];
  PetscBool      flg;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&args,NULL,help));
  PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD,&size));
  PetscCheck(size == 4,PETSC_COMM_WORLD,PETSC_ERR_USER,"This example requires 4 processes");
  PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD,&rank));
  PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
  PetscCall(MatCreate(PETSC_COMM_SELF,&aux));
  PetscCall(ISCreate(PETSC_COMM_SELF,&is));
  PetscCall(PetscStrcpy(dir,"."));
  PetscCall(PetscOptionsGetString(NULL,NULL,"-load_dir",dir,sizeof(dir),NULL));
  /* loading matrices */
  PetscCall(PetscSNPrintf(name,sizeof(name),"%s/sizes_%d_%d.dat",dir,rank,size));
  PetscCall(PetscViewerBinaryOpen(PETSC_COMM_SELF,name,FILE_MODE_READ,&viewer));
  PetscCall(ISCreate(PETSC_COMM_SELF,&sizes));
  PetscCall(ISLoad(sizes,viewer));
  PetscCall(ISGetIndices(sizes,&idx));
  PetscCall(MatSetSizes(A,idx[0],idx[1],idx[2],idx[3]));
  PetscCall(ISRestoreIndices(sizes,&idx));
  PetscCall(ISDestroy(&sizes));
  PetscCall(PetscViewerDestroy(&viewer));
  PetscCall(PetscSNPrintf(name,sizeof(name),"%s/A.dat",dir));
  PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD,name,FILE_MODE_READ,&viewer));
  PetscCall(MatLoad(A,viewer));
  PetscCall(PetscViewerDestroy(&viewer));
  PetscCall(PetscSNPrintf(name,sizeof(name),"%s/is_%d_%d.dat",dir,rank,size));
  PetscCall(PetscViewerBinaryOpen(PETSC_COMM_SELF,name,FILE_MODE_READ,&viewer));
  PetscCall(ISLoad(is,viewer));
  PetscCall(PetscViewerDestroy(&viewer));
  PetscCall(ISSetBlockSize(is,2));
  PetscCall(PetscSNPrintf(name,sizeof(name),"%s/Neumann_%d_%d.dat",dir,rank,size));
  PetscCall(PetscViewerBinaryOpen(PETSC_COMM_SELF,name,FILE_MODE_READ,&viewer));
  PetscCall(MatLoad(aux,viewer));
  PetscCall(PetscViewerDestroy(&viewer));
  PetscCall(MatSetBlockSizesFromMats(aux,A,A));
  PetscCall(MatSetOption(A,MAT_SYMMETRIC,PETSC_TRUE));
  PetscCall(MatSetOption(aux,MAT_SYMMETRIC,PETSC_TRUE));
  /* ready for testing */
  PetscCall(KSPCreate(PETSC_COMM_WORLD,&ksp));
  PetscCall(KSPSetOperators(ksp,A,A));
  PetscCall(KSPGetPC(ksp,&pc));
  PetscCall(PCSetType(pc,PCHPDDM));
  PetscCall(MatDuplicate(aux,MAT_DO_NOT_COPY_VALUES,&B)); /* duplicate so that MatStructure is SAME_NONZERO_PATTERN */
  PetscCall(MatGetDiagonalBlock(A,&a));
  PetscCall(MatGetOwnershipRange(A,&rstart,&rend));
  PetscCall(ISGetLocalSize(is,&m));
  PetscCall(MatCreateSeqAIJ(PETSC_COMM_SELF,rend-rstart,m,1,NULL,&P));
  for (m = rstart; m < rend; ++m) {
    PetscCall(ISLocate(is,m,&location));
    PetscCheck(location >= 0,PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"IS of the auxiliary Mat does not include all local rows of A");
    PetscCall(MatSetValue(P,m-rstart,location,1.0,INSERT_VALUES));
  }
  PetscCall(MatAssemblyBegin(P,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(P,MAT_FINAL_ASSEMBLY));
  PetscCall(MatPtAP(a,P,MAT_INITIAL_MATRIX,1.0,&X));
  PetscCall(MatDestroy(&P));
  PetscCall(MatAXPY(B,1.0,X,SUBSET_NONZERO_PATTERN));
  PetscCall(MatDestroy(&X));
  PetscCall(MatSetOption(B,MAT_SYMMETRIC,PETSC_TRUE));
  PetscCall(EPSCreate(PETSC_COMM_SELF,&eps));
  PetscCall(EPSSetOperators(eps,aux,B));
  PetscCall(EPSSetProblemType(eps,EPS_GHEP));
  PetscCall(EPSSetTarget(eps,0.0));
  PetscCall(EPSSetWhichEigenpairs(eps,EPS_TARGET_MAGNITUDE));
  PetscCall(EPSGetST(eps,&st));
  PetscCall(STSetType(st,STSINVERT));
  PetscCall(EPSSetFromOptions(eps));
  PetscCall(EPSGetDimensions(eps,&nev,NULL,NULL));
  PetscCall(EPSSolve(eps));
  PetscCall(EPSGetConverged(eps,&nconv));
  nev = PetscMin(nev,nconv);
  PetscCall(ISGetLocalSize(is,&m));
  PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,m,nev,NULL,&P));
  for (m = 0; m < nev; ++m) {
    PetscCall(MatDenseGetColumnVecWrite(P,m,&b));
    PetscCall(EPSGetEigenvector(eps,m,b,NULL));
    PetscCall(MatDenseRestoreColumnVecWrite(P,m,&b));
  }
  PetscCall(EPSDestroy(&eps));
  PetscCall(MatDestroy(&B));
  PetscCall(PCHPDDMSetDeflationMat(pc,is,P));
  PetscCall(MatDestroy(&P));
  PetscCall(MatDestroy(&aux));
  PetscCall(KSPSetFromOptions(ksp));
  PetscCall(ISDestroy(&is));
  PetscCall(MatCreateVecs(A,NULL,&b));
  PetscCall(VecSet(b,1.0));
  PetscCall(KSPSolve(ksp,b,b));
  PetscCall(PetscObjectTypeCompare((PetscObject)pc,PCHPDDM,&flg));
  if (flg) {
    PetscCall(PCHPDDMGetSTShareSubKSP(pc,&flg));
    /* since EPSSolve() is called outside PCSetUp_HPDDM() and there is not mechanism (yet) to pass the underlying ST, */
    /* PCHPDDMGetSTShareSubKSP() should not return PETSC_TRUE when using PCHPDDMSetDeflationMat()                     */
    PetscCheck(!flg,PETSC_COMM_SELF,PETSC_ERR_PLIB,"PCHPDDMGetSTShareSubKSP() should not return PETSC_TRUE");
  }
  PetscCall(VecDestroy(&b));
  PetscCall(KSPDestroy(&ksp));
  PetscCall(MatDestroy(&A));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   build:
      requires: hpddm

   test:
      requires: double !complex !defined(PETSC_USE_64BIT_INDICES) defined(PETSC_HAVE_DYNAMIC_LIBRARIES)
      nsize: 4
      filter: grep -v "    total: nonzeros=" | grep -v "    rows=" | grep -v "       factor fill ratio given " | grep -v "      using I-node" | sed -e "s/Linear solve converged due to CONVERGED_RTOL iterations 5/Linear solve converged due to CONVERGED_RTOL iterations 4/g" -e "s/amount of overlap = 1/user-defined overlap/g"
      # similar output as ex76 with -pc_hpddm_levels_eps_nev val instead of just -eps_nev val
      args: -ksp_rtol 1e-3 -ksp_max_it 10 -ksp_error_if_not_converged -ksp_converged_reason -ksp_view -ksp_type hpddm -pc_type hpddm -pc_hpddm_levels_1_sub_pc_type lu -pc_hpddm_define_subdomains {{false true}shared output} -pc_hpddm_levels_1_pc_type asm -pc_hpddm_levels_1_sub_pc_type lu -pc_hpddm_levels_1_st_share_sub_ksp -eps_nev 10 -load_dir ${DATAFILESPATH}/matrices/hpddm/GENEO -matload_block_size 2

TEST*/
