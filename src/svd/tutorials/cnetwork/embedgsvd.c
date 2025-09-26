/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Computes a GSVD associated with a complex network analysis problem.\n\n"
  "The command line options are:\n"
  "  -file <filename>, text file in TSV format containing undirected graph.\n"
  "  -beta <beta>, optional decay parameter.\n\n";

/*
    Computes the Katz high-order proximity embedding of a graph via a partial GSVD.

    [1] M. Ou et al. "Asymmetric transitivity preserving graph embedding" Proc. of
        ACM SIGKDD 2016 - https://doi.org/10.1145/2939672.2939751
*/

#include <slepcsvd.h>
#include "network.h"

PetscErrorCode SpectralRadius(Mat A,PetscReal *rho)
{
  EPS         eps;
  PetscInt    nconv;
  PetscScalar kr,ki;

  PetscFunctionBeginUser;
  PetscCall(EPSCreate(PETSC_COMM_WORLD,&eps));
  PetscCall(EPSSetOperators(eps,A,NULL));
  PetscCall(EPSSetProblemType(eps,EPS_NHEP));
  PetscCall(EPSSetDimensions(eps,1,PETSC_DETERMINE,PETSC_DETERMINE));
  PetscCall(EPSSetWhichEigenpairs(eps,EPS_LARGEST_MAGNITUDE));
  PetscCall(EPSSetFromOptions(eps));
  PetscCall(EPSSolve(eps));
  PetscCall(EPSGetConverged(eps,&nconv));
  *rho = 0.0;
  if (nconv) {
    PetscCall(EPSGetEigenvalue(eps,0,&kr,&ki));
    *rho = SlepcAbsEigenvalue(kr,ki);
  }
  PetscCall(EPSDestroy(&eps));
  PetscFunctionReturn(PETSC_SUCCESS);
}

int main(int argc,char **argv)
{
  Graph          graph;           /* directed graph */
  Mat            A,B;             /* matrices */
  SVD            svd;             /* singular value problem solver context */
  PetscReal      beta=0.0,rho=0.0;
  char           filename[PETSC_MAX_PATH_LEN];
  PetscBool      flg,terse;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Load graph and build adjacency matrix
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(PetscOptionsGetString(NULL,NULL,"-graph",filename,sizeof(filename),&flg));
  PetscCheck(flg,PETSC_COMM_WORLD,PETSC_ERR_USER_INPUT,"Must use -graph <name> to indicate a file containing the graph in TSV format");

  PetscCall(GraphCreate(PETSC_COMM_WORLD,&graph));
  PetscCall(GraphPreload(graph,filename));
  PetscCheck(graph->type==GRAPH_DIRECTED,PETSC_COMM_SELF,PETSC_ERR_FILE_UNEXPECTED,"A directed graph is required");
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nGraph with %" PetscInt_FMT " vertices and %" PetscInt_FMT " edges.\n\n",graph->nvertices,graph->nedges));
  if (graph->weight!=GRAPH_WEIGHT_UNWEIGHTED) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"WARNING: ignoring weights in input graph\n"));

  PetscCall(GraphPreallocate(graph,filename));
  PetscCall(GraphLoadUnweighted(graph,filename));

  PetscCall(MatViewFromOptions(graph->adjacency,NULL,"-adjacency"));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Build high-order proximity matrices (Katz): A=I-beta*G, B=beta*G
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(PetscOptionsGetReal(NULL,NULL,"-beta",&beta,NULL));
  if (beta==0.0) {  /* compute decay parameter beta */
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Computing spectral radius...\n"));
    PetscCall(SpectralRadius(graph->adjacency,&rho));
    if (rho) {
      beta = 0.8*rho;
      PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Done, rho=%g, setting decay parameter beta=%g\n\n",(double)rho,(double)beta));
    } else {
      beta = 0.5;
      PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Failed, using default decay parameter beta=%g\n\n",(double)beta));
    }
  }
  PetscCall(MatDuplicate(graph->adjacency,MAT_COPY_VALUES,&A));
  PetscCall(MatScale(A,-beta));
  PetscCall(MatShift(A,1.0));
  PetscCall(MatDuplicate(graph->adjacency,MAT_COPY_VALUES,&B));
  PetscCall(MatScale(B,beta));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
           Create the singular value solver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(SVDCreate(PETSC_COMM_WORLD,&svd));
  PetscCall(SVDSetOperators(svd,A,B));
  PetscCall(SVDSetProblemType(svd,SVD_GENERALIZED));
  PetscCall(SVDSetType(svd,SVDTRLANCZOS));
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
  PetscCall(GraphDestroy(&graph));
  PetscCall(MatDestroy(&A));
  PetscCall(MatDestroy(&B));
  PetscCall(SlepcFinalize());
  return 0;
}
/*TEST

   build:
      depends: network.o

   test:
      args: -graph ${SLEPC_DIR}/share/slepc/datafiles/graphs/out.moreno_taro_taro -svd_nsv 4 -terse
      filter: sed -e 's/4.38031/4.38032/' | sed -e 's/4.38033/4.38032/' | sed -e 's/3.75089/3.7509/' | sed -e 's/3.00071/3.00072/'

TEST*/
