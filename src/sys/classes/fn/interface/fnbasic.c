/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Basic FN routines
*/

#include <slepc/private/fnimpl.h>      /*I "slepcfn.h" I*/
#include <slepcblaslapack.h>

PetscFunctionList FNList = NULL;
PetscBool         FNRegisterAllCalled = PETSC_FALSE;
PetscClassId      FN_CLASSID = 0;
PetscLogEvent     FN_Evaluate = 0;
static PetscBool  FNPackageInitialized = PETSC_FALSE;

const char *FNParallelTypes[] = {"REDUNDANT","SYNCHRONIZED","FNParallelType","FN_PARALLEL_",NULL};

/*@C
   FNFinalizePackage - This function destroys everything in the SLEPc interface
   to the `FN` package. It is called from `SlepcFinalize()`.

   Level: developer

.seealso: [](sec:fn), `SlepcFinalize()`, `FNInitializePackage()`
@*/
PetscErrorCode FNFinalizePackage(void)
{
  PetscFunctionBegin;
  PetscCall(PetscFunctionListDestroy(&FNList));
  FNPackageInitialized = PETSC_FALSE;
  FNRegisterAllCalled  = PETSC_FALSE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   FNInitializePackage - This function initializes everything in the `FN` package.
   It is called from `PetscDLLibraryRegister_slepc()` when using dynamic libraries, and
   on the first call to `FNCreate()` when using shared or static libraries.

   Note:
   This function never needs to be called by SLEPc users.

   Level: developer

.seealso: [](sec:fn), `FN`, `SlepcInitialize()`, `FNFinalizePackage()`
@*/
PetscErrorCode FNInitializePackage(void)
{
  char           logList[256];
  PetscBool      opt,pkg;
  PetscClassId   classids[1];

  PetscFunctionBegin;
  if (FNPackageInitialized) PetscFunctionReturn(PETSC_SUCCESS);
  FNPackageInitialized = PETSC_TRUE;
  /* Register Classes */
  PetscCall(PetscClassIdRegister("Math Function",&FN_CLASSID));
  /* Register Constructors */
  PetscCall(FNRegisterAll());
  /* Register Events */
  PetscCall(PetscLogEventRegister("FNEvaluate",FN_CLASSID,&FN_Evaluate));
  /* Process Info */
  classids[0] = FN_CLASSID;
  PetscCall(PetscInfoProcessClass("fn",1,&classids[0]));
  /* Process summary exclusions */
  PetscCall(PetscOptionsGetString(NULL,NULL,"-log_exclude",logList,sizeof(logList),&opt));
  if (opt) {
    PetscCall(PetscStrInList("fn",logList,',',&pkg));
    if (pkg) PetscCall(PetscLogEventDeactivateClass(FN_CLASSID));
  }
  /* Register package finalizer */
  PetscCall(PetscRegisterFinalize(FNFinalizePackage));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   FNCreate - Creates an `FN` context.

   Collective

   Input Parameter:
.  comm - MPI communicator

   Output Parameter:
.  newfn - location to put the `FN` context

   Level: beginner

.seealso: [](sec:fn), `FNDestroy()`, `FN`
@*/
PetscErrorCode FNCreate(MPI_Comm comm,FN *newfn)
{
  FN             fn;

  PetscFunctionBegin;
  PetscAssertPointer(newfn,2);
  PetscCall(FNInitializePackage());
  PetscCall(SlepcHeaderCreate(fn,FN_CLASSID,"FN","Math Function","FN",comm,FNDestroy,FNView));

  fn->alpha    = 1.0;
  fn->beta     = 1.0;
  fn->method   = 0;

  fn->nw       = 0;
  fn->cw       = 0;
  fn->data     = NULL;

  *newfn = fn;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   FNSetOptionsPrefix - Sets the prefix used for searching for all
   `FN` options in the database.

   Logically Collective

   Input Parameters:
+  fn - the math function context
-  prefix - the prefix string to prepend to all `FN` option requests

   Notes:
   A hyphen (-) must NOT be given at the beginning of the prefix name.
   The first character of all runtime options is AUTOMATICALLY the
   hyphen.

   Level: advanced

.seealso: [](sec:fn), `FNAppendOptionsPrefix()`
@*/
PetscErrorCode FNSetOptionsPrefix(FN fn,const char prefix[])
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(fn,FN_CLASSID,1);
  PetscCall(PetscObjectSetOptionsPrefix((PetscObject)fn,prefix));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   FNAppendOptionsPrefix - Appends to the prefix used for searching for all
   `FN` options in the database.

   Logically Collective

   Input Parameters:
+  fn - the math function context
-  prefix - the prefix string to prepend to all `FN` option requests

   Notes:
   A hyphen (-) must NOT be given at the beginning of the prefix name.
   The first character of all runtime options is AUTOMATICALLY the hyphen.

   Level: advanced

.seealso: [](sec:fn), `FNSetOptionsPrefix()`
@*/
PetscErrorCode FNAppendOptionsPrefix(FN fn,const char prefix[])
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(fn,FN_CLASSID,1);
  PetscCall(PetscObjectAppendOptionsPrefix((PetscObject)fn,prefix));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   FNGetOptionsPrefix - Gets the prefix used for searching for all
   `FN` options in the database.

   Not Collective

   Input Parameter:
.  fn - the math function context

   Output Parameter:
.  prefix - pointer to the prefix string used is returned

   Level: advanced

.seealso: [](sec:fn), `FNSetOptionsPrefix()`, `FNAppendOptionsPrefix()`
@*/
PetscErrorCode FNGetOptionsPrefix(FN fn,const char *prefix[])
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(fn,FN_CLASSID,1);
  PetscAssertPointer(prefix,2);
  PetscCall(PetscObjectGetOptionsPrefix((PetscObject)fn,prefix));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   FNSetType - Selects the type for the `FN` object.

   Logically Collective

   Input Parameters:
+  fn   - the math function context
-  type - a known type

   Options Database Key:
.  -fn_type \<type\> - sets the `FN` type

   Note:
   The default is `FNRATIONAL`, which includes polynomials as a particular
   case as well as simple functions such as $f(x)=x$ and $f(x)=constant$.

   Level: intermediate

.seealso: [](sec:fn), `FNGetType()`
@*/
PetscErrorCode FNSetType(FN fn,FNType type)
{
  PetscErrorCode (*r)(FN);
  PetscBool      match;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(fn,FN_CLASSID,1);
  PetscAssertPointer(type,2);

  PetscCall(PetscObjectTypeCompare((PetscObject)fn,type,&match));
  if (match) PetscFunctionReturn(PETSC_SUCCESS);

  PetscCall(PetscFunctionListFind(FNList,type,&r));
  PetscCheck(r,PetscObjectComm((PetscObject)fn),PETSC_ERR_ARG_UNKNOWN_TYPE,"Unable to find requested FN type %s",type);

  PetscTryTypeMethod(fn,destroy);
  PetscCall(PetscMemzero(fn->ops,sizeof(struct _FNOps)));

  PetscCall(PetscObjectChangeTypeName((PetscObject)fn,type));
  PetscCall((*r)(fn));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   FNGetType - Gets the `FN` type name (as a string) from the `FN` context.

   Not Collective

   Input Parameter:
.  fn - the math function context

   Output Parameter:
.  type - name of the math function

   Level: intermediate

.seealso: [](sec:fn), `FNSetType()`
@*/
PetscErrorCode FNGetType(FN fn,FNType *type)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(fn,FN_CLASSID,1);
  PetscAssertPointer(type,2);
  *type = ((PetscObject)fn)->type_name;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   FNSetScale - Sets the scaling parameters that define the matematical function.

   Logically Collective

   Input Parameters:
+  fn    - the math function context
.  alpha - inner scaling (argument)
-  beta  - outer scaling (result)

   Notes:
   Given a function $f(x)$ specified by the `FN` type, the scaling parameters can
   be used to realize the function $\beta f(\alpha x)$. So when these values are given,
   the procedure for function evaluation will first multiply the argument by $\alpha$,
   then evaluate the function itself, and finally scale the result by $\beta$.
   Likewise, these values are also considered when evaluating the derivative.

   If you want to provide only one of the two scaling factors, set the other
   one to 1.0.

   Level: intermediate

.seealso: [](sec:fn), `FNGetScale()`, `FNEvaluateFunction()`
@*/
PetscErrorCode FNSetScale(FN fn,PetscScalar alpha,PetscScalar beta)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(fn,FN_CLASSID,1);
  PetscValidLogicalCollectiveScalar(fn,alpha,2);
  PetscValidLogicalCollectiveScalar(fn,beta,3);
  PetscCheck(PetscAbsScalar(alpha)!=0.0 && PetscAbsScalar(beta)!=0.0,PetscObjectComm((PetscObject)fn),PETSC_ERR_ARG_WRONG,"Scaling factors must be nonzero");
  fn->alpha = alpha;
  fn->beta  = beta;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   FNGetScale - Gets the scaling parameters that define the matematical function.

   Not Collective

   Input Parameter:
.  fn    - the math function context

   Output Parameters:
+  alpha - inner scaling (argument)
-  beta  - outer scaling (result)

   Level: intermediate

.seealso: [](sec:fn), `FNSetScale()`
@*/
PetscErrorCode FNGetScale(FN fn,PetscScalar *alpha,PetscScalar *beta)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(fn,FN_CLASSID,1);
  if (alpha) *alpha = fn->alpha;
  if (beta)  *beta  = fn->beta;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   FNSetMethod - Selects the method to be used to evaluate functions of matrices.

   Logically Collective

   Input Parameters:
+  fn   - the math function context
-  meth - an index identifying the method

   Options Database Key:
.  -fn_method \<meth\> - sets the method

   Notes:
   In some `FN` types there are more than one algorithm available for computing
   matrix functions. In that case, this function allows choosing the wanted method.

   If `meth` is currently set to 0 (the default) and the input argument `A` of
   `FNEvaluateFunctionMat()` is a symmetric/Hermitian matrix, then the computation
   is done via the eigendecomposition of `A`, rather than with the general algorithm.

   Level: intermediate

.seealso: [](sec:fn), `FNGetMethod()`, `FNEvaluateFunctionMat()`
@*/
PetscErrorCode FNSetMethod(FN fn,PetscInt meth)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(fn,FN_CLASSID,1);
  PetscValidLogicalCollectiveInt(fn,meth,2);
  PetscCheck(meth>=0,PetscObjectComm((PetscObject)fn),PETSC_ERR_ARG_OUTOFRANGE,"The method must be a non-negative integer");
  PetscCheck(meth<=FN_MAX_SOLVE,PetscObjectComm((PetscObject)fn),PETSC_ERR_ARG_OUTOFRANGE,"Too large value for the method");
  fn->method = meth;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   FNGetMethod - Gets the method currently used in the `FN`.

   Not Collective

   Input Parameter:
.  fn - the math function context

   Output Parameter:
.  meth - identifier of the method

   Level: intermediate

.seealso: [](sec:fn), `FNSetMethod()`
@*/
PetscErrorCode FNGetMethod(FN fn,PetscInt *meth)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(fn,FN_CLASSID,1);
  PetscAssertPointer(meth,2);
  *meth = fn->method;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   FNSetParallel - Selects the mode of operation in parallel runs.

   Logically Collective

   Input Parameters:
+  fn    - the math function context
-  pmode - the parallel mode

   Options Database Key:
.  -fn_parallel \<pmode\> - sets the parallel mode, either `redundant` or `synchronized`

   Notes:
   This is relevant only when the function is evaluated on a matrix, with
   either `FNEvaluateFunctionMat()` or `FNEvaluateFunctionMatVec()`.

   In the `redundant` parallel mode, all processes will make the computation
   redundantly, starting from the same data, and producing the same result.
   This result may be slightly different in the different processes if using a
   multithreaded BLAS library, which may cause issues in ill-conditioned problems.

   In the `synchronized` parallel mode, only the first MPI process performs the
   computation and then the computed matrix is broadcast to the other
   processes in the communicator. This communication is done automatically at
   the end of `FNEvaluateFunctionMat()` or `FNEvaluateFunctionMatVec()`.

   Level: advanced

.seealso: [](sec:fn), `FNEvaluateFunctionMat()` or `FNEvaluateFunctionMatVec()`, `FNGetParallel()`
@*/
PetscErrorCode FNSetParallel(FN fn,FNParallelType pmode)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(fn,FN_CLASSID,1);
  PetscValidLogicalCollectiveEnum(fn,pmode,2);
  fn->pmode = pmode;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   FNGetParallel - Gets the mode of operation in parallel runs.

   Not Collective

   Input Parameter:
.  fn - the math function context

   Output Parameter:
.  pmode - the parallel mode

   Level: advanced

.seealso: [](sec:fn), `FNSetParallel()`
@*/
PetscErrorCode FNGetParallel(FN fn,FNParallelType *pmode)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(fn,FN_CLASSID,1);
  PetscAssertPointer(pmode,2);
  *pmode = fn->pmode;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   FNEvaluateFunction - Computes the value of the function $f(x)$ for a given $x$.

   Not Collective

   Input Parameters:
+  fn - the math function context
-  x  - the value where the function must be evaluated

   Output Parameter:
.  y  - the result of $f(x)$

   Note:
   Scaling factors are taken into account, so the actual function evaluation
   will return $\beta f(\alpha x)$.

   Level: intermediate

.seealso: [](sec:fn), `FNEvaluateDerivative()`, `FNEvaluateFunctionMat()`, `FNSetScale()`
@*/
PetscErrorCode FNEvaluateFunction(FN fn,PetscScalar x,PetscScalar *y)
{
  PetscScalar    xf,yf;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(fn,FN_CLASSID,1);
  PetscValidType(fn,1);
  PetscAssertPointer(y,3);
  PetscCall(PetscLogEventBegin(FN_Evaluate,fn,0,0,0));
  xf = fn->alpha*x;
  PetscUseTypeMethod(fn,evaluatefunction,xf,&yf);
  *y = fn->beta*yf;
  PetscCall(PetscLogEventEnd(FN_Evaluate,fn,0,0,0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   FNEvaluateDerivative - Computes the value of the derivative $f'(x)$ for a given $x$.

   Not Collective

   Input Parameters:
+  fn - the math function context
-  x  - the value where the derivative must be evaluated

   Output Parameter:
.  y  - the result of $f'(x)$

   Note:
   Scaling factors are taken into account, so the actual derivative evaluation will
   return $\alpha\beta f'(\alpha x)$.

   Level: intermediate

.seealso: [](sec:fn), `FNEvaluateFunction()`, `FNSetScale()`
@*/
PetscErrorCode FNEvaluateDerivative(FN fn,PetscScalar x,PetscScalar *y)
{
  PetscScalar    xf,yf;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(fn,FN_CLASSID,1);
  PetscValidType(fn,1);
  PetscAssertPointer(y,3);
  PetscCall(PetscLogEventBegin(FN_Evaluate,fn,0,0,0));
  xf = fn->alpha*x;
  PetscUseTypeMethod(fn,evaluatederivative,xf,&yf);
  *y = fn->alpha*fn->beta*yf;
  PetscCall(PetscLogEventEnd(FN_Evaluate,fn,0,0,0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode FNEvaluateFunctionMat_Sym_Private(FN fn,const PetscScalar *As,PetscScalar *Bs,PetscInt m,PetscBool firstonly)
{
  PetscInt       i,j;
  PetscBLASInt   n,k,ld,lwork,info;
  PetscScalar    *Q,*W,*work,adummy,a,x,y,one=1.0,zero=0.0;
  PetscReal      *eig,dummy;
#if defined(PETSC_USE_COMPLEX)
  PetscReal      *rwork,rdummy;
#endif

  PetscFunctionBegin;
  PetscCall(PetscBLASIntCast(m,&n));
  ld = n;
  k = firstonly? 1: n;

  /* workspace query and memory allocation */
  lwork = -1;
#if defined(PETSC_USE_COMPLEX)
  PetscCallBLAS("LAPACKsyev",LAPACKsyev_("V","L",&n,&adummy,&ld,&dummy,&a,&lwork,&rdummy,&info));
  PetscCall(PetscBLASIntCast((PetscInt)PetscRealPart(a),&lwork));
  PetscCall(PetscMalloc5(m,&eig,m*m,&Q,m*k,&W,lwork,&work,PetscMax(1,3*m-2),&rwork));
#else
  PetscCallBLAS("LAPACKsyev",LAPACKsyev_("V","L",&n,&adummy,&ld,&dummy,&a,&lwork,&info));
  PetscCall(PetscBLASIntCast((PetscInt)a,&lwork));
  PetscCall(PetscMalloc4(m,&eig,m*m,&Q,m*k,&W,lwork,&work));
#endif

  /* compute eigendecomposition */
  for (j=0;j<n;j++) for (i=j;i<n;i++) Q[i+j*ld] = As[i+j*ld];
#if defined(PETSC_USE_COMPLEX)
  PetscCallBLAS("LAPACKsyev",LAPACKsyev_("V","L",&n,Q,&ld,eig,work,&lwork,rwork,&info));
#else
  PetscCallBLAS("LAPACKsyev",LAPACKsyev_("V","L",&n,Q,&ld,eig,work,&lwork,&info));
#endif
  SlepcCheckLapackInfo("syev",info);

  /* W = f(Lambda)*Q' */
  for (i=0;i<n;i++) {
    x = fn->alpha*eig[i];
    PetscUseTypeMethod(fn,evaluatefunction,x,&y);  /* y = f(x) */
    for (j=0;j<k;j++) W[i+j*ld] = PetscConj(Q[j+i*ld])*fn->beta*y;
  }
  /* Bs = Q*W */
  PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&n,&k,&n,&one,Q,&ld,W,&ld,&zero,Bs,&ld));
#if defined(PETSC_USE_COMPLEX)
  PetscCall(PetscFree5(eig,Q,W,work,rwork));
#else
  PetscCall(PetscFree4(eig,Q,W,work));
#endif
  PetscCall(PetscLogFlops(9.0*n*n*n+2.0*n*n*n));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   FNEvaluateFunctionMat_Sym_Default - given a symmetric matrix A,
   compute the matrix function as f(A)=Q*f(D)*Q' where the spectral
   decomposition of A is A=Q*D*Q'
*/
static PetscErrorCode FNEvaluateFunctionMat_Sym_Default(FN fn,Mat A,Mat B)
{
  PetscInt          m;
  const PetscScalar *As;
  PetscScalar       *Bs;

  PetscFunctionBegin;
  PetscCall(MatDenseGetArrayRead(A,&As));
  PetscCall(MatDenseGetArray(B,&Bs));
  PetscCall(MatGetSize(A,&m,NULL));
  PetscCall(FNEvaluateFunctionMat_Sym_Private(fn,As,Bs,m,PETSC_FALSE));
  PetscCall(MatDenseRestoreArrayRead(A,&As));
  PetscCall(MatDenseRestoreArray(B,&Bs));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode FNEvaluateFunctionMat_Basic(FN fn,Mat A,Mat F)
{
  PetscBool iscuda;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)A,MATSEQDENSECUDA,&iscuda));
  if (iscuda && !fn->ops->evaluatefunctionmatcuda[fn->method]) PetscCall(PetscInfo(fn,"The method %" PetscInt_FMT " is not implemented for CUDA, falling back to CPU version\n",fn->method));
  if (iscuda && fn->ops->evaluatefunctionmatcuda[fn->method]) PetscUseTypeMethod(fn,evaluatefunctionmatcuda[fn->method],A,F);
  else if (fn->ops->evaluatefunctionmat[fn->method]) PetscUseTypeMethod(fn,evaluatefunctionmat[fn->method],A,F);
  else {
    PetscCheck(fn->method,PetscObjectComm((PetscObject)fn),PETSC_ERR_SUP,"Matrix functions not implemented in this FN type");
    PetscCheck(!fn->method,PetscObjectComm((PetscObject)fn),PETSC_ERR_ARG_OUTOFRANGE,"The specified method number does not exist for this FN type");
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode FNEvaluateFunctionMat_Private(FN fn,Mat A,Mat B,PetscBool sync)
{
  PetscBool      set,flg,symm=PETSC_FALSE,iscuda,hasspecificmeth;
  PetscInt       m,n;
  PetscMPIInt    size,rank,n2;
  PetscScalar    *pF;
  Mat            M,F;

  PetscFunctionBegin;
  /* destination matrix */
  F = B?B:A;

  /* check symmetry of A */
  PetscCall(MatIsHermitianKnown(A,&set,&flg));
  symm = set? flg: PETSC_FALSE;

  PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)fn),&size));
  PetscCallMPI(MPI_Comm_rank(PetscObjectComm((PetscObject)fn),&rank));
  if (size==1 || fn->pmode==FN_PARALLEL_REDUNDANT || (fn->pmode==FN_PARALLEL_SYNCHRONIZED && !rank)) {
    PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
    PetscCall(PetscObjectTypeCompare((PetscObject)A,MATSEQDENSECUDA,&iscuda));
    hasspecificmeth = ((iscuda && fn->ops->evaluatefunctionmatcuda[fn->method]) || (!iscuda && fn->method && fn->ops->evaluatefunctionmat[fn->method]))? PETSC_TRUE: PETSC_FALSE;
    if (!hasspecificmeth && symm && !fn->method) {  /* prefer diagonalization */
      PetscCall(PetscInfo(fn,"Computing matrix function via diagonalization\n"));
      PetscCall(FNEvaluateFunctionMat_Sym_Default(fn,A,F));
    } else {
      /* scale argument */
      if (fn->alpha!=(PetscScalar)1.0) {
        PetscCall(FN_AllocateWorkMat(fn,A,&M));
        PetscCall(MatScale(M,fn->alpha));
      } else M = A;
      PetscCall(FNEvaluateFunctionMat_Basic(fn,M,F));
      if (fn->alpha!=(PetscScalar)1.0) PetscCall(FN_FreeWorkMat(fn,&M));
      /* scale result */
      PetscCall(MatScale(F,fn->beta));
    }
    PetscCall(PetscFPTrapPop());
  }
  if (size>1 && fn->pmode==FN_PARALLEL_SYNCHRONIZED && sync) {  /* synchronize */
    PetscCall(MatGetSize(A,&m,&n));
    PetscCall(MatDenseGetArray(F,&pF));
    PetscCall(PetscMPIIntCast(n*n,&n2));
    PetscCallMPI(MPI_Bcast(pF,n2,MPIU_SCALAR,0,PetscObjectComm((PetscObject)fn)));
    PetscCall(MatDenseRestoreArray(F,&pF));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   FNEvaluateFunctionMat - Computes the value of the function $f(A)$ for a given
   matrix $A$, where the result is also a matrix.

   Logically Collective

   Input Parameters:
+  fn - the math function context
-  A  - matrix on which the function must be evaluated

   Output Parameter:
.  B  - (optional) matrix resulting from evaluating $f(A)$

   Notes:
   Matrix `A` must be a square sequential dense `Mat`, with all entries equal on
   all processes (otherwise each process will compute different results).
   If matrix `B` is provided, it must also be a square sequential dense `Mat`, and
   both matrices must have the same dimensions. If `B` is `NULL` (or `B`=`A`) then
   the function will perform an in-place computation, overwriting `A` with $f(A)$.

   If `A` is known to be real symmetric or complex Hermitian then it is
   recommended to set the appropriate flag with `MatSetOption()`, because
   symmetry can sometimes be exploited by the algorithm.

   Scaling factors are taken into account, so the actual function evaluation
   will return $\beta f(\alpha A)$.

   Level: advanced

.seealso: [](sec:fn), `FNEvaluateFunction()`, `FNEvaluateFunctionMatVec()`, `FNSetMethod()`
@*/
PetscErrorCode FNEvaluateFunctionMat(FN fn,Mat A,Mat B)
{
  PetscBool      inplace=PETSC_FALSE;
  PetscInt       m,n,n1;
  MatType        type;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(fn,FN_CLASSID,1);
  PetscValidHeaderSpecific(A,MAT_CLASSID,2);
  PetscValidType(fn,1);
  PetscValidType(A,2);
  if (B) {
    PetscValidHeaderSpecific(B,MAT_CLASSID,3);
    PetscValidType(B,3);
  } else inplace = PETSC_TRUE;
  PetscCheckTypeNames(A,MATSEQDENSE,MATSEQDENSECUDA);   //SlepcMatCheckSeq(A);
  PetscCall(MatGetSize(A,&m,&n));
  PetscCheck(m==n,PetscObjectComm((PetscObject)fn),PETSC_ERR_ARG_SIZ,"Mat A is not square (has %" PetscInt_FMT " rows, %" PetscInt_FMT " cols)",m,n);
  if (!inplace) {
    PetscCall(MatGetType(A,&type));
    PetscCheckTypeName(B,type);
    n1 = n;
    PetscCall(MatGetSize(B,&m,&n));
    PetscCheck(m==n,PetscObjectComm((PetscObject)fn),PETSC_ERR_ARG_SIZ,"Mat B is not square (has %" PetscInt_FMT " rows, %" PetscInt_FMT " cols)",m,n);
    PetscCheck(n1==n,PetscObjectComm((PetscObject)fn),PETSC_ERR_ARG_SIZ,"Matrices A and B must have the same dimension");
  }

  /* evaluate matrix function */
  PetscCall(PetscLogEventBegin(FN_Evaluate,fn,0,0,0));
  PetscCall(FNEvaluateFunctionMat_Private(fn,A,B,PETSC_TRUE));
  PetscCall(PetscLogEventEnd(FN_Evaluate,fn,0,0,0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   FNEvaluateFunctionMatVec_Default - computes the full matrix f(A)
   and then copies the first column.
*/
static PetscErrorCode FNEvaluateFunctionMatVec_Default(FN fn,Mat A,Vec v)
{
  Mat            F;

  PetscFunctionBegin;
  PetscCall(FN_AllocateWorkMat(fn,A,&F));
  PetscCall(FNEvaluateFunctionMat_Basic(fn,A,F));
  PetscCall(MatGetColumnVector(F,v,0));
  PetscCall(FN_FreeWorkMat(fn,&F));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   FNEvaluateFunctionMatVec_Sym_Default - given a symmetric matrix A,
   compute the matrix function as f(A)=Q*f(D)*Q' where the spectral
   decomposition of A is A=Q*D*Q'. Only the first column is computed.
*/
static PetscErrorCode FNEvaluateFunctionMatVec_Sym_Default(FN fn,Mat A,Vec v)
{
  PetscInt          m;
  const PetscScalar *As;
  PetscScalar       *vs;

  PetscFunctionBegin;
  PetscCall(MatDenseGetArrayRead(A,&As));
  PetscCall(VecGetArray(v,&vs));
  PetscCall(MatGetSize(A,&m,NULL));
  PetscCall(FNEvaluateFunctionMat_Sym_Private(fn,As,vs,m,PETSC_TRUE));
  PetscCall(MatDenseRestoreArrayRead(A,&As));
  PetscCall(VecRestoreArray(v,&vs));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode FNEvaluateFunctionMatVec_Private(FN fn,Mat A,Vec v,PetscBool sync)
{
  PetscBool      set,flg,symm=PETSC_FALSE,iscuda,hasspecificmeth;
  PetscInt       m,n;
  Mat            M;
  PetscMPIInt    size,rank,n_;
  PetscScalar    *pv;

  PetscFunctionBegin;
  /* check symmetry of A */
  PetscCall(MatIsHermitianKnown(A,&set,&flg));
  symm = set? flg: PETSC_FALSE;

  /* evaluate matrix function */
  PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)fn),&size));
  PetscCallMPI(MPI_Comm_rank(PetscObjectComm((PetscObject)fn),&rank));
  if (size==1 || fn->pmode==FN_PARALLEL_REDUNDANT || (fn->pmode==FN_PARALLEL_SYNCHRONIZED && !rank)) {
    PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
    PetscCall(PetscObjectTypeCompare((PetscObject)A,MATSEQDENSECUDA,&iscuda));
    hasspecificmeth = ((iscuda && fn->ops->evaluatefunctionmatcuda[fn->method]) || (!iscuda && fn->method && fn->ops->evaluatefunctionmat[fn->method]))? PETSC_TRUE: PETSC_FALSE;
    if (!hasspecificmeth && symm && !fn->method) {  /* prefer diagonalization */
      PetscCall(PetscInfo(fn,"Computing matrix function via diagonalization\n"));
      PetscCall(FNEvaluateFunctionMatVec_Sym_Default(fn,A,v));
    } else {
      /* scale argument */
      if (fn->alpha!=(PetscScalar)1.0) {
        PetscCall(FN_AllocateWorkMat(fn,A,&M));
        PetscCall(MatScale(M,fn->alpha));
      } else M = A;
      if (iscuda && fn->ops->evaluatefunctionmatveccuda[fn->method]) PetscUseTypeMethod(fn,evaluatefunctionmatveccuda[fn->method],M,v);
      else if (fn->ops->evaluatefunctionmatvec[fn->method]) PetscUseTypeMethod(fn,evaluatefunctionmatvec[fn->method],M,v);
      else PetscCall(FNEvaluateFunctionMatVec_Default(fn,M,v));
      if (fn->alpha!=(PetscScalar)1.0) PetscCall(FN_FreeWorkMat(fn,&M));
      /* scale result */
      PetscCall(VecScale(v,fn->beta));
    }
    PetscCall(PetscFPTrapPop());
  }

  /* synchronize */
  if (size>1 && fn->pmode==FN_PARALLEL_SYNCHRONIZED && sync) {
    PetscCall(MatGetSize(A,&m,&n));
    PetscCall(VecGetArray(v,&pv));
    PetscCall(PetscMPIIntCast(n,&n_));
    PetscCallMPI(MPI_Bcast(pv,n_,MPIU_SCALAR,0,PetscObjectComm((PetscObject)fn)));
    PetscCall(VecRestoreArray(v,&pv));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   FNEvaluateFunctionMatVec - Computes the first column of the matrix $f(A)$
   for a given matrix $A$.

   Logically Collective

   Input Parameters:
+  fn - the math function context
-  A  - matrix on which the function must be evaluated

   Output Parameter:
.  v  - vector to hold the first column of $f(A)$

   Notes:
   This operation is similar to `FNEvaluateFunctionMat()` but returns only
   the first column of $f(A)$, hence saving computations in most cases.

   Level: advanced

.seealso: [](sec:fn), `FNEvaluateFunction()`, `FNEvaluateFunctionMat()`, `FNSetMethod()`
@*/
PetscErrorCode FNEvaluateFunctionMatVec(FN fn,Mat A,Vec v)
{
  PetscInt       m,n;
  PetscBool      iscuda;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(fn,FN_CLASSID,1);
  PetscValidHeaderSpecific(A,MAT_CLASSID,2);
  PetscValidHeaderSpecific(v,VEC_CLASSID,3);
  PetscValidType(fn,1);
  PetscValidType(A,2);
  PetscValidType(v,3);
  PetscCheckTypeNames(A,MATSEQDENSE,MATSEQDENSECUDA);   //SlepcMatCheckSeq(A);
  PetscCall(MatGetSize(A,&m,&n));
  PetscCheck(m==n,PetscObjectComm((PetscObject)fn),PETSC_ERR_ARG_SIZ,"Mat A is not square (has %" PetscInt_FMT " rows, %" PetscInt_FMT " cols)",m,n);
  PetscCall(PetscObjectTypeCompare((PetscObject)A,MATSEQDENSECUDA,&iscuda));
  PetscCheckTypeName(v,iscuda?VECSEQCUDA:VECSEQ);
  PetscCall(VecGetSize(v,&m));
  PetscCheck(m==n,PetscObjectComm((PetscObject)fn),PETSC_ERR_ARG_SIZ,"Matrix A and vector v must have the same size");
  PetscCall(PetscLogEventBegin(FN_Evaluate,fn,0,0,0));
  PetscCall(FNEvaluateFunctionMatVec_Private(fn,A,v,PETSC_TRUE));
  PetscCall(PetscLogEventEnd(FN_Evaluate,fn,0,0,0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   FNSetFromOptions - Sets `FN` options from the options database.

   Collective

   Input Parameter:
.  fn - the math function context

   Note:
   To see all options, run your program with the `-help` option.

   Level: beginner

.seealso: [](sec:fn), `FNSetOptionsPrefix()`
@*/
PetscErrorCode FNSetFromOptions(FN fn)
{
  char           type[256];
  PetscScalar    array[2];
  PetscInt       k,meth;
  PetscBool      flg;
  FNParallelType pmode;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(fn,FN_CLASSID,1);
  PetscCall(FNRegisterAll());
  PetscObjectOptionsBegin((PetscObject)fn);
    PetscCall(PetscOptionsFList("-fn_type","Math function type","FNSetType",FNList,(char*)(((PetscObject)fn)->type_name?((PetscObject)fn)->type_name:FNRATIONAL),type,sizeof(type),&flg));
    if (flg) PetscCall(FNSetType(fn,type));
    else if (!((PetscObject)fn)->type_name) PetscCall(FNSetType(fn,FNRATIONAL));

    k = 2;
    array[0] = 0.0; array[1] = 0.0;
    PetscCall(PetscOptionsScalarArray("-fn_scale","Scale factors (one or two scalar values separated with a comma without spaces)","FNSetScale",array,&k,&flg));
    if (flg) {
      if (k<2) array[1] = 1.0;
      PetscCall(FNSetScale(fn,array[0],array[1]));
    }

    PetscCall(PetscOptionsInt("-fn_method","Method to be used for computing matrix functions","FNSetMethod",fn->method,&meth,&flg));
    if (flg) PetscCall(FNSetMethod(fn,meth));

    PetscCall(PetscOptionsEnum("-fn_parallel","Operation mode in parallel runs","FNSetParallel",FNParallelTypes,(PetscEnum)fn->pmode,(PetscEnum*)&pmode,&flg));
    if (flg) PetscCall(FNSetParallel(fn,pmode));

    PetscTryTypeMethod(fn,setfromoptions,PetscOptionsObject);
    PetscCall(PetscObjectProcessOptionsHandlers((PetscObject)fn,PetscOptionsObject));
  PetscOptionsEnd();
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   FNView - Prints the `FN` data structure.

   Collective

   Input Parameters:
+  fn - the math function context
-  viewer - optional visualization context

   Note:
   The available visualization contexts include
+     `PETSC_VIEWER_STDOUT_SELF` - standard output (default)
-     `PETSC_VIEWER_STDOUT_WORLD` - synchronized standard output where only the
         first process opens the file; all other processes send their data to the
         first one to print

   The user can open an alternative visualization context with `PetscViewerASCIIOpen()`
   to output to a specified file.

   Use `FNViewFromOptions()` to allow the user to select many different `PetscViewerType`
   and formats from the options database.

   Level: beginner

.seealso: [](sec:fn), `FNCreate()`, `FNViewFromOptions()`
@*/
PetscErrorCode FNView(FN fn,PetscViewer viewer)
{
  PetscBool      isascii;
  PetscMPIInt    size;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(fn,FN_CLASSID,1);
  if (!viewer) PetscCall(PetscViewerASCIIGetStdout(PetscObjectComm((PetscObject)fn),&viewer));
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,2);
  PetscCheckSameComm(fn,1,viewer,2);
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) {
    PetscCall(PetscObjectPrintClassNamePrefixType((PetscObject)fn,viewer));
    PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)fn),&size));
    if (size>1) PetscCall(PetscViewerASCIIPrintf(viewer,"  parallel operation mode: %s\n",FNParallelTypes[fn->pmode]));
    PetscCall(PetscViewerASCIIPushTab(viewer));
    PetscTryTypeMethod(fn,view,viewer);
    PetscCall(PetscViewerASCIIPopTab(viewer));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   FNViewFromOptions - View (print) an `FN` object based on values in the options database.

   Collective

   Input Parameters:
+  fn   - the math function context
.  obj  - optional object that provides the options prefix used to query the options database
-  name - command line option

   Level: intermediate

.seealso: [](sec:fn), `FNView()`, `FNCreate()`, `PetscObjectViewFromOptions()`
@*/
PetscErrorCode FNViewFromOptions(FN fn,PetscObject obj,const char name[])
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(fn,FN_CLASSID,1);
  PetscCall(PetscObjectViewFromOptions((PetscObject)fn,obj,name));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   FNDuplicate - Duplicates a math function, copying all parameters, possibly with a
   different communicator.

   Collective

   Input Parameters:
+  fn   - the math function context
-  comm - MPI communicator

   Output Parameter:
.  newfn - location to put the new `FN` context

   Note:
   In order to use the same MPI communicator as in the original object,
   use `PetscObjectComm`((`PetscObject`)`fn`).

   Level: developer

.seealso: [](sec:fn), `FNCreate()`
@*/
PetscErrorCode FNDuplicate(FN fn,MPI_Comm comm,FN *newfn)
{
  FNType         type;
  PetscScalar    alpha,beta;
  PetscInt       meth;
  FNParallelType ptype;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(fn,FN_CLASSID,1);
  PetscValidType(fn,1);
  PetscAssertPointer(newfn,3);
  PetscCall(FNCreate(comm,newfn));
  PetscCall(FNGetType(fn,&type));
  PetscCall(FNSetType(*newfn,type));
  PetscCall(FNGetScale(fn,&alpha,&beta));
  PetscCall(FNSetScale(*newfn,alpha,beta));
  PetscCall(FNGetMethod(fn,&meth));
  PetscCall(FNSetMethod(*newfn,meth));
  PetscCall(FNGetParallel(fn,&ptype));
  PetscCall(FNSetParallel(*newfn,ptype));
  PetscTryTypeMethod(fn,duplicate,comm,newfn);
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   FNDestroy - Destroys an `FN` context that was created with `FNCreate()`.

   Collective

   Input Parameter:
.  fn - the math function context

   Level: beginner

.seealso: [](sec:fn), `FNCreate()`
@*/
PetscErrorCode FNDestroy(FN *fn)
{
  PetscInt       i;

  PetscFunctionBegin;
  if (!*fn) PetscFunctionReturn(PETSC_SUCCESS);
  PetscValidHeaderSpecific(*fn,FN_CLASSID,1);
  if (--((PetscObject)*fn)->refct > 0) { *fn = NULL; PetscFunctionReturn(PETSC_SUCCESS); }
  PetscTryTypeMethod(*fn,destroy);
  for (i=0;i<(*fn)->nw;i++) PetscCall(MatDestroy(&(*fn)->W[i]));
  PetscCall(PetscHeaderDestroy(fn));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   FNRegister - Adds a mathematical function to the `FN` package.

   Not Collective

   Input Parameters:
+  name - name of a new user-defined `FN`
-  function - routine to create the context

   Notes:
   `FNRegister()` may be called multiple times to add several user-defined functions.

   Level: advanced

.seealso: [](sec:fn), `FNRegisterAll()`
@*/
PetscErrorCode FNRegister(const char *name,PetscErrorCode (*function)(FN))
{
  PetscFunctionBegin;
  PetscCall(FNInitializePackage());
  PetscCall(PetscFunctionListAdd(&FNList,name,function));
  PetscFunctionReturn(PETSC_SUCCESS);
}
