/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   A polynomial filter based on a truncated Chebyshev series.
*/

#include <slepc/private/stimpl.h>
#include "filter.h"

static PetscErrorCode STComputeOperator_Filter_Chebyshev(ST st,Mat *G)
{
  PetscFunctionBegin;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode STDestroy_Filter_Chebyshev(ST st)
{
  ST_FILTER     *ctx = (ST_FILTER*)st->data;
  CHEBYSHEV_CTX cheby = (CHEBYSHEV_CTX)ctx->data;

  PetscFunctionBegin;
  PetscCall(PetscFree(cheby));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode STCreate_Filter_Chebyshev(ST st)
{
  ST_FILTER     *ctx = (ST_FILTER*)st->data;
  CHEBYSHEV_CTX cheby;

  PetscFunctionBegin;
  PetscCall(PetscNew(&cheby));
  ctx->data = (void*)cheby;

  ctx->computeoperator = STComputeOperator_Filter_Chebyshev;
  ctx->destroy         = STDestroy_Filter_Chebyshev;
  SETERRQ(PetscObjectComm((PetscObject)st),PETSC_ERR_SUP,"Not implemented yet");
  PetscFunctionReturn(PETSC_SUCCESS);
}
