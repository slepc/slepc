/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Modification of the *temp* implementation of the BLOPEX multivector in order
   to wrap created PETSc vectors as multivectors.
*/

#include <slepc/private/bvimpl.h>
#include <stdlib.h>
#include <interpreter.h>
#include <temp_multivector.h>
#include "blopex.h"

static void *mv_TempMultiVectorCreateFromBV(void *ii_,BlopexInt n,void *sample)
{
  BV                      bv = (BV)sample;
  Vec                     v;
  PetscInt                i,l,k,nc,useconstr=PETSC_FALSE,flg;
  mv_TempMultiVector      *x;
  mv_InterfaceInterpreter *ii = (mv_InterfaceInterpreter*)ii_;

  x = (mv_TempMultiVector*)malloc(sizeof(mv_TempMultiVector));
  if (!x) SETERRABORT(PETSC_COMM_SELF,PETSC_ERR_PLIB,"Allocation for x failed");

  x->interpreter = ii;
  x->numVectors  = n;

  x->vector = (void**)calloc(n,sizeof(void*));
  if (!x->vector) SETERRABORT(PETSC_COMM_SELF,PETSC_ERR_PLIB,"Allocation for x->vector failed");

  x->ownsVectors = 1;
  x->mask = NULL;
  x->ownsMask = 0;

  PetscCallAbort(PETSC_COMM_SELF,BVGetActiveColumns(bv,&l,&k));
  PetscCallAbort(PETSC_COMM_SELF,PetscObjectComposedDataGetInt((PetscObject)bv,slepc_blopex_useconstr,useconstr,flg));
  if (!l && flg && useconstr) {
    PetscCallAbort(PETSC_COMM_SELF,BVGetNumConstraints(bv,&nc));
    l = -nc;
  }
  if (n != k-l) SETERRABORT(PETSC_COMM_SELF,PETSC_ERR_PLIB,"BV active columns plus constraints do not match argument n");
  for (i=0;i<n;i++) {
    PetscCallAbort(PETSC_COMM_SELF,BVGetColumn(bv,l+i,&v));
    PetscCallAbort(PETSC_COMM_SELF,PetscObjectReference((PetscObject)v));
    x->vector[i] = (void*)v;
    PetscCallAbort(PETSC_COMM_SELF,BVRestoreColumn(bv,l+i,&v));
  }
  return x;
}

/*
    Create an InterfaceInterpreter using the PETSc implementation
    but overloading CreateMultiVector that doesn't create any
    new vector.
*/
int SLEPCSetupInterpreter(mv_InterfaceInterpreter *i)
{
  PETSCSetupInterpreter(i);
  i->CreateMultiVector = mv_TempMultiVectorCreateFromBV;

  return 0;
}
