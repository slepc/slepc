/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#pragma once

#include <slepcsys.h>
#include <petsc/private/petscimpl.h>

/* SUBMANSEC = sys */

SLEPC_INTERN PetscBool SlepcBeganPetsc;

/* SlepcSwap - swap two variables a,b of the same type using a temporary variable t */
#define SlepcSwap(a,b,t) do {t=a;a=b;b=t;} while (0)

/*MC
    SlepcHeaderCreate - Creates a SLEPc object

    Input Parameters:
+   classid - the classid associated with this object
.   class_name - string name of class; should be static
.   descr - string containing short description; should be static
.   mansec - string indicating section in manual pages; should be static
.   comm - the MPI Communicator
.   destroy - the destroy routine for this object
-   view - the view routine for this object

    Output Parameter:
.   h - the newly created object

    Note:
    This is equivalent to PetscHeaderCreate but makes sure that SlepcInitialize
    has been called.

    Level: developer
M*/
#define SlepcHeaderCreate(h,classid,class_name,descr,mansec,comm,destroy,view) \
    ((PetscErrorCode)((!SlepcInitializeCalled && \
                       PetscError(comm,__LINE__,PETSC_FUNCTION_NAME,__FILE__,PETSC_ERR_ORDER,PETSC_ERROR_INITIAL, \
                                  "Must call SlepcInitialize instead of PetscInitialize to use SLEPc classes")) || \
                      PetscHeaderCreate(h,classid,class_name,descr,mansec,comm,destroy,view)))

/* context for monitors of type XXXMonitorConverged */
struct _n_SlepcConvMon {
  void     *ctx;
  PetscInt oldnconv;  /* previous value of nconv */
};

/* context for structured eigenproblem matrices created via MatCreateXXX */
struct _n_SlepcMatStruct {
  PetscInt    cookie;    /* identify which structured matrix */
  PetscScalar s;         /* in BSE sign of the bottom part of the vector */
};
typedef struct _n_SlepcMatStruct* SlepcMatStruct;

#define SLEPC_MAT_STRUCT_BSE 88101

/*
  SlepcCheckMatStruct - Check that a given Mat is a structured matrix of the wanted type.

  Returns true/false in flg if it is given, otherwise yields an error if the check fails.
  If cookie==0 it will check for any type.
*/
static inline PetscErrorCode SlepcCheckMatStruct(Mat A,PetscInt cookie,PetscBool *flg)
{
  PetscContainer container;
  SlepcMatStruct mctx;

  PetscFunctionBegin;
  if (flg) *flg = PETSC_FALSE;
  PetscCall(PetscObjectQuery((PetscObject)A,"SlepcMatStruct",(PetscObject*)&container));
  if (flg && !container) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCheck(container,PetscObjectComm((PetscObject)A),PETSC_ERR_ARG_WRONG,"The Mat is not a structured matrix");
  if (cookie) {
    PetscCall(PetscContainerGetPointer(container,(void**)&mctx));
    if (flg && (!mctx || mctx->cookie!=cookie)) PetscFunctionReturn(PETSC_SUCCESS);
    PetscCheck(mctx && mctx->cookie==cookie,PetscObjectComm((PetscObject)A),PETSC_ERR_ARG_WRONG,"The type of structured matrix is different from the expected one");
  }
  if (flg) *flg = PETSC_TRUE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  SlepcPrintEigenvalueASCII - Print an eigenvalue on an ASCII viewer.
*/
static inline PetscErrorCode SlepcPrintEigenvalueASCII(PetscViewer viewer,PetscScalar eigr,PetscScalar eigi)
{
  PetscReal      re,im;

  PetscFunctionBegin;
#if defined(PETSC_USE_COMPLEX)
  re = PetscRealPart(eigr);
  im = PetscImaginaryPart(eigr);
  (void)eigi;
#else
  re = eigr;
  im = eigi;
#endif
  /* print zero instead of tiny value */
  if (PetscAbs(im) && PetscAbs(re)/PetscAbs(im)<PETSC_SMALL) re = 0.0;
  if (PetscAbs(re) && PetscAbs(im)/PetscAbs(re)<PETSC_SMALL) im = 0.0;
  /* print as real if imaginary part is zero */
  if (im!=(PetscReal)0.0) PetscCall(PetscViewerASCIIPrintf(viewer,"%.5f%+.5fi",(double)re,(double)im));
  else PetscCall(PetscViewerASCIIPrintf(viewer,"%.5f",(double)re));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  SlepcViewEigenvector - Outputs an eigenvector xr,xi to a viewer.
  In complex scalars only xr is written.
  The name of xr,xi is set before writing, based on the label, the index, and the name of obj.
*/
static inline PetscErrorCode SlepcViewEigenvector(PetscViewer viewer,Vec xr,Vec xi,const char *label,PetscInt index,PetscObject obj)
{
  size_t         count;
  char           vname[30];
  const char     *pname;

  PetscFunctionBegin;
  PetscCall(PetscObjectGetName(obj,&pname));
  PetscCall(PetscSNPrintfCount(vname,sizeof(vname),"%s%s",&count,label,PetscDefined(USE_COMPLEX)?"":"r"));
  count--;
  PetscCall(PetscSNPrintf(vname+count,sizeof(vname)-count,"%" PetscInt_FMT "_%s",index,pname));
  PetscCall(PetscObjectSetName((PetscObject)xr,vname));
  PetscCall(VecView(xr,viewer));
#if !defined(PETSC_USE_COMPLEX)
  vname[count-1] = 'i';
  PetscCall(PetscObjectSetName((PetscObject)xi,vname));
  PetscCall(VecView(xi,viewer));
#else
  (void)xi;
#endif
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* Macros for strings with different value in real and complex */
#if defined(PETSC_USE_COMPLEX)
#define SLEPC_STRING_HERMITIAN "hermitian"
#else
#define SLEPC_STRING_HERMITIAN "symmetric"
#endif

/* Private functions that are shared by several classes */
SLEPC_SINGLE_LIBRARY_INTERN PetscErrorCode SlepcBasisReference_Private(PetscInt,Vec*,PetscInt*,Vec**);
SLEPC_SINGLE_LIBRARY_INTERN PetscErrorCode SlepcBasisDestroy_Private(PetscInt*,Vec**);
SLEPC_SINGLE_LIBRARY_INTERN PetscErrorCode SlepcMonitorMakeKey_Internal(const char[],PetscViewerType,PetscViewerFormat,char[]);
SLEPC_SINGLE_LIBRARY_INTERN PetscErrorCode PetscViewerAndFormatCreate_Internal(PetscViewer,PetscViewerFormat,void*,PetscViewerAndFormat**);

SLEPC_INTERN PetscErrorCode SlepcCitationsInitialize(void);
SLEPC_INTERN PetscErrorCode SlepcInitialize_DynamicLibraries(void);
SLEPC_INTERN PetscErrorCode SlepcInitialize_Packages(void);

/* Macro to check a sequential Mat (including GPU) */
#if !defined(PETSC_USE_DEBUG)
#define SlepcMatCheckSeq(h) do {(void)(h);} while (0)
#else
#if defined(PETSC_HAVE_CUDA)
#define SlepcMatCheckSeq(h) do { PetscCheckTypeNames((h),MATSEQDENSE,MATSEQDENSECUDA); } while (0)
#elif defined(PETSC_HAVE_HIP)
#define SlepcMatCheckSeq(h) do { PetscCheckTypeNames((h),MATSEQDENSE,MATSEQDENSEHIP); } while (0)
#else
#define SlepcMatCheckSeq(h) do { PetscCheckTypeName((h),MATSEQDENSE); } while (0)
#endif
#endif

/* Definitions needed to work with GPU kernels */
#if defined(PETSC_HAVE_CUPM)
#include <petscdevice_cupm.h>

#define X_AXIS 0
#define Y_AXIS 1

#define SLEPC_TILE_SIZE_X  32
#define SLEPC_BLOCK_SIZE_X 128
#define SLEPC_TILE_SIZE_Y  32
#define SLEPC_BLOCK_SIZE_Y 128

static inline PetscErrorCode SlepcKernelSetGrid1D(PetscInt rows,dim3 *dimGrid,dim3 *dimBlock,PetscInt *dimGrid_xcount)
{
  int card;
#if defined(PETSC_HAVE_CUDA)
  struct cudaDeviceProp devprop;
#elif defined(PETSC_HAVE_HIP)
  hipDeviceProp_t devprop;
#endif

  PetscFunctionBegin;
#if defined(PETSC_HAVE_CUDA)
  PetscCallCUDA(cudaGetDevice(&card));
  PetscCallCUDA(cudaGetDeviceProperties(&devprop,card));
#elif defined(PETSC_HAVE_HIP)
  PetscCallHIP(hipGetDevice(&card));
  PetscCallHIP(hipGetDeviceProperties(&devprop,card));
#endif
  *dimGrid_xcount = 1;

  /* X axis */
  dimGrid->x  = 1;
  dimBlock->x = SLEPC_BLOCK_SIZE_X;
  if (rows>SLEPC_BLOCK_SIZE_X) dimGrid->x = (rows+SLEPC_BLOCK_SIZE_X-1)/SLEPC_BLOCK_SIZE_X;
  else dimBlock->x = rows;
  if (dimGrid->x>(unsigned)devprop.maxGridSize[X_AXIS]) {
    *dimGrid_xcount = (dimGrid->x+(devprop.maxGridSize[X_AXIS]-1))/devprop.maxGridSize[X_AXIS];
    dimGrid->x = devprop.maxGridSize[X_AXIS];
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static inline PetscErrorCode SlepcKernelSetGrid2DTiles(PetscInt rows,PetscInt cols,dim3 *dimGrid,dim3 *dimBlock,PetscInt *dimGrid_xcount,PetscInt *dimGrid_ycount)
{
  int card;
#if defined(PETSC_HAVE_CUDA)
  struct cudaDeviceProp devprop;
#elif defined(PETSC_HAVE_HIP)
  hipDeviceProp_t devprop;
#endif

  PetscFunctionBegin;
#if defined(PETSC_HAVE_CUDA)
  PetscCallCUDA(cudaGetDevice(&card));
  PetscCallCUDA(cudaGetDeviceProperties(&devprop,card));
#elif defined(PETSC_HAVE_HIP)
  PetscCallHIP(hipGetDevice(&card));
  PetscCallHIP(hipGetDeviceProperties(&devprop,card));
#endif
  *dimGrid_xcount = *dimGrid_ycount = 1;

  /* X axis */
  dimGrid->x  = 1;
  dimBlock->x = SLEPC_BLOCK_SIZE_X;
  if (rows>SLEPC_BLOCK_SIZE_X*SLEPC_TILE_SIZE_X) dimGrid->x = (rows+SLEPC_BLOCK_SIZE_X*SLEPC_TILE_SIZE_X-1)/(SLEPC_BLOCK_SIZE_X*SLEPC_TILE_SIZE_X);
  else dimBlock->x = (rows+SLEPC_TILE_SIZE_X-1)/SLEPC_TILE_SIZE_X;
  if (dimGrid->x>(unsigned)devprop.maxGridSize[X_AXIS]) {
    *dimGrid_xcount = (dimGrid->x+(devprop.maxGridSize[X_AXIS]-1))/devprop.maxGridSize[X_AXIS];
    dimGrid->x = devprop.maxGridSize[X_AXIS];
  }

  /* Y axis */
  dimGrid->y  = 1;
  dimBlock->y = SLEPC_BLOCK_SIZE_Y;
  if (cols>SLEPC_BLOCK_SIZE_Y*SLEPC_TILE_SIZE_Y) dimGrid->y = (cols+SLEPC_BLOCK_SIZE_Y*SLEPC_TILE_SIZE_Y-1)/(SLEPC_BLOCK_SIZE_Y*SLEPC_TILE_SIZE_Y);
  else dimBlock->y = (cols+SLEPC_TILE_SIZE_Y-1)/SLEPC_TILE_SIZE_Y;
  if (dimGrid->y>(unsigned)devprop.maxGridSize[Y_AXIS]) {
    *dimGrid_ycount = (dimGrid->y+(devprop.maxGridSize[Y_AXIS]-1))/devprop.maxGridSize[Y_AXIS];
    dimGrid->y = devprop.maxGridSize[Y_AXIS];
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}
#undef X_AXIS
#undef Y_AXIS
#endif
