/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   User interface for the basis vectors object in SLEPc
*/

#pragma once

#include <slepcsys.h>

/* SUBMANSEC = BV */

SLEPC_EXTERN PetscErrorCode BVInitializePackage(void);
SLEPC_EXTERN PetscErrorCode BVFinalizePackage(void);

/*S
    BV - Basis vectors, SLEPc object representing a collection of vectors
    that typically constitute a basis of a subspace.

    Level: beginner

.seealso:  BVCreate()
S*/
typedef struct _p_BV* BV;

/*J
    BVType - String with the name of the type of BV. Each type differs in
    the way data is stored internally.

    Level: beginner

.seealso: BVSetType(), BV
J*/
typedef const char *BVType;
#define BVMAT        "mat"
#define BVSVEC       "svec"
#define BVVECS       "vecs"
#define BVCONTIGUOUS "contiguous"
#define BVTENSOR     "tensor"

/* Logging support */
SLEPC_EXTERN PetscClassId BV_CLASSID;

/*E
    BVOrthogType - Determines the method used in the orthogonalization
    of vectors

    Level: advanced

.seealso: BVSetOrthogonalization(), BVGetOrthogonalization(), BVOrthogonalizeColumn(), BVOrthogRefineType
E*/
typedef enum { BV_ORTHOG_CGS,
               BV_ORTHOG_MGS } BVOrthogType;
SLEPC_EXTERN const char *BVOrthogTypes[];

/*E
    BVOrthogRefineType - Determines what type of refinement to use
    during orthogonalization of vectors

    Level: advanced

.seealso: BVSetOrthogonalization(), BVGetOrthogonalization(), BVOrthogonalizeColumn()
E*/
typedef enum { BV_ORTHOG_REFINE_IFNEEDED,
               BV_ORTHOG_REFINE_NEVER,
               BV_ORTHOG_REFINE_ALWAYS } BVOrthogRefineType;
SLEPC_EXTERN const char *BVOrthogRefineTypes[];

/*E
    BVOrthogBlockType - Determines the method used in block
    orthogonalization (simultaneous orthogonalization of a set of vectors)

    Level: advanced

.seealso: BVSetOrthogonalization(), BVGetOrthogonalization(), BVOrthogonalize()
E*/
typedef enum { BV_ORTHOG_BLOCK_GS,
               BV_ORTHOG_BLOCK_CHOL,
               BV_ORTHOG_BLOCK_TSQR,
               BV_ORTHOG_BLOCK_TSQRCHOL,
               BV_ORTHOG_BLOCK_SVQB     } BVOrthogBlockType;
SLEPC_EXTERN const char *BVOrthogBlockTypes[];

/*E
   BVMatMultType - Different ways of performing the BVMatMult() operation

   Notes:
   Allowed values are
+  BV_MATMULT_VECS - perform a matrix-vector multiply per each column
.  BV_MATMULT_MAT - carry out a Mat-Mat product with a dense matrix
-  BV_MATMULT_MAT_SAVE - this case is deprecated

   The default is BV_MATMULT_MAT except in the case of BVVECS.

   Level: advanced

.seealso: BVSetMatMultMethod(), BVMatMult()
E*/
typedef enum { BV_MATMULT_VECS,
               BV_MATMULT_MAT,
               BV_MATMULT_MAT_SAVE } BVMatMultType;
SLEPC_EXTERN const char *BVMatMultTypes[];

SLEPC_EXTERN PetscErrorCode BVCreate(MPI_Comm,BV*);
SLEPC_EXTERN PetscErrorCode BVDestroy(BV*);
SLEPC_EXTERN PetscErrorCode BVSetType(BV,BVType);
SLEPC_EXTERN PetscErrorCode BVGetType(BV,BVType*);
SLEPC_EXTERN PetscErrorCode BVSetSizes(BV,PetscInt,PetscInt,PetscInt);
SLEPC_EXTERN PetscErrorCode BVSetSizesFromVec(BV,Vec,PetscInt);
SLEPC_EXTERN PetscErrorCode BVSetLeadingDimension(BV,PetscInt);
SLEPC_EXTERN PetscErrorCode BVGetLeadingDimension(BV,PetscInt*);
SLEPC_EXTERN PetscErrorCode BVGetSizes(BV,PetscInt*,PetscInt*,PetscInt*);
SLEPC_EXTERN PetscErrorCode BVResize(BV,PetscInt,PetscBool);
SLEPC_EXTERN PetscErrorCode BVSetFromOptions(BV);
SLEPC_EXTERN PetscErrorCode BVView(BV,PetscViewer);
SLEPC_EXTERN PetscErrorCode BVViewFromOptions(BV,PetscObject,const char[]);

SLEPC_EXTERN PetscErrorCode BVGetColumn(BV,PetscInt,Vec*);
SLEPC_EXTERN PetscErrorCode BVRestoreColumn(BV,PetscInt,Vec*);
SLEPC_EXTERN PetscErrorCode BVGetSplit(BV,BV*,BV*);
SLEPC_EXTERN PetscErrorCode BVRestoreSplit(BV,BV*,BV*);
SLEPC_EXTERN PetscErrorCode BVGetSplitRows(BV,IS,IS,BV*,BV*);
SLEPC_EXTERN PetscErrorCode BVRestoreSplitRows(BV,IS,IS,BV*,BV*);
SLEPC_EXTERN PetscErrorCode BVGetArray(BV,PetscScalar*[]);
SLEPC_EXTERN PetscErrorCode BVRestoreArray(BV,PetscScalar*[]);
SLEPC_EXTERN PetscErrorCode BVGetArrayRead(BV,const PetscScalar*[]);
SLEPC_EXTERN PetscErrorCode BVRestoreArrayRead(BV,const PetscScalar*[]);
SLEPC_EXTERN PetscErrorCode BVCreateVec(BV,Vec*);
SLEPC_EXTERN PetscErrorCode BVCreateVecEmpty(BV,Vec*);
SLEPC_EXTERN PetscErrorCode BVSetVecType(BV,VecType);
SLEPC_EXTERN PetscErrorCode BVGetVecType(BV,VecType*);
SLEPC_EXTERN PetscErrorCode BVSetActiveColumns(BV,PetscInt,PetscInt);
SLEPC_EXTERN PetscErrorCode BVGetActiveColumns(BV,PetscInt*,PetscInt*);
SLEPC_EXTERN PetscErrorCode BVInsertVec(BV,PetscInt,Vec);
SLEPC_EXTERN PetscErrorCode BVInsertVecs(BV,PetscInt,PetscInt*,Vec*,PetscBool);
SLEPC_EXTERN PetscErrorCode BVInsertConstraints(BV,PetscInt*,Vec*);
SLEPC_EXTERN PetscErrorCode BVSetNumConstraints(BV,PetscInt);
SLEPC_EXTERN PetscErrorCode BVGetNumConstraints(BV,PetscInt*);
SLEPC_EXTERN PetscErrorCode BVSetDefiniteTolerance(BV,PetscReal);
SLEPC_EXTERN PetscErrorCode BVGetDefiniteTolerance(BV,PetscReal*);
SLEPC_EXTERN PetscErrorCode BVDuplicate(BV,BV*);
SLEPC_EXTERN PetscErrorCode BVDuplicateResize(BV,PetscInt,BV*);
SLEPC_EXTERN PetscErrorCode BVCopy(BV,BV);
SLEPC_EXTERN PetscErrorCode BVCopyVec(BV,PetscInt,Vec);
SLEPC_EXTERN PetscErrorCode BVCopyColumn(BV,PetscInt,PetscInt);
SLEPC_EXTERN PetscErrorCode BVSetMatrix(BV,Mat,PetscBool);
SLEPC_EXTERN PetscErrorCode BVGetMatrix(BV,Mat*,PetscBool*);
SLEPC_EXTERN PetscErrorCode BVApplyMatrix(BV,Vec,Vec);
SLEPC_EXTERN PetscErrorCode BVApplyMatrixBV(BV,BV);
SLEPC_EXTERN PetscErrorCode BVGetCachedBV(BV,BV*);
SLEPC_EXTERN PetscErrorCode BVSetSignature(BV,Vec);
SLEPC_EXTERN PetscErrorCode BVGetSignature(BV,Vec);
SLEPC_EXTERN PetscErrorCode BVSetBufferVec(BV,Vec);
SLEPC_EXTERN PetscErrorCode BVGetBufferVec(BV,Vec*);

SLEPC_EXTERN PetscErrorCode BVMult(BV,PetscScalar,PetscScalar,BV,Mat);
SLEPC_EXTERN PetscErrorCode BVMultVec(BV,PetscScalar,PetscScalar,Vec,PetscScalar*);
SLEPC_EXTERN PetscErrorCode BVMultColumn(BV,PetscScalar,PetscScalar,PetscInt,PetscScalar*);
SLEPC_EXTERN PetscErrorCode BVMultInPlace(BV,Mat,PetscInt,PetscInt);
SLEPC_EXTERN PetscErrorCode BVMultInPlaceHermitianTranspose(BV,Mat,PetscInt,PetscInt);
PETSC_DEPRECATED_FUNCTION(3, 16, 0, "BVMultInPlaceHermitianTranspose()", ) static inline PetscErrorCode BVMultInPlaceTranspose(BV bv,Mat A,PetscInt s,PetscInt e) {return BVMultInPlaceHermitianTranspose(bv,A,s,e);}
SLEPC_EXTERN PetscErrorCode BVMatMult(BV,Mat,BV);
SLEPC_EXTERN PetscErrorCode BVMatMultTranspose(BV,Mat,BV);
SLEPC_EXTERN PetscErrorCode BVMatMultHermitianTranspose(BV,Mat,BV);
SLEPC_EXTERN PetscErrorCode BVMatMultColumn(BV,Mat,PetscInt);
SLEPC_EXTERN PetscErrorCode BVMatMultTransposeColumn(BV,Mat,PetscInt);
SLEPC_EXTERN PetscErrorCode BVMatMultHermitianTransposeColumn(BV,Mat,PetscInt);
SLEPC_EXTERN PetscErrorCode BVMatProject(BV,Mat,BV,Mat);
SLEPC_EXTERN PetscErrorCode BVMatArnoldi(BV,Mat,Mat,PetscInt,PetscInt*,PetscReal*,PetscBool*);
SLEPC_EXTERN PetscErrorCode BVMatLanczos(BV,Mat,Mat,PetscInt,PetscInt*,PetscReal*,PetscBool*);

SLEPC_EXTERN PetscErrorCode BVDot(BV,BV,Mat);
SLEPC_EXTERN PetscErrorCode BVDotVec(BV,Vec,PetscScalar*);
SLEPC_EXTERN PetscErrorCode BVDotVecBegin(BV,Vec,PetscScalar*);
SLEPC_EXTERN PetscErrorCode BVDotVecEnd(BV,Vec,PetscScalar*);
SLEPC_EXTERN PetscErrorCode BVDotColumn(BV,PetscInt,PetscScalar*);
SLEPC_EXTERN PetscErrorCode BVDotColumnBegin(BV,PetscInt,PetscScalar*);
SLEPC_EXTERN PetscErrorCode BVDotColumnEnd(BV,PetscInt,PetscScalar*);
SLEPC_EXTERN PetscErrorCode BVScale(BV,PetscScalar);
SLEPC_EXTERN PetscErrorCode BVScaleColumn(BV,PetscInt,PetscScalar);
SLEPC_EXTERN PetscErrorCode BVNorm(BV,NormType,PetscReal*);
SLEPC_EXTERN PetscErrorCode BVNormVec(BV,Vec,NormType,PetscReal*);
SLEPC_EXTERN PetscErrorCode BVNormVecBegin(BV,Vec,NormType,PetscReal*);
SLEPC_EXTERN PetscErrorCode BVNormVecEnd(BV,Vec,NormType,PetscReal*);
SLEPC_EXTERN PetscErrorCode BVNormColumn(BV,PetscInt,NormType,PetscReal*);
SLEPC_EXTERN PetscErrorCode BVNormColumnBegin(BV,PetscInt,NormType,PetscReal*);
SLEPC_EXTERN PetscErrorCode BVNormColumnEnd(BV,PetscInt,NormType,PetscReal*);
SLEPC_EXTERN PetscErrorCode BVNormalize(BV,PetscScalar*);
SLEPC_EXTERN PetscErrorCode BVSetRandom(BV);
SLEPC_EXTERN PetscErrorCode BVSetRandomNormal(BV);
SLEPC_EXTERN PetscErrorCode BVSetRandomSign(BV);
SLEPC_EXTERN PetscErrorCode BVSetRandomColumn(BV,PetscInt);
SLEPC_EXTERN PetscErrorCode BVSetRandomCond(BV,PetscReal);
SLEPC_EXTERN PetscErrorCode BVSetRandomContext(BV,PetscRandom);
SLEPC_EXTERN PetscErrorCode BVGetRandomContext(BV,PetscRandom*);

SLEPC_EXTERN PetscErrorCode BVSetOrthogonalization(BV,BVOrthogType,BVOrthogRefineType,PetscReal,BVOrthogBlockType);
SLEPC_EXTERN PetscErrorCode BVGetOrthogonalization(BV,BVOrthogType*,BVOrthogRefineType*,PetscReal*,BVOrthogBlockType*);
SLEPC_EXTERN PetscErrorCode BVOrthogonalize(BV,Mat);
SLEPC_EXTERN PetscErrorCode BVOrthogonalizeVec(BV,Vec,PetscScalar*,PetscReal*,PetscBool*);
SLEPC_EXTERN PetscErrorCode BVOrthogonalizeColumn(BV,PetscInt,PetscScalar*,PetscReal*,PetscBool*);
SLEPC_EXTERN PetscErrorCode BVOrthonormalizeColumn(BV,PetscInt,PetscBool,PetscReal*,PetscBool*);
SLEPC_EXTERN PetscErrorCode BVOrthogonalizeSomeColumn(BV,PetscInt,PetscBool*,PetscScalar*,PetscReal*,PetscBool*);
SLEPC_EXTERN PetscErrorCode BVBiorthogonalizeColumn(BV,BV,PetscInt);
SLEPC_EXTERN PetscErrorCode BVBiorthonormalizeColumn(BV,BV,PetscInt,PetscReal*);
SLEPC_EXTERN PetscErrorCode BVSetMatMultMethod(BV,BVMatMultType);
SLEPC_EXTERN PetscErrorCode BVGetMatMultMethod(BV,BVMatMultType*);

SLEPC_EXTERN PetscErrorCode BVCreateFromMat(Mat,BV*);
SLEPC_EXTERN PetscErrorCode BVCreateMat(BV,Mat*);
SLEPC_EXTERN PetscErrorCode BVGetMat(BV,Mat*);
SLEPC_EXTERN PetscErrorCode BVRestoreMat(BV,Mat*);

SLEPC_EXTERN PetscErrorCode BVScatter(BV,BV,VecScatter,Vec);
SLEPC_EXTERN PetscErrorCode BVSumQuadrature(BV,BV,PetscInt,PetscInt,PetscInt,PetscScalar*,PetscScalar*,VecScatter,PetscSubcomm,PetscInt,PetscBool);
SLEPC_EXTERN PetscErrorCode BVDotQuadrature(BV,BV,PetscScalar*,PetscInt,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscSubcomm,PetscInt,PetscBool);
SLEPC_EXTERN PetscErrorCode BVTraceQuadrature(BV,BV,PetscInt,PetscInt,PetscScalar*,VecScatter,PetscSubcomm,PetscInt,PetscBool,PetscReal*);

/*E
   BVSVDMethod - Different methods for computing the SVD of a BV

   Notes:
   Allowed values are
+  BV_SVD_METHOD_REFINE - based on the SVD of the cross product matrix S'*S, with refinement
.  BV_SVD_METHOD_QR     - based on the SVD of the triangular factor of qr(S)
-  BV_SVD_METHOD_QR_CAA - variant of QR intended for use in cammunication-avoiding Arnoldi

   Level: developer

.seealso: BVSVDAndRank()
E*/
typedef enum { BV_SVD_METHOD_REFINE,
               BV_SVD_METHOD_QR,
               BV_SVD_METHOD_QR_CAA } BVSVDMethod;
SLEPC_EXTERN const char *BVSVDMethods[];

SLEPC_EXTERN PetscErrorCode BVSVDAndRank(BV,PetscInt,PetscInt,PetscReal,BVSVDMethod,PetscScalar*,PetscReal*,PetscInt*);
SLEPC_EXTERN PetscErrorCode BVCISSResizeBases(BV,BV,BV,PetscInt,PetscInt,PetscInt,PetscInt);

SLEPC_EXTERN PetscErrorCode BVCreateTensor(BV,PetscInt,BV*);
SLEPC_EXTERN PetscErrorCode BVTensorBuildFirstColumn(BV,PetscInt);
SLEPC_EXTERN PetscErrorCode BVTensorCompress(BV,PetscInt);
SLEPC_EXTERN PetscErrorCode BVTensorGetDegree(BV,PetscInt*);
SLEPC_EXTERN PetscErrorCode BVTensorGetFactors(BV,BV*,Mat*);
SLEPC_EXTERN PetscErrorCode BVTensorRestoreFactors(BV,BV*,Mat*);

SLEPC_EXTERN PetscErrorCode BVSetOptionsPrefix(BV,const char*);
SLEPC_EXTERN PetscErrorCode BVAppendOptionsPrefix(BV,const char*);
SLEPC_EXTERN PetscErrorCode BVGetOptionsPrefix(BV,const char*[]);

SLEPC_EXTERN PetscFunctionList BVList;
SLEPC_EXTERN PetscErrorCode BVRegister(const char[],PetscErrorCode(*)(BV));
