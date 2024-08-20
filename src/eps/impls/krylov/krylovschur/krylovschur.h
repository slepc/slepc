/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Private header for Krylov-Schur
*/

#pragma once

SLEPC_INTERN PetscErrorCode EPSSolve_KrylovSchur_Default(EPS);
SLEPC_INTERN PetscErrorCode EPSSolve_KrylovSchur_TwoSided(EPS);
SLEPC_INTERN PetscErrorCode EPSSolve_KrylovSchur_Slice(EPS);
SLEPC_INTERN PetscErrorCode EPSSetUp_KrylovSchur_Slice(EPS);
SLEPC_INTERN PetscErrorCode EPSReset_KrylovSchur_Slice(EPS);
SLEPC_INTERN PetscErrorCode EPSDestroy_KrylovSchur_Slice(EPS);
SLEPC_INTERN PetscErrorCode EPSSolve_KrylovSchur_Indefinite(EPS);
SLEPC_INTERN PetscErrorCode EPSSetUp_KrylovSchur_BSE(EPS);
SLEPC_INTERN PetscErrorCode EPSSolve_KrylovSchur_BSE_Shao(EPS);
SLEPC_INTERN PetscErrorCode EPSSolve_KrylovSchur_BSE_Gruning(EPS);
SLEPC_INTERN PetscErrorCode EPSSolve_KrylovSchur_BSE_ProjectedBSE(EPS);
SLEPC_INTERN PetscErrorCode EPSSetUp_KrylovSchur_Hamilt(EPS);
SLEPC_INTERN PetscErrorCode EPSSolve_KrylovSchur_Hamilt(EPS);
SLEPC_INTERN PetscErrorCode EPSGetArbitraryValues(EPS,PetscScalar*,PetscScalar*);
SLEPC_INTERN PetscErrorCode EPSKrylovSchurGetChildEPS(EPS,EPS*);

/* Structure characterizing a shift in spectrum slicing */
typedef struct _n_shift *EPS_shift;
struct _n_shift {
  PetscReal     value;
  PetscInt      inertia;
  PetscBool     comp[2];      /* Shows completion of subintervals (left and right) */
  EPS_shift     neighb[2];    /* Adjacent shifts */
  PetscInt      index;        /* Index in eig where found values are stored */
  PetscInt      neigs;        /* Number of values found */
  PetscReal     ext[2];       /* Limits for accepted values */
  PetscInt      nsch[2];      /* Number of missing values for each subinterval */
  PetscInt      nconv[2];     /* Converged on each side (accepted or not) */
  PetscBool     rep;          /* Present shift is reused */
};

/* Structure for storing the state of spectrum slicing */
struct _n_SR {
  PetscReal     int0,int1;    /* Extremes of the interval */
  PetscInt      dir;          /* Determines the order of values in eig (+1 incr, -1 decr) */
  PetscBool     hasEnd;       /* Tells whether the interval has an end */
  PetscInt      inertia0,inertia1;
  PetscScalar   *back;
  PetscInt      numEigs;      /* Number of eigenvalues in the interval */
  PetscInt      indexEig;
  EPS_shift     sPres;        /* Present shift */
  EPS_shift     *pending;     /* Pending shifts array */
  PetscInt      nPend;        /* Number of pending shifts */
  PetscInt      maxPend;      /* Size of "pending" array */
  PetscInt      *idxDef;      /* For deflation */
  PetscInt      nMAXCompl;
  PetscInt      iterCompl;
  PetscInt      itsKs;        /* Krylovschur restarts */
  PetscInt      nleap;
  EPS_shift     s0;           /* Initial shift */
  PetscScalar   *S;           /* Matrix for projected problem */
  PetscInt      nS;
  EPS_shift     sPrev;
  PetscInt      nv;           /* position of restart vector */
  BV            V;            /* working basis (for subsolve) */
  BV            Vnext;        /* temporary working basis during change of shift */
  PetscScalar   *eigr,*eigi;  /* eigenvalues (for subsolve) */
  PetscReal     *errest;      /* error estimates (for subsolve) */
  PetscInt      *perm;        /* permutation (for subsolve) */
};
typedef struct _n_SR *EPS_SR;

typedef struct {
  PetscReal        keep;               /* restart parameter */
  PetscBool        lock;               /* locking/non-locking variant */
  /* the following are used only in spectrum slicing */
  EPS_SR           sr;                 /* spectrum slicing context */
  PetscInt         nev;                /* number of eigenvalues to compute */
  PetscInt         ncv;                /* number of basis vectors */
  PetscInt         mpd;                /* maximum dimension of projected problem */
  PetscInt         npart;              /* number of partitions of subcommunicator */
  PetscBool        detect;             /* check for zeros during factorizations */
  PetscReal        *subintervals;      /* partition of global interval */
  PetscBool        subintset;          /* subintervals set by user */
  PetscMPIInt      *nconv_loc;         /* converged eigenpairs for each subinterval */
  EPS              eps;                /* additional eps for slice runs */
  PetscBool        global;             /* flag distinguishing global from local eps */
  PetscReal        *shifts;            /* array containing global shifts */
  PetscInt         *inertias;          /* array containing global inertias */
  PetscInt         nshifts;            /* elements in the arrays of shifts and inertias */
  PetscSubcomm     subc;               /* context for subcommunicators */
  MPI_Comm         commrank;           /* group processes with same rank in subcommunicators */
  PetscBool        commset;            /* flag indicating that commrank was created */
  PetscObjectState Astate,Bstate;      /* state of subcommunicator matrices */
  PetscObjectId    Aid,Bid;            /* Id of subcommunicator matrices */
  IS               isrow,iscol;        /* index sets used in update of subcomm mats */
  Mat              *submata,*submatb;  /* seq matrices used in update of subcomm mats */
  /* the following are used only in filter */
  PetscBool        estimatedrange;     /* the filter range was not set by the user */
  /* the following are used only for BSE problem type */
  EPSKrylovSchurBSEType bse;           /* the BSE method */
} EPS_KRYLOVSCHUR;
