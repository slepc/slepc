/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#pragma once

typedef struct {
  PetscErrorCode  (*computeoperator)(ST,Mat*);
  PetscErrorCode  (*getthreshold)(ST,PetscReal*);
  PetscErrorCode  (*destroy)(ST);

  STFilterType    type;            /* the filter method */
  PetscReal       inta,intb;       /* bounds of the interval of desired eigenvalues */
  PetscReal       left,right;      /* approximate left and right bounds of the interval containing all eigenvalues */
  PetscInt        polyDegree;      /* degree of the polynomial filter */
  STFilterDamping damping;         /* the type of damping */
  Mat             T;               /* the matrix used to build the filter */
  PetscBool       filtch;          /* filter parameters have changed since last setup */
  Mat             *W;              /* work matrices for the matrix-matrix application of the filter */
  PetscInt        nW;              /* number of W matrices */
  void            *data;           /* placeholder for method-specific stuff */
} ST_FILTER;

/* FILTLAN */

/* IntervalOptions structure used by GetIntervals */
struct _n_FILTLAN_IOP {
  PetscReal   intervalWeights[5];  /* weight of the subintervals (5 in mid-pass, 3 in high-pass) */
  PetscReal   transIntervalRatio;  /* the (relative) length of transition interval */
  PetscBool   reverseInterval;     /* whether to reverse the interval or not; effective only for mid-pass filters */
  PetscReal   initialPlateau;      /* initial plateau relative to the length of interval */
  PetscReal   plateauShrinkRate;   /* the rate at which the plateau shrinks at each iteration */
  PetscReal   initialShiftStep;    /* initial shift step relative to the length of interval */
  PetscReal   shiftStepExpanRate;  /* the rate at which the shift step expands */
  PetscInt    maxInnerIter;        /* maximum number of inner iterations to determine the (transition) intervals */
  PetscReal   yLimitTol;           /* tolerance allowed for |p(inta)-p(intb)| in a mid-pass filter p(x) */
  PetscInt    maxOuterIter;        /* maximum number of outer iterations to determine the (transition) intervals */
  PetscInt    numGridPoints;       /* number of grid points, used to measure the maximum p(z) for z not in the interval*/
  PetscReal   yBottomLine;         /* p(x) should be greater than this value for x in the interval */
  PetscReal   yRippleLimit;        /* limit of height of ripples relative to the bottom of polynomial values */
};
typedef struct _n_FILTLAN_IOP *FILTLAN_IOP;

/* PolynomialFilterInfo structure filled by GetIntervals */
struct _n_FILTLAN_PFI {
  PetscInt    filterType;          /* 1 = high-pass filter (or low-pass filter with conversion); 2 = mid-pass filter */
  PetscInt    filterOK;            /* 0 = no acceptable found; 1 = OK filter found; 2 = optimal filter found */
  PetscReal   filterQualityIndex;  /* between 0.0 and 1.0; the higher the better quality of the filter */
  PetscInt    numIter;             /* number of iterations to get the (transition) intervals */
  PetscInt    totalNumIter;        /* total number of iterations performed */
  PetscReal   yLimit;              /* lowest polynomial value P(z) for z in the interval [a1,b1] of desired eigenvalues */
  PetscReal   ySummit;             /* height of (highest, if more than one) summit in interval [a1,b1] of desired evals */
  PetscInt    numLeftSteps;        /* number of steps moving leftward */
  PetscReal   yLeftSummit;         /* height of highest summit in the left-hand side of the interval of desired evals */
  PetscReal   yLeftBottom;         /* height of lowest bottom in the left-hand side of the interval of desired evals */
  PetscReal   yLimitGap;           /* |p(a1)-p(b1)|, where [a1,b1] is the interval of desired eigenvalues */
  PetscInt    numRightSteps;       /* number of steps moving rightward */
  PetscReal   yRightSummit;        /* height of highest summit in the right-hand side of the interval of desired evals */
  PetscReal   yRightBottom;        /* height of lowest bottom in the right-hand side of the interval of desired evals */
};
typedef struct _n_FILTLAN_PFI *FILTLAN_PFI;

struct _n_FILTLAN_CTX {
  PetscInt      baseDegree;        /* left and right degrees of the base filter for each interval */
  FILTLAN_IOP   opts;              /* interval options */
  FILTLAN_PFI   filterInfo;        /* polynomial filter info */
  PetscReal     frame[4];          /* outer and inner intervals:
                                      [frame[0],frame[3]] (tightly) contains all eigenvalues
                                      [frame[1],frame[2]] is the interval in which the eigenvalues are sought */
  PetscReal     intervals[6];      /* the range of eigenvalues is partitioned into intervals which determine
                                      the base filter approximated by a polynomial filter;
                                      the j-th interval is [intervals(j),intervals(j+1)) */
  PetscReal     intervals2[6];     /* modified intervals */
  PetscReal     *baseFilter;       /* coefficients of the base filter (a piecewise polynomial) */
};
typedef struct _n_FILTLAN_CTX *FILTLAN_CTX;

SLEPC_INTERN PetscErrorCode STCreate_Filter_FILTLAN(ST);

/* Chebyshev series */

struct _n_CHEBYSHEV_CTX {
  PetscReal   *damping_coeffs;
  PetscReal   *coeffs;
};
typedef struct _n_CHEBYSHEV_CTX *CHEBYSHEV_CTX;

SLEPC_INTERN PetscErrorCode STCreate_Filter_Chebyshev(ST);
