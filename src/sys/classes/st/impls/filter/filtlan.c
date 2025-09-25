/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   This file is an adaptation of several subroutines from FILTLAN, the
   Filtered Lanczos Package, authored by Haw-ren Fang and Yousef Saad.

   More information at:
   https://www-users.cs.umn.edu/~saad/software/filtlan

   References:

       [1] H. Fang and Y. Saad, "A filtered Lanczos procedure for extreme and interior
           eigenvalue problems", SIAM J. Sci. Comput. 34(4):A2220-A2246, 2012.
*/

#include <slepc/private/stimpl.h>
#include "filter.h"

static PetscErrorCode FILTLAN_FilteredConjugateResidualPolynomial(PetscReal*,PetscReal*,PetscInt,PetscReal*,PetscInt,PetscReal*,PetscInt);
static PetscReal FILTLAN_PiecewisePolynomialEvaluationInChebyshevBasis(PetscReal*,PetscInt,PetscReal*,PetscInt,PetscReal);
static PetscErrorCode FILTLAN_ExpandNewtonPolynomialInChebyshevBasis(PetscInt,PetscReal,PetscReal,PetscReal*,PetscReal*,PetscReal*,PetscReal*);

/* ////////////////////////////////////////////////////////////////////////////
   //    Newton - Hermite Polynomial Interpolation
   //////////////////////////////////////////////////////////////////////////// */

/*
   FILTLAN function NewtonPolynomial

   build P(z) by Newton's divided differences in the form
       P(z) = a(1) + a(2)*(z-x(1)) + a(3)*(z-x(1))*(z-x(2)) + ... + a(n)*(z-x(1))*...*(z-x(n-1)),
   such that P(x(i)) = y(i) for i=1,...,n, where
       x,y are input vectors of length n, and a is the output vector of length n
   if x(i)==x(j) for some i!=j, then it is assumed that the derivative of P(z) is to be zero at x(i),
       and the Hermite polynomial interpolation is applied
   in general, if there are k x(i)'s with the same value x0, then
       the j-th order derivative of P(z) is zero at z=x0 for j=1,...,k-1
*/
static inline PetscErrorCode FILTLAN_NewtonPolynomial(PetscInt n,PetscReal *x,PetscReal *y,PetscReal *sa,PetscReal *sf)
{
  PetscReal      d,*sx=x,*sy=y;
  PetscInt       j,k;

  PetscFunctionBegin;
  PetscCall(PetscArraycpy(sf,sy,n));

  /* apply Newton's finite difference method */
  sa[0] = sf[0];
  for (j=1;j<n;j++) {
    for (k=n-1;k>=j;k--) {
      d = sx[k]-sx[k-j];
      if (PetscUnlikely(d == 0.0)) sf[k] = 0.0;  /* assume that the derivative is 0.0 and apply the Hermite interpolation */
      else sf[k] = (sf[k]-sf[k-1]) / d;
    }
    sa[j] = sf[j];
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   FILTLAN function HermiteBaseFilterInChebyshevBasis

   compute a base filter P(z) which is a continuous, piecewise polynomial P(z) expanded
   in a basis of `translated' (i.e. scale-and-shift) Chebyshev polynomials in each interval

   The base filter P(z) equals P_j(z) for z in the j-th interval [intv(j), intv(j+1)), where
   P_j(z) a Hermite interpolating polynomial

   input:
   intv is a vector which defines the intervals; the j-th interval is [intv(j), intv(j+1))
   HiLowFlags determines the shape of the base filter P(z)
   Consider the j-th interval [intv(j), intv(j+1)]
   HighLowFlag[j-1]==1,  P(z)==1 for z in [intv(j), intv(j+1)]
                   ==0,  P(z)==0 for z in [intv(j), intv(j+1)]
                   ==-1, [intv(j), intv(j+1)] is a transition interval;
                         P(intv(j)) and P(intv(j+1)) are defined such that P(z) is continuous
   baseDeg is the degree of smoothness of the Hermite (piecewise polynomial) interpolation
   to be precise, the i-th derivative of P(z) is zero, i.e. d^{i}P(z)/dz^i==0, at all interval
   end points z=intv(j) for i=1,...,baseDeg

   output:
   P(z) expanded in a basis of `translated' (scale-and-shift) Chebyshev polynomials
   to be precise, for z in the j-th interval [intv(j),intv(j+1)), P(z) equals
       P_j(z) = pp(1,j)*S_0(z) + pp(2,j)*S_1(z) + ... + pp(n,j)*S_{n-1}(z),
   where S_i(z) is the `translated' Chebyshev polynomial in that interval,
       S_i((z-c)/h) = T_i(z),  c = (intv(j)+intv(j+1))) / 2,  h = (intv(j+1)-intv(j)) / 2,
   with T_i(z) the Chebyshev polynomial of the first kind,
       T_0(z) = 1, T_1(z) = z, and T_i(z) = 2*z*T_{i-1}(z) - T_{i-2}(z) for i>=2
   the return matrix is the matrix of Chebyshev coefficients pp just described

   note that the degree of P(z) in each interval is (at most) 2*baseDeg+1, with 2*baseDeg+2 coefficients
   let n be the length of intv; then there are n-1 intervals
   therefore the return matrix pp is of size (2*baseDeg+2)-by-(n-1)
*/
static PetscErrorCode FILTLAN_HermiteBaseFilterInChebyshevBasis(PetscReal *baseFilter,PetscReal *intv,PetscInt npoints,const PetscInt *HighLowFlags,PetscInt baseDeg)
{
  PetscInt       m,ii,jj;
  PetscReal      flag,flag0,flag2,aa,bb,*px,*py,*sx,*sy,*pp,*qq,*sq,*sbf,*work,*currentPoint = intv;
  const PetscInt *hilo = HighLowFlags;

  PetscFunctionBegin;
  m = 2*baseDeg+2;
  jj = npoints-1;  /* jj is initialized as the number of intervals */
  PetscCall(PetscMalloc5(m,&px,m,&py,m,&pp,m,&qq,m,&work));
  sbf = baseFilter;

  while (jj--) {  /* the main loop to compute the Chebyshev coefficients */

    flag = (PetscReal)(*hilo++);  /* get the flag of the current interval */
    if (flag == -1.0) {  /* flag == -1 means that the current interval is a transition polynomial */

      flag2 = (PetscReal)(*hilo);  /* get flag2, the flag of the next interval */
      flag0 = 1.0-flag2;       /* the flag of the previous interval is 1-flag2 */

      /* two pointers for traversing x[] and y[] */
      sx = px;
      sy = py;

      /* find the current interval [aa,bb] */
      aa = *currentPoint++;
      bb = *currentPoint;

      /* now left-hand side */
      ii = baseDeg+1;
      while (ii--) {
        *sy++ = flag0;
        *sx++ = aa;
      }

      /* now right-hand side */
      ii = baseDeg+1;
      while (ii--) {
        *sy++ = flag2;
        *sx++ = bb;
      }

      /* build a Newton polynomial (indeed, the generalized Hermite interpolating polynomial) with x[] and y[] */
      PetscCall(FILTLAN_NewtonPolynomial(m,px,py,pp,work));

      /* pp contains coefficients of the Newton polynomial P(z) in the current interval [aa,bb], where
         P(z) = pp(1) + pp(2)*(z-px(1)) + pp(3)*(z-px(1))*(z-px(2)) + ... + pp(n)*(z-px(1))*...*(z-px(n-1)) */

      /* translate the Newton coefficients to the Chebyshev coefficients */
      PetscCall(FILTLAN_ExpandNewtonPolynomialInChebyshevBasis(m,aa,bb,pp,px,qq,work));
      /* qq contains coefficients of the polynomial in [aa,bb] in the `translated' Chebyshev basis */

      /* copy the Chebyshev coefficients to baseFilter
         OCTAVE/MATLAB: B(:,j) = qq, where j = (npoints-1)-jj and B is the return matrix */
      sq = qq;
      ii = 2*baseDeg+2;
      while (ii--) *sbf++ = *sq++;

    } else {

      /* a constant polynomial P(z)=flag, where either flag==0 or flag==1
       OCTAVE/MATLAB: B(1,j) = flag, where j = (npoints-1)-jj and B is the return matrix */
      *sbf++ = flag;

      /* the other coefficients are all zero, since P(z) is a constant
       OCTAVE/MATLAB: B(1,j) = 0, where j = (npoints-1)-jj and B is the return matrix */
      ii = 2*baseDeg+1;
      while (ii--) *sbf++ = 0.0;

      /* for the next point */
      currentPoint++;
    }
  }
  PetscCall(PetscFree5(px,py,pp,qq,work));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* ////////////////////////////////////////////////////////////////////////////
   //    Base Filter
   //////////////////////////////////////////////////////////////////////////// */

/*
   FILTLAN function GetIntervals

   this routine determines the intervals (including the transition one(s)) by an iterative process

   frame is a vector consisting of 4 ordered elements:
       [frame(1),frame(4)] is the interval which (tightly) contains all eigenvalues, and
       [frame(2),frame(3)] is the interval in which the eigenvalues are sought
   baseDeg is the left-and-right degree of the base filter for each interval
   polyDeg is the (maximum possible) degree of s(z), with z*s(z) the polynomial filter
   intv is the output; the j-th interval is [intv(j),intv(j+1))
   opts is a collection of interval options

   the base filter P(z) is a piecewise polynomial from Hermite interpolation with degree baseDeg
   at each end point of intervals

   the polynomial filter Q(z) is in the form z*s(z), i.e. Q(0)==0, such that ||(1-P(z))-z*s(z)||_w is
   minimized with s(z) a polynomial of degree up to polyDeg

   the resulting polynomial filter Q(z) satisfies Q(x)>=Q(y) for x in [frame[1],frame[2]], and
   y in [frame[0],frame[3]] but not in [frame[1],frame[2]]

   the routine fills a PolynomialFilterInfo struct which gives some information of the polynomial filter
*/
static PetscErrorCode FILTLAN_GetIntervals(PetscReal *intervals,PetscReal *frame,PetscInt polyDeg,PetscInt baseDeg,FILTLAN_IOP opts,FILTLAN_PFI filterInfo)
{
  PetscReal       intv[6],x,y,z1,z2,c,c1,c2,fc,fc2,halfPlateau,leftDelta,rightDelta,gridSize;
  PetscReal       yLimit,ySummit,yLeftLimit,yRightLimit,bottom,qIndex,*baseFilter,*polyFilter;
  PetscReal       yLimitGap=0.0,yLeftSummit=0.0,yLeftBottom=0.0,yRightSummit=0.0,yRightBottom=0.0;
  PetscInt        i,ii,npoints,numIter,numLeftSteps=1,numRightSteps=1,numMoreLooked=0;
  PetscBool       leftBoundaryMet=PETSC_FALSE,rightBoundaryMet=PETSC_FALSE,stepLeft,stepRight;
  const PetscReal a=frame[0],a1=frame[1],b1=frame[2],b=frame[3];
  const PetscInt  HighLowFlags[5] = { 1, -1, 0, -1, 1 };  /* if filterType is 1, only first 3 elements will be used */
  const PetscInt  numLookMore = 2*(PetscInt)(0.5+(PetscLogReal(2.0)/PetscLogReal(opts->shiftStepExpanRate)));

  PetscFunctionBegin;
  PetscCheck(a<=a1 && a1<=b1 && b1<=b,PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"Values in the frame vector should be non-decreasing");
  PetscCheck(a1!=b1,PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"The range of wanted eigenvalues cannot be of size zero");
  filterInfo->filterType = 2;      /* mid-pass filter, for interior eigenvalues */
  if (b == b1) {
    PetscCheck(a!=a1,PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"A polynomial filter should not cover all eigenvalues");
    filterInfo->filterType = 1;    /* high-pass filter, for largest eigenvalues */
  } else PetscCheck(a!=a1,PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"filterType==3 for smallest eigenvalues should be pre-converted to filterType==1 for largest eigenvalues");

  /* the following recipe follows Yousef Saad (2005, 2006) with a few minor adaptations / enhancements */
  halfPlateau = 0.5*(b1-a1)*opts->initialPlateau;    /* half length of the "plateau" of the (dual) base filter */
  leftDelta = (b1-a1)*opts->initialShiftStep;        /* initial left shift */
  rightDelta = leftDelta;                            /* initial right shift */
  opts->numGridPoints = PetscMax(opts->numGridPoints,(PetscInt)(2.0*(b-a)/halfPlateau));
  gridSize = (b-a) / (PetscReal)opts->numGridPoints;

  for (i=0;i<6;i++) intv[i] = 0.0;
  if (filterInfo->filterType == 2) {  /* for interior eigenvalues */
    npoints = 6;
    intv[0] = a;
    intv[5] = b;
    /* intv[1], intv[2], intv[3], and intv[4] to be determined */
  } else { /* filterType == 1 (or 3 with conversion), for extreme eigenvalues */
    npoints = 4;
    intv[0] = a;
    intv[3] = b;
    /* intv[1], and intv[2] to be determined */
  }
  z1 = a1 - leftDelta;
  z2 = b1 + rightDelta;
  filterInfo->filterOK = 0;  /* not yet found any OK filter */

  /* allocate matrices */
  PetscCall(PetscMalloc2((polyDeg+2)*(npoints-1),&polyFilter,(2*baseDeg+2)*(npoints-1),&baseFilter));

  /* initialize the intervals, mainly for the case opts->maxOuterIter == 0 */
  intervals[0] = intv[0];
  intervals[3] = intv[3];
  intervals[5] = intv[5];
  intervals[1] = z1;
  if (filterInfo->filterType == 2) {  /* a mid-pass filter for interior eigenvalues */
    intervals[4] = z2;
    c = (a1+b1) / 2.0;
    intervals[2] = c - halfPlateau;
    intervals[3] = c + halfPlateau;
  } else {  /* filterType == 1 (or 3 with conversion) for extreme eigenvalues */
    intervals[2] = z1 + (b1-z1)*opts->transIntervalRatio;
  }

  /* the main loop */
  for (numIter=1;numIter<=opts->maxOuterIter;numIter++) {
    if (z1 <= a) {  /* outer loop updates z1 and z2 */
      z1 = a;
      leftBoundaryMet = PETSC_TRUE;
    }
    if (filterInfo->filterType == 2) {  /* a <= z1 < (a1) */
      if (z2 >= b) {  /* a mid-pass filter for interior eigenvalues */
        z2 = b;
        rightBoundaryMet = PETSC_TRUE;
      }
      /* a <= z1 < c-h < c+h < z2 <= b, where h is halfPlateau */
      /* [z1, c-h] and [c+h, z2] are transition interval */
      intv[1] = z1;
      intv[4] = z2;
      c1 = z1 + halfPlateau;
      intv[2] = z1;           /* i.e. c1 - halfPlateau */
      intv[3] = c1 + halfPlateau;
      PetscCall(FILTLAN_HermiteBaseFilterInChebyshevBasis(baseFilter,intv,6,HighLowFlags,baseDeg));
      PetscCall(FILTLAN_FilteredConjugateResidualPolynomial(polyFilter,baseFilter,2*baseDeg+2,intv,6,opts->intervalWeights,polyDeg));
      /* fc1 = FILTLAN_PiecewisePolynomialEvaluationInChebyshevBasis(polyFilter,polyDeg+2,intv,npoints,b1) - FILTLAN_PiecewisePolynomialEvaluationInChebyshevBasis(polyFilter,polyDeg+2,intv,npoints,a1); */
      c2 = z2 - halfPlateau;
      intv[2] = c2 - halfPlateau;
      intv[3] = z2;           /* i.e. c2 + halfPlateau */
      PetscCall(FILTLAN_HermiteBaseFilterInChebyshevBasis(baseFilter,intv,6,HighLowFlags,baseDeg));
      PetscCall(FILTLAN_FilteredConjugateResidualPolynomial(polyFilter,baseFilter,2*baseDeg+2,intv,6,opts->intervalWeights,polyDeg));
      fc2 = FILTLAN_PiecewisePolynomialEvaluationInChebyshevBasis(polyFilter,polyDeg+2,intv,npoints,b1) - FILTLAN_PiecewisePolynomialEvaluationInChebyshevBasis(polyFilter,polyDeg+2,intv,npoints,a1);
      yLimitGap = PETSC_MAX_REAL;
      ii = opts->maxInnerIter;
      while (ii-- && !(yLimitGap <= opts->yLimitTol)) {
        /* recursive bisection to get c such that p(a1) are p(b1) approximately the same */
        c = (c1+c2) / 2.0;
        intv[2] = c - halfPlateau;
        intv[3] = c + halfPlateau;
        PetscCall(FILTLAN_HermiteBaseFilterInChebyshevBasis(baseFilter,intv,6,HighLowFlags,baseDeg));
        PetscCall(FILTLAN_FilteredConjugateResidualPolynomial(polyFilter,baseFilter,2*baseDeg+2,intv,6,opts->intervalWeights,polyDeg));
        fc = FILTLAN_PiecewisePolynomialEvaluationInChebyshevBasis(polyFilter,polyDeg+2,intv,npoints,b1) - FILTLAN_PiecewisePolynomialEvaluationInChebyshevBasis(polyFilter,polyDeg+2,intv,npoints,a1);
        if (fc*fc2 < 0.0) {
          c1 = c;
          /* fc1 = fc; */
        } else {
          c2 = c;
          fc2 = fc;
        }
        yLimitGap = PetscAbsReal(fc);
      }
    } else {  /* filterType == 1 (or 3 with conversion) for extreme eigenvalues */
      intv[1] = z1;
      intv[2] = z1 + (b1-z1)*opts->transIntervalRatio;
      PetscCall(FILTLAN_HermiteBaseFilterInChebyshevBasis(baseFilter,intv,4,HighLowFlags,baseDeg));
      PetscCall(FILTLAN_FilteredConjugateResidualPolynomial(polyFilter,baseFilter,2*baseDeg+2,intv,4,opts->intervalWeights,polyDeg));
    }
    /* polyFilter contains the coefficients of the polynomial filter which approximates phi(x)
       expanded in the `translated' Chebyshev basis */
    /* psi(x) = 1.0 - phi(x) is the dual base filter approximated by a polynomial in the form x*p(x) */
    yLeftLimit  = 1.0 - FILTLAN_PiecewisePolynomialEvaluationInChebyshevBasis(polyFilter,polyDeg+2,intv,npoints,a1);
    yRightLimit = 1.0 - FILTLAN_PiecewisePolynomialEvaluationInChebyshevBasis(polyFilter,polyDeg+2,intv,npoints,b1);
    yLimit  = (yLeftLimit < yRightLimit) ? yLeftLimit : yRightLimit;
    ySummit = (yLeftLimit > yRightLimit) ? yLeftLimit : yRightLimit;
    x = a1;
    while ((x+=gridSize) < b1) {
      y = 1.0 - FILTLAN_PiecewisePolynomialEvaluationInChebyshevBasis(polyFilter,polyDeg+2,intv,npoints,x);
      if (y < yLimit)  yLimit  = y;
      if (y > ySummit) ySummit = y;
    }
    /* now yLimit is the minimum of x*p(x) for x in [a1, b1] */
    stepLeft  = PETSC_FALSE;
    stepRight = PETSC_FALSE;
    if ((yLimit < yLeftLimit && yLimit < yRightLimit) || yLimit < opts->yBottomLine) {
      /* very bad, step to see what will happen */
      stepLeft = PETSC_TRUE;
      if (filterInfo->filterType == 2) stepRight = PETSC_TRUE;
    } else if (filterInfo->filterType == 2) {
      if (yLeftLimit < yRightLimit) {
        if (yRightLimit-yLeftLimit > opts->yLimitTol) stepLeft = PETSC_TRUE;
      } else if (yLeftLimit-yRightLimit > opts->yLimitTol) stepRight = PETSC_TRUE;
    }
    if (!stepLeft) {
      yLeftBottom = yLeftLimit;
      x = a1;
      while ((x-=gridSize) >= a) {
        y = 1.0 - FILTLAN_PiecewisePolynomialEvaluationInChebyshevBasis(polyFilter,polyDeg+2,intv,npoints,x);
        if (y < yLeftBottom) yLeftBottom = y;
        else if (y > yLeftBottom) break;
      }
      yLeftSummit = yLeftBottom;
      while ((x-=gridSize) >= a) {
        y = 1.0 - FILTLAN_PiecewisePolynomialEvaluationInChebyshevBasis(polyFilter,polyDeg+2,intv,npoints,x);
        if (y > yLeftSummit) {
          yLeftSummit = y;
          if (yLeftSummit > yLimit*opts->yRippleLimit) {
            stepLeft = PETSC_TRUE;
            break;
          }
        }
        if (y < yLeftBottom) yLeftBottom = y;
      }
    }
    if (filterInfo->filterType == 2 && !stepRight) {
      yRightBottom = yRightLimit;
      x = b1;
      while ((x+=gridSize) <= b) {
        y = 1.0 - FILTLAN_PiecewisePolynomialEvaluationInChebyshevBasis(polyFilter,polyDeg+2,intv,npoints,x);
        if (y < yRightBottom) yRightBottom = y;
        else if (y > yRightBottom) break;
      }
      yRightSummit = yRightBottom;
      while ((x+=gridSize) <= b) {
        y = 1.0 - FILTLAN_PiecewisePolynomialEvaluationInChebyshevBasis(polyFilter,polyDeg+2,intv,npoints,x);
        if (y > yRightSummit) {
          yRightSummit = y;
          if (yRightSummit > yLimit*opts->yRippleLimit) {
            stepRight = PETSC_TRUE;
            break;
          }
        }
        if (y < yRightBottom) yRightBottom = y;
      }
    }
    if (!stepLeft && !stepRight) {
      if (filterInfo->filterType == 2) bottom = PetscMin(yLeftBottom,yRightBottom);
      else bottom = yLeftBottom;
      qIndex = 1.0 - (yLimit-bottom) / (ySummit-bottom);
      if (filterInfo->filterOK == 0 || filterInfo->filterQualityIndex < qIndex) {
        /* found the first OK filter or a better filter */
        for (i=0;i<6;i++) intervals[i] = intv[i];
        filterInfo->filterOK           = 1;
        filterInfo->filterQualityIndex = qIndex;
        filterInfo->numIter            = numIter;
        filterInfo->yLimit             = yLimit;
        filterInfo->ySummit            = ySummit;
        filterInfo->numLeftSteps       = numLeftSteps;
        filterInfo->yLeftSummit        = yLeftSummit;
        filterInfo->yLeftBottom        = yLeftBottom;
        if (filterInfo->filterType == 2) {
          filterInfo->yLimitGap        = yLimitGap;
          filterInfo->numRightSteps    = numRightSteps;
          filterInfo->yRightSummit     = yRightSummit;
          filterInfo->yRightBottom     = yRightBottom;
        }
        numMoreLooked = 0;
      } else if (++numMoreLooked == numLookMore) {
        /* filter has been optimal */
        filterInfo->filterOK = 2;
        break;
      }
      /* try stepping further to see whether it can improve */
      stepLeft = PETSC_TRUE;
      if (filterInfo->filterType == 2) stepRight = PETSC_TRUE;
    }
    /* check whether we can really "step" */
    if (leftBoundaryMet) {
      if (filterInfo->filterType == 1 || rightBoundaryMet) break;  /* cannot step further, so break the loop */
      if (stepLeft) {
        /* cannot step left, so try stepping right */
        stepLeft  = PETSC_FALSE;
        stepRight = PETSC_TRUE;
      }
    }
    if (rightBoundaryMet && stepRight) {
      /* cannot step right, so try stepping left */
      stepRight = PETSC_FALSE;
      stepLeft  = PETSC_TRUE;
    }
    /* now "step" */
    if (stepLeft) {
      numLeftSteps++;
      if (filterInfo->filterType == 2) leftDelta *= opts->shiftStepExpanRate; /* expand the step for faster convergence */
      z1 -= leftDelta;
    }
    if (stepRight) {
      numRightSteps++;
      rightDelta *= opts->shiftStepExpanRate;  /* expand the step for faster convergence */
      z2 += rightDelta;
    }
    if (filterInfo->filterType == 2) {
      /* shrink the "plateau" of the (dual) base filter */
      if (stepLeft && stepRight) halfPlateau /= opts->plateauShrinkRate;
      else halfPlateau /= PetscSqrtReal(opts->plateauShrinkRate);
    }
  }
  PetscCheck(filterInfo->filterOK,PETSC_COMM_SELF,PETSC_ERR_CONV_FAILED,"STFILTER cannot get the filter specified; please adjust your filter parameters (e.g. increasing the polynomial degree)");

  filterInfo->totalNumIter = numIter;
  PetscCall(PetscFree2(polyFilter,baseFilter));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* ////////////////////////////////////////////////////////////////////////////
   //    Chebyshev Polynomials
   //////////////////////////////////////////////////////////////////////////// */

/*
   FILTLAN function ExpandNewtonPolynomialInChebyshevBasis

   translate the coefficients of a Newton polynomial to the coefficients
   in a basis of the `translated' (scale-and-shift) Chebyshev polynomials

   input:
   a Newton polynomial defined by vectors a and x:
       P(z) = a(1) + a(2)*(z-x(1)) + a(3)*(z-x(1))*(z-x(2)) + ... + a(n)*(z-x(1))*...*(z-x(n-1))
   the interval [aa,bb] defines the `translated' Chebyshev polynomials S_i(z) = T_i((z-c)/h),
       where c=(aa+bb)/2 and h=(bb-aa)/2, and T_i is the Chebyshev polynomial of the first kind
   note that T_i is defined by T_0(z)=1, T_1(z)=z, and T_i(z)=2*z*T_{i-1}(z)+T_{i-2}(z) for i>=2

   output:
   a vector q containing the Chebyshev coefficients:
       P(z) = q(1)*S_0(z) + q(2)*S_1(z) + ... + q(n)*S_{n-1}(z)
*/
static inline PetscErrorCode FILTLAN_ExpandNewtonPolynomialInChebyshevBasis(PetscInt n,PetscReal aa,PetscReal bb,PetscReal *a,PetscReal *x,PetscReal *q,PetscReal *q2)
{
  PetscInt  m,mm;
  PetscReal *sa,*sx,*sq,*sq2,c,c2,h,h2;

  PetscFunctionBegin;
  sa = a+n;    /* pointers for traversing a and x */
  sx = x+n-1;
  *q = *--sa;  /* set q[0] = a(n) */

  c = (aa+bb)/2.0;
  h = (bb-aa)/2.0;
  h2 = h/2.0;

  for (m=1;m<=n-1;m++) {  /* the main loop for translation */

    /* compute q2[0:m-1] = (c-x[n-m-1])*q[0:m-1] */
    mm = m;
    sq = q;
    sq2 = q2;
    c2 = c-(*--sx);
    while (mm--) *(sq2++) = c2*(*sq++);
    *sq2 = 0.0;         /* q2[m] = 0.0 */
    *(q2+1) += h*(*q);  /* q2[1] = q2[1] + h*q[0] */

    /* compute q2[0:m-2] = q2[0:m-2] + h2*q[1:m-1] */
    mm = m-1;
    sq2 = q2;
    sq = q+1;
    while (mm--) *(sq2++) += h2*(*sq++);

    /* compute q2[2:m] = q2[2:m] + h2*q[1:m-1] */
    mm = m-1;
    sq2 = q2+2;
    sq = q+1;
    while (mm--) *(sq2++) += h2*(*sq++);

    /* compute q[0:m] = q2[0:m] */
    mm = m+1;
    sq2 = q2;
    sq = q;
    while (mm--) *sq++ = *sq2++;
    *q += (*--sa);      /* q[0] = q[0] + p[n-m-1] */
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   FILTLAN function PolynomialEvaluationInChebyshevBasis

   evaluate P(z) at z=z0, where P(z) is a polynomial expanded in a basis of
   the `translated' (i.e. scale-and-shift) Chebyshev polynomials

   input:
   c is a vector of Chebyshev coefficients which defines the polynomial
       P(z) = c(1)*S_0(z) + c(2)*S_1(z) + ... + c(n)*S_{n-1}(z),
   where S_i is the `translated' Chebyshev polynomial S_i((z-c)/h) = T_i(z), with
       c = (intv(j)+intv(j+1)) / 2,  h = (intv(j+1)-intv(j)) / 2
   note that T_i(z) is the Chebyshev polynomial of the first kind,
       T_0(z) = 1, T_1(z) = z, and T_i(z) = 2*z*T_{i-1}(z) - T_{i-2}(z) for i>=2

   output:
   the evaluated value of P(z) at z=z0
*/
static inline PetscReal FILTLAN_PolynomialEvaluationInChebyshevBasis(PetscReal *pp,PetscInt m,PetscInt idx,PetscReal z0,PetscReal aa,PetscReal bb)
{
  PetscInt  ii,deg1;
  PetscReal y,zz,t0,t1,t2,*sc;

  PetscFunctionBegin;
  deg1 = m;
  if (aa==-1.0 && bb==1.0) zz = z0;  /* treat it as a special case to reduce rounding errors */
  else zz = (aa==bb)? 0.0 : -1.0+2.0*(z0-aa)/(bb-aa);

  /* compute y = P(z0), where we utilize the Chebyshev recursion */
  sc = pp+(idx-1)*m;   /* sc points to column idx of pp */
  y = *sc++;
  t1 = 1.0; t2 = zz;
  ii = deg1-1;
  while (ii--) {
    /* Chebyshev recursion: T_0(zz)=1, T_1(zz)=zz, and T_{i+1}(zz) = 2*zz*T_i(zz) + T_{i-1}(zz) for i>=2
       the values of T_{i+1}(zz), T_i(zz), T_{i-1}(zz) are stored in t0, t1, t2, respectively */
    t0 = 2*zz*t1 - t2;
    /* it also works for the base case / the first iteration, where t0 equals 2*zz*1-zz == zz which is T_1(zz) */
    t2 = t1;
    t1 = t0;
    y += (*sc++)*t0;
  }
  PetscFunctionReturn(y);
}

#define basisTranslated PETSC_TRUE
/*
   FILTLAN function PiecewisePolynomialEvaluationInChebyshevBasis

   evaluate P(z) at z=z0, where P(z) is a piecewise polynomial expanded
   in a basis of the (optionally translated, i.e. scale-and-shift) Chebyshev polynomials for each interval

   input:
   intv is a vector which defines the intervals; the j-th interval is [intv(j), intv(j+1))
   pp is a matrix of Chebyshev coefficients which defines a piecewise polynomial P(z)
   in a basis of the `translated' Chebyshev polynomials in each interval
   the polynomial P_j(z) in the j-th interval, i.e. when z is in [intv(j), intv(j+1)), is defined by the j-th column of pp:
       if basisTranslated == false, then
           P_j(z) = pp(1,j)*T_0(z) + pp(2,j)*T_1(z) + ... + pp(n,j)*T_{n-1}(z),
       where T_i(z) is the Chebyshev polynomial of the first kind,
           T_0(z) = 1, T_1(z) = z, and T_i(z) = 2*z*T_{i-1}(z) - T_{i-2}(z) for i>=2
       if basisTranslated == true, then
           P_j(z) = pp(1,j)*S_0(z) + pp(2,j)*S_1(z) + ... + pp(n,j)*S_{n-1}(z),
       where S_i is the `translated' Chebyshev polynomial S_i((z-c)/h) = T_i(z), with
           c = (intv(j)+intv(j+1)) / 2,  h = (intv(j+1)-intv(j)) / 2

   output:
   the evaluated value of P(z) at z=z0

   note that if z0 falls below the first interval, then the polynomial in the first interval will be used to evaluate P(z0)
             if z0 flies over  the last  interval, then the polynomial in the last  interval will be used to evaluate P(z0)
*/
static inline PetscReal FILTLAN_PiecewisePolynomialEvaluationInChebyshevBasis(PetscReal *pp,PetscInt m,PetscReal *intv,PetscInt npoints,PetscReal z0)
{
  PetscReal *sintv,aa,bb,resul;
  PetscInt  idx;

  PetscFunctionBegin;
  /* determine the interval which contains z0 */
  sintv = &intv[1];
  idx = 1;
  if (npoints>2 && z0 >= *sintv) {
    sintv++;
    while (++idx < npoints-1) {
      if (z0 < *sintv) break;
      sintv++;
    }
  }
  /* idx==1 if npoints<=2; otherwise idx satisfies:
         intv(idx) <= z0 < intv(idx+1),  if 2 <= idx <= npoints-2
         z0 < intv(idx+1),               if idx == 1
         intv(idx) <= z0,                if idx == npoints-1
     in addition, sintv points to &intv(idx+1) */
  if (basisTranslated) {
    /* the basis consists of `translated' Chebyshev polynomials */
    /* find the interval of concern, [aa,bb] */
    aa = *(sintv-1);
    bb = *sintv;
    resul = FILTLAN_PolynomialEvaluationInChebyshevBasis(pp,m,idx,z0,aa,bb);
  } else {
    /* the basis consists of standard Chebyshev polynomials, with interval [-1.0,1.0] for integration */
    resul = FILTLAN_PolynomialEvaluationInChebyshevBasis(pp,m,idx,z0,-1.0,1.0);
  }
  PetscFunctionReturn(resul);
}

/*
   FILTLAN function PiecewisePolynomialInnerProductInChebyshevBasis

   compute the weighted inner product of two piecewise polynomials expanded
   in a basis of `translated' (i.e. scale-and-shift) Chebyshev polynomials for each interval

   pp and qq are two matrices of Chebyshev coefficients which define the piecewise polynomials P(z) and Q(z), respectively
   for z in the j-th interval, P(z) equals
       P_j(z) = pp(1,j)*S_0(z) + pp(2,j)*S_1(z) + ... + pp(n,j)*S_{n-1}(z),
   and Q(z) equals
       Q_j(z) = qq(1,j)*S_0(z) + qq(2,j)*S_1(z) + ... + qq(n,j)*S_{n-1}(z),
   where S_i(z) is the `translated' Chebyshev polynomial in that interval,
       S_i((z-c)/h) = T_i(z),  c = (aa+bb)) / 2,  h = (bb-aa) / 2,
   with T_i(z) the Chebyshev polynomial of the first kind
       T_0(z) = 1, T_1(z) = z, and T_i(z) = 2*z*T_{i-1}(z) - T_{i-2}(z) for i>=2

   the (scaled) j-th interval inner product is defined by
       <P_j,Q_j> = (Pi/2)*(pp(1,j)*qq(1,j) + sum_{k} pp(k,j)*qq(k,j)),
   which comes from the property
       <T_0,T_0>=pi, <T_i,T_i>=pi/2 for i>=1, and <T_i,T_j>=0 for i!=j

   the weighted inner product is <P,Q> = sum_{j} intervalWeights(j)*<P_j,Q_j>,
   which is the return value

   note that for unit weights, pass an empty vector of intervalWeights (i.e. of length 0)
*/
static inline PetscReal FILTLAN_PiecewisePolynomialInnerProductInChebyshevBasis(PetscReal *pp,PetscInt prows,PetscInt pcols,PetscInt ldp,PetscReal *qq,PetscInt qrows,PetscInt qcols,PetscInt ldq,const PetscReal *intervalWeights)
{
  PetscInt  nintv,deg1,i,k;
  PetscReal *sp,*sq,ans=0.0,ans2;

  PetscFunctionBegin;
  deg1 = PetscMin(prows,qrows);  /* number of effective coefficients, one more than the effective polynomial degree */
  if (PetscUnlikely(!deg1)) PetscFunctionReturn(0.0);
  nintv = PetscMin(pcols,qcols);  /* number of intervals */

  /* scaled by intervalWeights(i) in the i-th interval (we assume intervalWeights[] are always provided).
     compute ans = sum_{i=1,...,nintv} intervalWeights(i)*[ pp(1,i)*qq(1,i) + sum_{k=1,...,deg} pp(k,i)*qq(k,i) ] */
  for (i=0;i<nintv;i++) {   /* compute ans2 = pp(1,i)*qq(1,i) + sum_{k=1,...,deg} pp(k,i)*qq(k,i) */
    sp = pp+i*ldp;
    sq = qq+i*ldq;
    ans2 = (*sp) * (*sq);  /* the first term pp(1,i)*qq(1,i) is being added twice */
    for (k=0;k<deg1;k++) ans2 += (*sp++)*(*sq++);  /* add pp(k,i)*qq(k,i) */
    ans += ans2*intervalWeights[i];  /* compute ans += ans2*intervalWeights(i) */
  }
  PetscFunctionReturn(ans*PETSC_PI/2.0);
}

/*
   FILTLAN function PiecewisePolynomialInChebyshevBasisMultiplyX

   compute Q(z) = z*P(z), where P(z) and Q(z) are piecewise polynomials expanded
   in a basis of `translated' (i.e. scale-and-shift) Chebyshev polynomials for each interval

   P(z) and Q(z) are stored as matrices of Chebyshev coefficients pp and qq, respectively

   For z in the j-th interval, P(z) equals
       P_j(z) = pp(1,j)*S_0(z) + pp(2,j)*S_1(z) + ... + pp(n,j)*S_{n-1}(z),
   and Q(z) equals
       Q_j(z) = qq(1,j)*S_0(z) + qq(2,j)*S_1(z) + ... + qq(n,j)*S_{n-1}(z),
   where S_i(z) is the `translated' Chebyshev polynomial in that interval,
       S_i((z-c)/h) = T_i(z),  c = (intv(j)+intv(j+1))) / 2,  h = (intv(j+1)-intv(j)) / 2,
   with T_i(z) the Chebyshev polynomial of the first kind
       T_0(z) = 1, T_1(z) = z, and T_i(z) = 2*z*T_{i-1}(z) - T_{i-2}(z) for i>=2

   the returned matrix is qq which represents Q(z) = z*P(z)
*/
static inline PetscErrorCode FILTLAN_PiecewisePolynomialInChebyshevBasisMultiplyX(PetscReal *pp,PetscInt deg1,PetscInt ldp,PetscReal *intv,PetscInt nintv,PetscReal *qq,PetscInt ldq)
{
  PetscInt  i,j;
  PetscReal c,h,h2,tmp,*sp,*sq,*sp2,*sq2;

  PetscFunctionBegin;
  for (j=0;j<nintv;j++) {    /* consider interval between intv(j) and intv(j+1) */
    sp = pp+j*ldp;
    sq = qq+j*ldq;
    sp2 = sp;
    sq2 = sq+1;
    c = 0.5*(intv[j] + intv[j+1]);   /* compute c = (intv(j) + intv(j+1))/2 */
    h = 0.5*(intv[j+1] - intv[j]);   /* compute h = (intv(j+1) - intv(j))/2  and  h2 = h/2 */
    h2 = 0.5*h;
    i = deg1;
    while (i--) *sq++ = c*(*sp++);    /* compute q(1:deg1,j) = c*p(1:deg1,j) */
    *sq++ = 0.0;                      /* set q(deg1+1,j) = 0.0 */
    *(sq2++) += h*(*sp2++);           /* compute q(2,j) = q(2,j) + h*p(1,j) */
    i = deg1-1;
    while (i--) {       /* compute q(3:deg1+1,j) += h2*p(2:deg1,j) and then q(1:deg1-1,j) += h2*p(2:deg1,j) */
      tmp = h2*(*sp2++);
      *(sq2-2) += tmp;
      *(sq2++) += tmp;
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* ////////////////////////////////////////////////////////////////////////////
   //    Conjugate Residual Method in the Polynomial Space
   //////////////////////////////////////////////////////////////////////////// */

/*
    B := alpha*A + beta*B

    A,B are nxk
*/
static inline PetscErrorCode Mat_AXPY_BLAS(PetscInt n,PetscInt k,PetscReal alpha,const PetscReal *A,PetscInt lda,PetscReal beta,PetscReal *B,PetscInt ldb)
{
  PetscInt       i,j;

  PetscFunctionBegin;
  if (beta==(PetscReal)1.0) {
    for (j=0;j<k;j++) for (i=0;i<n;i++) B[i+j*ldb] += alpha*A[i+j*lda];
    PetscCall(PetscLogFlops(2.0*n*k));
  } else {
    for (j=0;j<k;j++) for (i=0;i<n;i++) B[i+j*ldb] = alpha*A[i+j*lda] + beta*B[i+j*ldb];
    PetscCall(PetscLogFlops(3.0*n*k));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   FILTLAN function FilteredConjugateResidualPolynomial

   ** Conjugate Residual Method in the Polynomial Space

   this routine employs a conjugate-residual-type algorithm in polynomial space to minimize ||P(z)-Q(z)||_w,
   where P(z), the base filter, is the input piecewise polynomial, and
         Q(z) is the output polynomial satisfying Q(0)==1, i.e. the constant term of Q(z) is 1
   niter is the number of conjugate-residual iterations; therefore, the degree of Q(z) is up to niter+1
   both P(z) and Q(z) are expanded in the `translated' (scale-and-shift) Chebyshev basis for each interval,
   and presented as matrices of Chebyshev coefficients, denoted by pp and qq, respectively

   input:
   intv is a vector which defines the intervals; the j-th interval is [intv(j),intv(j+1))
   w is a vector of Chebyshev weights; the weight of j-th interval is w(j)
       the interval weights define the inner product of two continuous functions and then
       the derived w-norm ||P(z)-Q(z)||_w
   pp is a matrix of Chebyshev coefficients which defines the piecewise polynomial P(z)
   to be specific, for z in [intv(j), intv(j+1)), P(z) equals
       P_j(z) = pp(1,j)*S_0(z) + pp(2,j)*S_1(z) + ... + pp(niter+2,j)*S_{niter+1}(z),
   where S_i(z) is the `translated' Chebyshev polynomial in that interval,
       S_i((z-c)/h) = T_i(z),  c = (intv(j)+intv(j+1))) / 2,  h = (intv(j+1)-intv(j)) / 2,
   with T_i(z) the Chebyshev polynomial of the first kind,
       T_0(z) = 1, T_1(z) = z, and T_i(z) = 2*z*T_{i-1}(z) - T_{i-2}(z) for i>=2

   output:
   the return matrix, denoted by qq, represents a polynomial Q(z) with degree up to 1+niter
   and satisfying Q(0)==1, such that ||P(z))-Q(z)||_w is minimized
   this polynomial Q(z) is expanded in the `translated' Chebyshev basis for each interval
   to be precise, considering z in [intv(j), intv(j+1)), Q(z) equals
       Q_j(z) = qq(1,j)*S_0(z) + qq(2,j)*S_1(z) + ... + qq(niter+2,j)*S_{niter+1}(z)

   note:
   1. since Q(0)==1, P(0)==1 is expected; if P(0)!=1, one can translate P(z)
      for example, if P(0)==0, one can use 1-P(z) as input instead of P(z)
   2. typically, the base filter, defined by pp and intv, is from Hermite interpolation
      in intervals [intv(j),intv(j+1)) for j=1,...,nintv, with nintv the number of intervals
*/
static PetscErrorCode FILTLAN_FilteredConjugateResidualPolynomial(PetscReal *cpol,PetscReal *baseFilter,PetscInt nbase,PetscReal *intv,PetscInt m,PetscReal *intervalWeights,PetscInt niter)
{
  PetscInt       i,j,srpol,scpol,sarpol,sppol,sappol,ld,nintv;
  PetscReal      rho,rho0,rho1,den,bet,alp,alp0,*ppol,*rpol,*appol,*arpol;

  PetscFunctionBegin;
  nintv = m-1;
  ld = niter+2;  /* leading dimension */
  PetscCall(PetscCalloc4(ld*nintv,&ppol,ld*nintv,&rpol,ld*nintv,&appol,ld*nintv,&arpol));
  PetscCall(PetscArrayzero(cpol,ld*nintv));
  /* initialize polynomial ppol to be 1 (i.e. multiplicative identity) in all intervals */
  sppol = 2;
  srpol = 2;
  scpol = 2;
  for (j=0;j<nintv;j++) {
    ppol[j*ld] = 1.0;
    rpol[j*ld] = 1.0;
    cpol[j*ld] = 1.0;
  }
  /* ppol is the initial p-polynomial (corresponding to the A-conjugate vector p in CG)
     rpol is the r-polynomial (corresponding to the residual vector r in CG)
     cpol is the "corrected" residual polynomial (result of this function) */
  sappol = 3;
  sarpol = 3;
  PetscCall(FILTLAN_PiecewisePolynomialInChebyshevBasisMultiplyX(ppol,sppol,ld,intv,nintv,appol,ld));
  for (i=0;i<3;i++) for (j=0;j<nintv;j++) arpol[i+j*ld] = appol[i+j*ld];
  rho = FILTLAN_PiecewisePolynomialInnerProductInChebyshevBasis(rpol,srpol,nintv,ld,arpol,sarpol,nintv,ld,intervalWeights);
  for (i=0;i<niter;i++) {
    den = FILTLAN_PiecewisePolynomialInnerProductInChebyshevBasis(appol,sappol,nintv,ld,appol,sappol,nintv,ld,intervalWeights);
    alp0 = rho/den;
    rho1 = FILTLAN_PiecewisePolynomialInnerProductInChebyshevBasis(baseFilter,nbase,nintv,nbase,appol,sappol,nintv,ld,intervalWeights);
    alp  = (rho-rho1)/den;
    srpol++;
    scpol++;
    PetscCall(Mat_AXPY_BLAS(srpol,nintv,-alp0,appol,ld,1.0,rpol,ld));
    PetscCall(Mat_AXPY_BLAS(scpol,nintv,-alp,appol,ld,1.0,cpol,ld));
    if (i+1 == niter) break;
    sarpol++;
    PetscCall(FILTLAN_PiecewisePolynomialInChebyshevBasisMultiplyX(rpol,srpol,ld,intv,nintv,arpol,ld));
    rho0 = rho;
    rho = FILTLAN_PiecewisePolynomialInnerProductInChebyshevBasis(rpol,srpol,nintv,ld,arpol,sarpol,nintv,ld,intervalWeights);
    bet  = rho / rho0;
    sppol++;
    sappol++;
    PetscCall(Mat_AXPY_BLAS(sppol,nintv,1.0,rpol,ld,bet,ppol,ld));
    PetscCall(Mat_AXPY_BLAS(sappol,nintv,1.0,arpol,ld,bet,appol,ld));
  }
  PetscCall(PetscFree4(ppol,rpol,appol,arpol));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   FILTLAN function FilteredConjugateResidualMatrixPolynomialVectorProduct

   this routine employs a conjugate-residual-type algorithm in polynomial space to compute
   x = x0 + s(A)*r0 with r0 = b - A*x0, such that ||(1-P(z))-z*s(z)||_w is minimized, where
   P(z) is a given piecewise polynomial, called the base filter,
   s(z) is a polynomial of degree up to niter, the number of conjugate-residual iterations,
   and b and x0 are given vectors

   note that P(z) is expanded in the `translated' (scale-and-shift) Chebyshev basis for each interval,
   and presented as a matrix of Chebyshev coefficients pp

   input:
   A is a sparse matrix
   x0, b are vectors
   niter is the number of conjugate-residual iterations
   intv is a vector which defines the intervals; the j-th interval is [intv(j),intv(j+1))
   w is a vector of Chebyshev weights; the weight of j-th interval is w(j)
       the interval weights define the inner product of two continuous functions and then
       the derived w-norm ||P(z)-Q(z)||_w
   pp is a matrix of Chebyshev coefficients which defines the piecewise polynomial P(z)
   to be specific, for z in [intv(j), intv(j+1)), P(z) equals
       P_j(z) = pp(1,j)*S_0(z) + pp(2,j)*S_1(z) + ... + pp(niter+2,j)*S_{niter+1}(z),
   where S_i(z) is the `translated' Chebyshev polynomial in that interval,
       S_i((z-c)/h) = T_i(z),  c = (intv(j)+intv(j+1))) / 2,  h = (intv(j+1)-intv(j)) / 2,
   with T_i(z) the Chebyshev polynomial of the first kind,
       T_0(z) = 1, T_1(z) = z, and T_i(z) = 2*z*T_{i-1}(z) - T_{i-2}(z) for i>=2
   tol is the tolerance; if the residual polynomial in z-norm is dropped by a factor lower
       than tol, then stop the conjugate-residual iteration

   output:
   the return vector is x = x0 + s(A)*r0 with r0 = b - A*x0, such that ||(1-P(z))-z*s(z)||_w is minimized,
   subject to that s(z) is a polynomial of degree up to niter, where P(z) is the base filter
   in short, z*s(z) approximates 1-P(z)

   note:
   1. since z*s(z) approximates 1-P(z), P(0)==1 is expected; if P(0)!=1, one can translate P(z)
      for example, if P(0)==0, one can use 1-P(z) as input instead of P(z)
   2. typically, the base filter, defined by pp and intv, is from Hermite interpolation
      in intervals [intv(j),intv(j+1)) for j=1,...,nintv, with nintv the number of intervals
   3. a common application is to compute R(A)*b, where R(z) approximates 1-P(z)
      in this case, one can set x0 = 0 and then the return vector is x = s(A)*b, where
      z*s(z) approximates 1-P(z); therefore, A*x is the wanted R(A)*b
*/
static PetscErrorCode FILTLAN_FilteredConjugateResidualMatrixPolynomialVectorProduct(Mat A,Vec b,Vec x,PetscReal *baseFilter,PetscInt nbase,PetscReal *intv,PetscInt nintv,PetscReal *intervalWeights,PetscInt niter,Vec *work)
{
  PetscInt       i,j,srpol,scpol,sarpol,sppol,sappol,ld;
  PetscReal      rho,rho0,rho00,rho1,den,bet,alp,alp0,*cpol,*ppol,*rpol,*appol,*arpol,tol=0.0;
  Vec            r=work[0],p=work[1],ap=work[2],w=work[3];
  PetscScalar    alpha;

  PetscFunctionBegin;
  ld = niter+3;  /* leading dimension */
  PetscCall(PetscCalloc5(ld*nintv,&ppol,ld*nintv,&rpol,ld*nintv,&cpol,ld*nintv,&appol,ld*nintv,&arpol));
  /* initialize polynomial ppol to be 1 (i.e. multiplicative identity) in all intervals */
  sppol = 2;
  srpol = 2;
  scpol = 2;
  for (j=0;j<nintv;j++) {
    ppol[j*ld] = 1.0;
    rpol[j*ld] = 1.0;
    cpol[j*ld] = 1.0;
  }
  /* ppol is the initial p-polynomial (corresponding to the A-conjugate vector p in CG)
     rpol is the r-polynomial (corresponding to the residual vector r in CG)
     cpol is the "corrected" residual polynomial */
  sappol = 3;
  sarpol = 3;
  PetscCall(FILTLAN_PiecewisePolynomialInChebyshevBasisMultiplyX(ppol,sppol,ld,intv,nintv,appol,ld));
  for (i=0;i<3;i++) for (j=0;j<nintv;j++) arpol[i+j*ld] = appol[i+j*ld];
  rho00 = FILTLAN_PiecewisePolynomialInnerProductInChebyshevBasis(rpol,srpol,nintv,ld,arpol,sarpol,nintv,ld,intervalWeights);
  rho = rho00;

  /* corrected CR in vector space */
  /* we assume x0 is always 0 */
  PetscCall(VecSet(x,0.0));
  PetscCall(VecCopy(b,r));     /* initial residual r = b-A*x0 */
  PetscCall(VecCopy(r,p));     /* p = r */
  PetscCall(MatMult(A,p,ap));  /* ap = A*p */

  for (i=0;i<niter;i++) {
    /* iteration in the polynomial space */
    den = FILTLAN_PiecewisePolynomialInnerProductInChebyshevBasis(appol,sappol,nintv,ld,appol,sappol,nintv,ld,intervalWeights);
    alp0 = rho/den;
    rho1 = FILTLAN_PiecewisePolynomialInnerProductInChebyshevBasis(baseFilter,nbase,nintv,nbase,appol,sappol,nintv,ld,intervalWeights);
    alp  = (rho-rho1)/den;
    srpol++;
    scpol++;
    PetscCall(Mat_AXPY_BLAS(srpol,nintv,-alp0,appol,ld,1.0,rpol,ld));
    PetscCall(Mat_AXPY_BLAS(scpol,nintv,-alp,appol,ld,1.0,cpol,ld));
    sarpol++;
    PetscCall(FILTLAN_PiecewisePolynomialInChebyshevBasisMultiplyX(rpol,srpol,ld,intv,nintv,arpol,ld));
    rho0 = rho;
    rho = FILTLAN_PiecewisePolynomialInnerProductInChebyshevBasis(rpol,srpol,nintv,ld,arpol,sarpol,nintv,ld,intervalWeights);

    /* update x in the vector space */
    alpha = alp;
    PetscCall(VecAXPY(x,alpha,p));   /* x += alp*p */
    if (rho < tol*rho00) break;

    /* finish the iteration in the polynomial space */
    bet = rho / rho0;
    sppol++;
    sappol++;
    PetscCall(Mat_AXPY_BLAS(sppol,nintv,1.0,rpol,ld,bet,ppol,ld));
    PetscCall(Mat_AXPY_BLAS(sappol,nintv,1.0,arpol,ld,bet,appol,ld));

    /* finish the iteration in the vector space */
    alpha = -alp0;
    PetscCall(VecAXPY(r,alpha,ap));    /* r -= alp0*ap */
    alpha = bet;
    PetscCall(VecAYPX(p,alpha,r));     /* p = r + bet*p */
    PetscCall(MatMult(A,r,w));         /* ap = A*r + bet*ap */
    PetscCall(VecAYPX(ap,alpha,w));
  }
  PetscCall(PetscFree5(ppol,rpol,cpol,appol,arpol));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   Gateway to FILTLAN for evaluating y=p(A)*x
*/
static PetscErrorCode MatMult_FILTLAN(Mat A,Vec x,Vec y)
{
  ST             st;
  ST_FILTER      *ctx;
  FILTLAN_CTX    filtlan;
  PetscInt       npoints;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(A,&st));
  ctx = (ST_FILTER*)st->data;
  filtlan = (FILTLAN_CTX)ctx->data;
  npoints = (filtlan->filterInfo->filterType == 2)? 6: 4;
  PetscCall(FILTLAN_FilteredConjugateResidualMatrixPolynomialVectorProduct(ctx->T,x,y,filtlan->baseFilter,2*filtlan->baseDegree+2,filtlan->intervals,npoints-1,filtlan->opts->intervalWeights,ctx->polyDegree,st->work));
  PetscCall(VecCopy(y,st->work[0]));
  PetscCall(MatMult(ctx->T,st->work[0],y));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* Block version of FILTLAN_FilteredConjugateResidualMatrixPolynomialVectorProduct */
static PetscErrorCode FILTLAN_FilteredConjugateResidualMatrixPolynomialVectorProductBlock(Mat A,Mat B,Mat C,PetscReal *baseFilter,PetscInt nbase,PetscReal *intv,PetscInt nintv,PetscReal *intervalWeights,PetscInt niter,Vec *work,Mat R,Mat P,Mat AP,Mat W)
{
  PetscInt       i,j,srpol,scpol,sarpol,sppol,sappol,ld;
  PetscReal      rho,rho0,rho00,rho1,den,bet,alp,alp0,*cpol,*ppol,*rpol,*appol,*arpol,tol=0.0;
  PetscScalar    alpha;

  PetscFunctionBegin;
  ld = niter+3;  /* leading dimension */
  PetscCall(PetscCalloc5(ld*nintv,&ppol,ld*nintv,&rpol,ld*nintv,&cpol,ld*nintv,&appol,ld*nintv,&arpol));
  /* initialize polynomial ppol to be 1 (i.e. multiplicative identity) in all intervals */
  sppol = 2;
  srpol = 2;
  scpol = 2;
  for (j=0;j<nintv;j++) {
    ppol[j*ld] = 1.0;
    rpol[j*ld] = 1.0;
    cpol[j*ld] = 1.0;
  }
  /* ppol is the initial p-polynomial (corresponding to the A-conjugate vector p in CG)
     rpol is the r-polynomial (corresponding to the residual vector r in CG)
     cpol is the "corrected" residual polynomial */
  sappol = 3;
  sarpol = 3;
  PetscCall(FILTLAN_PiecewisePolynomialInChebyshevBasisMultiplyX(ppol,sppol,ld,intv,nintv,appol,ld));
  for (i=0;i<3;i++) for (j=0;j<nintv;j++) arpol[i+j*ld] = appol[i+j*ld];
  rho00 = FILTLAN_PiecewisePolynomialInnerProductInChebyshevBasis(rpol,srpol,nintv,ld,arpol,sarpol,nintv,ld,intervalWeights);
  rho = rho00;

  /* corrected CR in vector space */
  /* we assume x0 is always 0 */
  PetscCall(MatZeroEntries(C));
  PetscCall(MatCopy(B,R,SAME_NONZERO_PATTERN));     /* initial residual r = b-A*x0 */
  PetscCall(MatCopy(R,P,SAME_NONZERO_PATTERN));     /* p = r */
  PetscCall(MatMatMult(A,P,MAT_REUSE_MATRIX,PETSC_DEFAULT,&AP));  /* ap = A*p */

  for (i=0;i<niter;i++) {
    /* iteration in the polynomial space */
    den = FILTLAN_PiecewisePolynomialInnerProductInChebyshevBasis(appol,sappol,nintv,ld,appol,sappol,nintv,ld,intervalWeights);
    alp0 = rho/den;
    rho1 = FILTLAN_PiecewisePolynomialInnerProductInChebyshevBasis(baseFilter,nbase,nintv,nbase,appol,sappol,nintv,ld,intervalWeights);
    alp  = (rho-rho1)/den;
    srpol++;
    scpol++;
    PetscCall(Mat_AXPY_BLAS(srpol,nintv,-alp0,appol,ld,1.0,rpol,ld));
    PetscCall(Mat_AXPY_BLAS(scpol,nintv,-alp,appol,ld,1.0,cpol,ld));
    sarpol++;
    PetscCall(FILTLAN_PiecewisePolynomialInChebyshevBasisMultiplyX(rpol,srpol,ld,intv,nintv,arpol,ld));
    rho0 = rho;
    rho = FILTLAN_PiecewisePolynomialInnerProductInChebyshevBasis(rpol,srpol,nintv,ld,arpol,sarpol,nintv,ld,intervalWeights);

    /* update x in the vector space */
    alpha = alp;
    PetscCall(MatAXPY(C,alpha,P,SAME_NONZERO_PATTERN));   /* x += alp*p */
    if (rho < tol*rho00) break;

    /* finish the iteration in the polynomial space */
    bet = rho / rho0;
    sppol++;
    sappol++;
    PetscCall(Mat_AXPY_BLAS(sppol,nintv,1.0,rpol,ld,bet,ppol,ld));
    PetscCall(Mat_AXPY_BLAS(sappol,nintv,1.0,arpol,ld,bet,appol,ld));

    /* finish the iteration in the vector space */
    alpha = -alp0;
    PetscCall(MatAXPY(R,alpha,AP,SAME_NONZERO_PATTERN));    /* r -= alp0*ap */
    alpha = bet;
    PetscCall(MatAYPX(P,alpha,R,SAME_NONZERO_PATTERN));     /* p = r + bet*p */
    PetscCall(MatMatMult(A,R,MAT_REUSE_MATRIX,PETSC_DEFAULT,&W));         /* ap = A*r + bet*ap */
    PetscCall(MatAYPX(AP,alpha,W,SAME_NONZERO_PATTERN));
  }
  PetscCall(PetscFree5(ppol,rpol,cpol,appol,arpol));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   Gateway to FILTLAN for evaluating C=p(A)*B
*/
static PetscErrorCode MatMatMult_FILTLAN(Mat A,Mat B,Mat C,void *pctx)
{
  ST             st;
  ST_FILTER      *ctx;
  FILTLAN_CTX    filtlan;
  PetscInt       i,m1,m2,npoints;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(A,&st));
  ctx = (ST_FILTER*)st->data;
  filtlan = (FILTLAN_CTX)ctx->data;
  npoints = (filtlan->filterInfo->filterType == 2)? 6: 4;
  if (ctx->nW) {  /* check if work matrices must be resized */
    PetscCall(MatGetSize(B,NULL,&m1));
    PetscCall(MatGetSize(ctx->W[0],NULL,&m2));
    if (m1!=m2) {
      PetscCall(MatDestroyMatrices(ctx->nW,&ctx->W));
      ctx->nW = 0;
    }
  }
  if (!ctx->nW) {  /* allocate work matrices */
    ctx->nW = 4;
    PetscCall(PetscMalloc1(ctx->nW,&ctx->W));
    for (i=0;i<ctx->nW;i++) PetscCall(MatDuplicate(B,MAT_DO_NOT_COPY_VALUES,&ctx->W[i]));
  }
  PetscCall(FILTLAN_FilteredConjugateResidualMatrixPolynomialVectorProductBlock(ctx->T,B,ctx->W[0],filtlan->baseFilter,2*filtlan->baseDegree+2,filtlan->intervals,npoints-1,filtlan->opts->intervalWeights,ctx->polyDegree,st->work,C,ctx->W[1],ctx->W[2],ctx->W[3]));
  PetscCall(MatMatMult(ctx->T,ctx->W[0],MAT_REUSE_MATRIX,PETSC_DEFAULT,&C));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   FILTLAN function PolynomialFilterInterface::setFilter

   Creates the shifted (and scaled) matrix and the base filter P(z).
   M is a shell matrix whose MatMult() applies the filter.
*/
static PetscErrorCode STComputeOperator_Filter_FILTLAN(ST st,Mat *G)
{
  ST_FILTER      *ctx = (ST_FILTER*)st->data;
  FILTLAN_CTX    filtlan = (FILTLAN_CTX)ctx->data;
  PetscInt       i,npoints,n,m,N,M;
  PetscReal      frame2[4];
  PetscScalar    alpha;
  const PetscInt HighLowFlags[5] = { 1, -1, 0, -1, 1 };

  PetscFunctionBegin;
  PetscCall(STSetWorkVecs(st,4));
  filtlan->frame[0] = ctx->left;
  filtlan->frame[1] = ctx->inta;
  filtlan->frame[2] = ctx->intb;
  filtlan->frame[3] = ctx->right;
  if (filtlan->frame[0] == filtlan->frame[1]) {  /* low pass filter, convert it to high pass filter */
    /* T = frame[3]*eye(n) - A */
    PetscCall(MatDestroy(&ctx->T));
    PetscCall(MatDuplicate(st->A[0],MAT_COPY_VALUES,&ctx->T));
    PetscCall(MatScale(ctx->T,-1.0));
    alpha = filtlan->frame[3];
    PetscCall(MatShift(ctx->T,alpha));
    for (i=0;i<4;i++) frame2[i] = filtlan->frame[3] - filtlan->frame[3-i];
  } else {  /* it can be a mid-pass filter or a high-pass filter */
    if (filtlan->frame[0] == 0.0) {
      PetscCall(PetscObjectReference((PetscObject)st->A[0]));
      PetscCall(MatDestroy(&ctx->T));
      ctx->T = st->A[0];
      for (i=0;i<4;i++) frame2[i] = filtlan->frame[i];
    } else {
      /* T = A - frame[0]*eye(n) */
      PetscCall(MatDestroy(&ctx->T));
      PetscCall(MatDuplicate(st->A[0],MAT_COPY_VALUES,&ctx->T));
      alpha = -filtlan->frame[0];
      PetscCall(MatShift(ctx->T,alpha));
      for (i=0;i<4;i++) frame2[i] = filtlan->frame[i] - filtlan->frame[0];
    }
  }

  /* no need to recompute filter if the parameters did not change */
  if (st->state==ST_STATE_INITIAL || ctx->filtch) {
    PetscCall(FILTLAN_GetIntervals(filtlan->intervals,frame2,ctx->polyDegree,filtlan->baseDegree,filtlan->opts,filtlan->filterInfo));
    /* translate the intervals back */
    if (filtlan->frame[0] == filtlan->frame[1]) {  /* low pass filter, convert it to high pass filter */
      for (i=0;i<4;i++) filtlan->intervals2[i] = filtlan->frame[3] - filtlan->intervals[3-i];
    } else {  /* it can be a mid-pass filter or a high-pass filter */
      if (filtlan->frame[0] == 0.0) {
        for (i=0;i<6;i++) filtlan->intervals2[i] = filtlan->intervals[i];
      } else {
        for (i=0;i<6;i++) filtlan->intervals2[i] = filtlan->intervals[i] + filtlan->frame[0];
      }
    }

    npoints = (filtlan->filterInfo->filterType == 2)? 6: 4;
    PetscCall(PetscFree(filtlan->baseFilter));
    PetscCall(PetscMalloc1((2*filtlan->baseDegree+2)*(npoints-1),&filtlan->baseFilter));
    PetscCall(FILTLAN_HermiteBaseFilterInChebyshevBasis(filtlan->baseFilter,filtlan->intervals,npoints,HighLowFlags,filtlan->baseDegree));
    PetscCall(PetscInfo(st,"Computed value of yLimit = %g\n",(double)filtlan->filterInfo->yLimit));
  }
  ctx->filtch = PETSC_FALSE;

  /* create shell matrix*/
  if (!*G) {
    PetscCall(MatGetSize(ctx->T,&N,&M));
    PetscCall(MatGetLocalSize(ctx->T,&n,&m));
    PetscCall(MatCreateShell(PetscObjectComm((PetscObject)st),n,m,N,M,st,G));
    PetscCall(MatShellSetOperation(*G,MATOP_MULT,(PetscErrorCodeFn*)MatMult_FILTLAN));
    PetscCall(MatShellSetMatProductOperation(*G,MATPRODUCT_AB,NULL,MatMatMult_FILTLAN,NULL,MATDENSE,MATDENSE));
    PetscCall(MatShellSetMatProductOperation(*G,MATPRODUCT_AB,NULL,MatMatMult_FILTLAN,NULL,MATDENSECUDA,MATDENSECUDA));
    PetscCall(MatShellSetMatProductOperation(*G,MATPRODUCT_AB,NULL,MatMatMult_FILTLAN,NULL,MATDENSEHIP,MATDENSEHIP));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode STFilterGetThreshold_Filter_FILTLAN(ST st,PetscReal *gamma)
{
  ST_FILTER   *ctx = (ST_FILTER*)st->data;
  FILTLAN_CTX filtlan = (FILTLAN_CTX)ctx->data;

  PetscFunctionBegin;
  *gamma = filtlan->filterInfo->yLimit;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode STDestroy_Filter_FILTLAN(ST st)
{
  ST_FILTER   *ctx = (ST_FILTER*)st->data;
  FILTLAN_CTX filtlan = (FILTLAN_CTX)ctx->data;

  PetscFunctionBegin;
  PetscCall(PetscFree(filtlan->opts));
  PetscCall(PetscFree(filtlan->filterInfo));
  PetscCall(PetscFree(filtlan->baseFilter));
  PetscCall(PetscFree(filtlan));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode STCreate_Filter_FILTLAN(ST st)
{
  ST_FILTER   *ctx = (ST_FILTER*)st->data;
  FILTLAN_CTX filtlan;
  FILTLAN_IOP iop;
  FILTLAN_PFI pfi;

  PetscFunctionBegin;
  PetscCall(PetscNew(&filtlan));
  ctx->data = (void*)filtlan;
  filtlan->baseDegree = 10;

  PetscCall(PetscNew(&iop));
  filtlan->opts           = iop;
  iop->intervalWeights[0] = 100.0;
  iop->intervalWeights[1] = 1.0;
  iop->intervalWeights[2] = 1.0;
  iop->intervalWeights[3] = 1.0;
  iop->intervalWeights[4] = 100.0;
  iop->transIntervalRatio = 0.6;
  iop->reverseInterval    = PETSC_FALSE;
  iop->initialPlateau     = 0.1;
  iop->plateauShrinkRate  = 1.5;
  iop->initialShiftStep   = 0.01;
  iop->shiftStepExpanRate = 1.5;
  iop->maxInnerIter       = 30;
  iop->yLimitTol          = 0.0001;
  iop->maxOuterIter       = 50;
  iop->numGridPoints      = 200;
  iop->yBottomLine        = 0.001;
  iop->yRippleLimit       = 0.8;

  PetscCall(PetscNew(&pfi));
  filtlan->filterInfo     = pfi;

  ctx->computeoperator = STComputeOperator_Filter_FILTLAN;
  ctx->getthreshold    = STFilterGetThreshold_Filter_FILTLAN;
  ctx->destroy         = STDestroy_Filter_FILTLAN;
  PetscFunctionReturn(PETSC_SUCCESS);
}
