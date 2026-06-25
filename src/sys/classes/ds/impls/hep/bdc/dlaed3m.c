/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   BDC - Block-divide and conquer (see description in README file)
*/

#include <slepc/private/dsimpl.h>
#include <slepcblaslapack.h>

PetscErrorCode BDC_dlaed3m_(const char *jobz,const char *defl,PetscBLASInt k,PetscBLASInt n,
        PetscBLASInt n1,PetscReal *d,PetscReal *q,PetscBLASInt ldq,
        PetscReal rho,PetscReal *dlamda,PetscReal *q2,PetscBLASInt *indx,
        PetscBLASInt *ctot,PetscReal *w,PetscReal *s,PetscBLASInt *oinfo,
        PetscBLASInt jobz_len,PetscBLASInt defl_len)
{
/*  -- Routine written in LAPACK version 3.0 style -- */
/* *************************************************** */
/*     Written by */
/*     Michael Moldaschl and Wilfried Gansterer */
/*     University of Vienna */
/*     last modification: March 16, 2014 */

/*     Small adaptations of original code written by */
/*     Wilfried Gansterer and Bob Ward, */
/*     Department of Computer Science, University of Tennessee */
/*     see https://doi.org/10.1137/S1064827501399432 */
/* *************************************************** */

/*  Purpose */
/*  ======= */

/*  DLAED3M finds the roots of the secular equation, as defined by the */
/*  values in D, W, and RHO, between 1 and K.  It makes the */
/*  appropriate calls to DLAED4 and then updates the eigenvectors by */
/*  multiplying the matrix of eigenvectors of the pair of eigensystems */
/*  being combined by the matrix of eigenvectors of the K-by-K system */
/*  which is solved here. */

/*  This code makes very mild assumptions about floating point */
/*  arithmetic. It will work on machines with a guard digit in */
/*  add/subtract, or on those binary machines without guard digits */
/*  which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2. */
/*  It could conceivably fail on hexadecimal or decimal machines */
/*  without guard digits, but we know of none. */

/*  Arguments */
/*  ========= */

/*  JOBZ    (input) CHARACTER*1 */
/*          = 'N':  Do not accumulate eigenvectors (not implemented); */
/*          = 'D':  Do accumulate eigenvectors in the divide-and-conquer */
/*                  process. */

/*  DEFL    (input) CHARACTER*1 */
/*          = '0':  No deflation happened in DSRTDF */
/*          = '1':  Some deflation happened in DSRTDF (and therefore some */
/*                  Givens rotations need to be applied to the computed */
/*                  eigenvector matrix Q) */

/*  K       (input) INTEGER */
/*          The number of terms in the rational function to be solved by */
/*          DLAED4. 0 <= K <= N. */

/*  N       (input) INTEGER */
/*          The number of rows and columns in the Q matrix. */
/*          N >= K (deflation may result in N>K). */

/*  N1      (input) INTEGER */
/*          The location of the last eigenvalue in the leading submatrix. */
/*          min(1,N) <= N1 <= max(1,N-1). */

/*  D       (output) DOUBLE PRECISION array, dimension (N) */
/*          D(I) contains the updated eigenvalues for */
/*          1 <= I <= K. */

/*  Q       (output) DOUBLE PRECISION array, dimension (LDQ,N) */
/*          Initially the first K columns are used as workspace. */
/*          On output the columns 1 to K contain */
/*          the updated eigenvectors. */

/*  LDQ     (input) INTEGER */
/*          The leading dimension of the array Q.  LDQ >= max(1,N). */

/*  RHO     (input) DOUBLE PRECISION */
/*          The value of the parameter in the rank one update equation. */
/*          RHO >= 0 required. */

/*  DLAMDA  (input/output) DOUBLE PRECISION array, dimension (K) */
/*          The first K elements of this array contain the old roots */
/*          of the deflated updating problem.  These are the poles */
/*          of the secular equation. May be changed on output by */
/*          having lowest order bit set to zero on Cray X-MP, Cray Y-MP, */
/*          Cray-2, or Cray C-90, as described above. */

/*  Q2      (input) DOUBLE PRECISION array, dimension (LDQ2, N) */
/*          The first K columns of this matrix contain the non-deflated */
/*          eigenvectors for the split problem. */

/*  INDX    (input) INTEGER array, dimension (N) */
/*          The permutation used to arrange the columns of the deflated */
/*          Q matrix into three groups (see DLAED2). */
/*          The rows of the eigenvectors found by DLAED4 must be likewise */
/*          permuted before the matrix multiply can take place. */

/*  CTOT    (input) INTEGER array, dimension (4) */
/*          A count of the total number of the various types of columns */
/*          in Q, as described in INDX.  The fourth column type is any */
/*          column which has been deflated. */

/*  W       (input/output) DOUBLE PRECISION array, dimension (K) */
/*          The first K elements of this array contain the components */
/*          of the deflation-adjusted updating vector. Destroyed on */
/*          output. */

/*  S       (workspace) DOUBLE PRECISION array, dimension */
/*          (MAX(CTOT(1)+CTOT(2),CTOT(2)+CTOT(3)) + 1)*K */
/*          Will contain parts of the eigenvectors of the repaired matrix */
/*          which will be multiplied by the previously accumulated */
/*          eigenvectors to update the system. This array is a major */
/*          source of workspace requirements ! */

/*  INFO    (output) INTEGER */
/*          = 0:  successful exit. */
/*          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/*          > 0:  if INFO = i, eigenpair i was not computed successfully */

/*  Further Details */
/*  =============== */

/*  Based on code written by */
/*     Wilfried Gansterer and Bob Ward, */
/*     Department of Computer Science, University of Tennessee */
/*  Based on the design of the LAPACK code DLAED3 with small modifications */
/*  (Note that in contrast to the original DLAED3, this routine */
/*  DOES NOT require that N1 <= N/2) */

/*  Based on contributions by */
/*     Jeff Rutter, Computer Science Division, University of California */
/*     at Berkeley, USA */
/*  Modified by Francoise Tisseur, University of Tennessee. */

/*  ===================================================================== */

  PetscReal    temp, done = 1.0, dzero = 0.0;
  PetscBLASInt i, j, n2, n12, ii, n23, iq2, i1, one=1;

  PetscFunctionBegin;
  *oinfo = 0;

  if (k < 0) *oinfo = -3;
  else if (n < k) *oinfo = -4;
  else if (n1 < PetscMin(1,n) || n1 > PetscMax(1,n)) *oinfo = -5;
  else if (ldq < PetscMax(1,n)) *oinfo = -8;
  else if (rho < 0.) *oinfo = -9;
  PetscCheck(!*oinfo,PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"Wrong argument %" PetscBLASInt_FMT " in DLAED3M",-(*oinfo));

  /* Quick return if possible */

  if (k == 0) PetscFunctionReturn(PETSC_SUCCESS);

  /* Modify values DLAMDA(i) to make sure all DLAMDA(i)-DLAMDA(j) can */
  /* be computed with high relative accuracy (barring over/underflow). */
  /* This is a problem on machines without a guard digit in */
  /* add/subtract (Cray XMP, Cray YMP, Cray C 90 and Cray 2). */
  /* The following code replaces DLAMDA(I) by 2*DLAMDA(I)-DLAMDA(I), */
  /* which on any of these machines zeros out the bottommost */
  /* bit of DLAMDA(I) if it is 1; this makes the subsequent */
  /* subtractions DLAMDA(I)-DLAMDA(J) unproblematic when cancellation */
  /* occurs. On binary machines with a guard digit (almost all */
  /* machines) it does not change DLAMDA(I) at all. On hexadecimal */
  /* and decimal machines with a guard digit, it slightly */
  /* changes the bottommost bits of DLAMDA(I). It does not account */
  /* for hexadecimal or decimal machines without guard digits */
  /* (we know of none). We use a subroutine call to compute */
  /* 2*DLAMBDA(I) to prevent optimizing compilers from eliminating */
  /* this code. */

  for (i = 0; i < k; ++i) {
    dlamda[i] = LAPACKlamc3_(&dlamda[i], &dlamda[i]) - dlamda[i];
  }

  for (j = 1; j <= k; ++j) {

    /* ....calling DLAED4 for eigenpair J.... */

    PetscCallLAPACKInfo("LAPACKlaed4",LAPACKlaed4_(&k, &j, dlamda, w, &q[(j-1)*ldq], &rho, &d[j-1], &info));

    if (j < k) {

      /* If the zero finder terminated properly, but the computed */
      /* eigenvalues are not ordered, issue an error statement */
      /* but continue computation. */

      PetscCheck(dlamda[j-1]<dlamda[j],PETSC_COMM_SELF,PETSC_ERR_FP,"DLAMDA(%" PetscBLASInt_FMT ") is greater or equal than DLAMDA(%" PetscBLASInt_FMT ")", j, j+1);
      PetscCheck(d[j-1]>=dlamda[j-1] && d[j-1]<=dlamda[j],PETSC_COMM_SELF,PETSC_ERR_FP,"DLAMDA(%" PetscBLASInt_FMT ") = %g D(%" PetscBLASInt_FMT ") = %g DLAMDA(%" PetscBLASInt_FMT ") = %g", j, (double)dlamda[j-1], j, (double)d[j-1], j+1, (double)dlamda[j]);
    }
  }

  if (k == 1) goto L110;

  if (k == 2) {

    /* permute the components of Q(:,J) (the information returned by DLAED4 */
    /* necessary to construct the eigenvectors) according to the permutation */
    /* stored in INDX, resulting from deflation */

    for (j = 0; j < k; ++j) {
      w[0] = q[0+j*ldq];
      w[1] = q[1+j*ldq];
      ii = indx[0];
      q[0+j*ldq] = w[ii-1];
      ii = indx[1];
      q[1+j*ldq] = w[ii-1];
    }
    goto L110;
  }

  /* ....K.GE.3.... */
  /* Compute updated W (used for computing the eigenvectors corresponding */
  /* to the previously computed eigenvalues). */

  PetscCallBLAS("BLAScopy",BLAScopy_(&k, w, &one, s, &one));

  /* Initialize W(I) = Q(I,I) */

  i1 = ldq + 1;
  PetscCallBLAS("BLAScopy",BLAScopy_(&k, q, &i1, w, &one));
  for (j = 0; j < k; ++j) {
    for (i = 0; i < j; ++i) {
      w[i] *= q[i+j*ldq] / (dlamda[i] - dlamda[j]);
    }
    for (i = j + 1; i < k; ++i) {
      w[i] *= q[i+j*ldq] / (dlamda[i] - dlamda[j]);
    }
  }
  for (i = 0; i < k; ++i) {
    temp = PetscSqrtReal(-w[i]);
    if (temp<0) temp = -temp;
    w[i] =  (s[i] >= 0) ? temp : -temp;
  }

  /* Compute eigenvectors of the modified rank-1 modification (using the */
  /* vector W). */

  for (j = 0; j < k; ++j) {
    for (i = 0; i < k; ++i) {
      s[i] = w[i] / q[i+j*ldq];
    }
    temp = BLASnrm2_(&k, s, &one);
    for (i = 0; i < k; ++i) {

      /* apply the permutation resulting from deflation as stored */
      /* in INDX */

      ii = indx[i];
      q[i+j*ldq] = s[ii-1] / temp;
    }
  }

/* ************************************************************************** */

  /* ....updating the eigenvectors.... */

L110:

  n2 = n - n1;
  n12 = ctot[0] + ctot[1];
  n23 = ctot[1] + ctot[2];
  if (*(unsigned char *)jobz == 'D') {

    /* Compute the updated eigenvectors. (NOTE that every call of */
    /* DGEMM requires three DISTINCT arrays) */

    /* copy Q(CTOT(1)+1:K,1:K) to S */

    for (j=0;j<k;j++) for (i=0;i<n23;i++) s[i+j*n23] = q[ctot[0]+i+j*ldq];
    iq2 = n1 * n12 + 1;

    if (n23 != 0) {

      /* multiply the second part of Q2 (the eigenvectors of the */
      /* lower block) with S and write the result into the lower part of */
      /* Q, i.e., Q(N1+1:N,1:K) */

      PetscCallBLAS("BLASgemm",BLASgemm_("N", "N", &n2, &k, &n23, &done,
                  &q2[iq2-1], &n2, s, &n23, &dzero, &q[n1], &ldq));
    } else {
      for (j=0;j<k;j++) for (i=0;i<n2;i++) q[n1+i+j*ldq] = 0.0;
    }

    /* copy Q(1:CTOT(1)+CTOT(2),1:K) to S */

    for (j=0;j<k;j++) for (i=0;i<n12;i++) s[i+j*n12] = q[i+j*ldq];

    if (n12 != 0) {

      /* multiply the first part of Q2 (the eigenvectors of the */
      /* upper block) with S and write the result into the upper part of */
      /* Q, i.e., Q(1:N1,1:K) */

      PetscCallBLAS("BLASgemm",BLASgemm_("N", "N", &n1, &k, &n12, &done,
                  q2, &n1, s, &n12, &dzero, q, &ldq));
    } else {
      for (j=0;j<k;j++) for (i=0;i<n1;i++) q[i+j*ldq] = 0.0;
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}
