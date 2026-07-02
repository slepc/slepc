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

PetscErrorCode BDC_dsbtdc_(const char *jobz,const char *jobacc,PetscBLASInt n,
        PetscBLASInt nblks,PetscBLASInt *ksizes,PetscReal *d,PetscBLASInt l1d,
        PetscBLASInt l2d,PetscReal *e,PetscBLASInt l1e,PetscBLASInt l2e,PetscReal tol,
        PetscReal tau1,PetscReal tau2,PetscReal *ev,PetscReal *z,PetscBLASInt ldz,
        PetscReal *work,PetscBLASInt lwork,PetscBLASInt *iwork,PetscBLASInt liwork,
        PetscReal *mingap,PetscBLASInt *mingapi,PetscBLASInt *oinfo,
        PetscBLASInt jobz_len,PetscBLASInt jobacc_len)
{
/*  -- Routine written in LAPACK Version 3.0 style -- */
/* *************************************************** */
/*     Written by */
/*     Michael Moldaschl and Wilfried Gansterer */
/*     University of Vienna */
/*     last modification: March 28, 2014 */

/*     Small adaptations of original code written by */
/*     Wilfried Gansterer and Bob Ward, */
/*     Department of Computer Science, University of Tennessee */
/*     see https://doi.org/10.1137/S1064827501399432 */
/* *************************************************** */

/*  Purpose */
/*  ======= */

/*  DSBTDC computes approximations to all eigenvalues and eigenvectors */
/*  of a symmetric block tridiagonal matrix using the divide and */
/*  conquer method with lower rank approximations to the subdiagonal blocks. */

/*  This code makes very mild assumptions about floating point */
/*  arithmetic. It will work on machines with a guard digit in */
/*  add/subtract, or on those binary machines without guard digits */
/*  which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2. */
/*  It could conceivably fail on hexadecimal or decimal machines */
/*  without guard digits, but we know of none.  See DLAED3M for details. */

/*  Arguments */
/*  ========= */

/*  JOBZ    (input) CHARACTER*1 */
/*          = 'N':  Compute eigenvalues only (not implemented); */
/*          = 'D':  Compute eigenvalues and eigenvectors. Eigenvectors */
/*                  are accumulated in the divide-and-conquer process. */

/*  JOBACC  (input) CHARACTER*1 */
/*          = 'A' ("automatic"): The accuracy parameters TAU1 and TAU2 */
/*                               are determined automatically from the */
/*                               parameter TOL according to the analytical */
/*                               bounds. In that case the input values of */
/*                               TAU1 and TAU2 are irrelevant (ignored). */
/*          = 'M' ("manual"): The input values of the accuracy parameters */
/*                            TAU1 and TAU2 are used. In that case the input */
/*                            value of the parameter TOL is irrelevant */
/*                            (ignored). */

/*  N       (input) INTEGER */
/*          The dimension of the symmetric block tridiagonal matrix. */
/*          N >= 1. */

/*  NBLKS   (input) INTEGER, 1 <= NBLKS <= N */
/*          The number of diagonal blocks in the matrix. */

/*  KSIZES  (input) INTEGER array, dimension (NBLKS) */
/*          The dimensions of the square diagonal blocks from top left */
/*          to bottom right.  KSIZES(I) >= 1 for all I, and the sum of */
/*          KSIZES(I) for I = 1 to NBLKS has to be equal to N. */

/*  D       (input) DOUBLE PRECISION array, dimension (L1D,L2D,NBLKS) */
/*          The lower triangular elements of the symmetric diagonal */
/*          blocks of the block tridiagonal matrix. The elements of the top */
/*          left diagonal block, which is of dimension KSIZES(1), have to */
/*          be placed in D(*,*,1); the elements of the next diagonal */
/*          block, which is of dimension KSIZES(2), have to be placed in */
/*          D(*,*,2); etc. */

/*  L1D     (input) INTEGER */
/*          The leading dimension of the array D.  L1D >= max(3,KMAX), */
/*          where KMAX is the dimension of the largest diagonal block, */
/*          i.e.,  KMAX = max_I (KSIZES(I)). */

/*  L2D     (input) INTEGER */
/*          The second dimension of the array D.  L2D >= max(3,KMAX), */
/*          where KMAX is as stated in L1D above. */

/*  E       (input) DOUBLE PRECISION array, dimension (L1E,L2E,NBLKS-1) */
/*          The elements of the subdiagonal blocks of the */
/*          block tridiagonal matrix. The elements of the top left */
/*          subdiagonal block, which is KSIZES(2) x KSIZES(1), have to be */
/*          placed in E(*,*,1); the elements of the next subdiagonal block, */
/*          which is KSIZES(3) x KSIZES(2), have to be placed in E(*,*,2); etc. */
/*          During runtime, the original contents of E(*,*,K) is */
/*          overwritten by the singular vectors and singular values of */
/*          the lower rank representation. */

/*  L1E     (input) INTEGER */
/*          The leading dimension of the array E.  L1E >= max(3,2*KMAX+1), */
/*          where KMAX is as stated in L1D above. The size of L1E enables */
/*          the storage of ALL singular vectors and singular values for */
/*          the corresponding off-diagonal block in E(*,*,K) and therefore */
/*          there are no restrictions on the rank of the approximation */
/*          (only the "natural" restriction */
/*          RANK(K) .LE. MIN(KSIZES(K),KSIZES(K+1))). */

/*  L2E     (input) INTEGER */
/*          The second dimension of the array E.  L2E >= max(3,2*KMAX+1), */
/*          where KMAX is as stated in L1D above. The size of L2E enables */
/*          the storage of ALL singular vectors and singular values for */
/*          the corresponding off-diagonal block in E(*,*,K) and therefore */
/*          there are no restrictions on the rank of the approximation */
/*          (only the "natural" restriction */
/*          RANK(K) .LE. MIN(KSIZES(K),KSIZES(K+1))). */

/*  TOL     (input) DOUBLE PRECISION, TOL.LE.TOLMAX */
/*          User specified tolerance for the residuals of the computed */
/*          eigenpairs. If (JOBACC.EQ.'A') then it is used to determine */
/*          TAU1 and TAU2; ignored otherwise. */
/*          If (TOL.LT.40*EPS .AND. JOBACC.EQ.'A') then TAU1 is set to machine */
/*          epsilon and TAU2 is set to the standard deflation tolerance from */
/*          LAPACK. */

/*  TAU1    (input) DOUBLE PRECISION, TAU1.LE.TOLMAX/2 */
/*          User specified tolerance for determining the rank of the */
/*          lower rank approximations to the off-diagonal blocks. */
/*          The rank for each off-diagonal block is determined such that */
/*          the resulting absolute eigenvalue error is less than or equal */
/*          to TAU1. */
/*          If (JOBACC.EQ.'A') then TAU1 is determined automatically from */
/*             TOL and the input value is ignored. */
/*          If (JOBACC.EQ.'M' .AND. TAU1.LT.20*EPS) then TAU1 is set to */
/*             machine epsilon. */

/*  TAU2    (input) DOUBLE PRECISION, TAU2.LE.TOLMAX/2 */
/*          User specified deflation tolerance for the routine DIBTDC. */
/*          If (1.0D-1.GT.TAU2.GT.20*EPS) then TAU2 is used as */
/*          the deflation tolerance in DSRTDF (EPS is the machine epsilon). */
/*          If (TAU2.LE.20*EPS) then the standard deflation tolerance from */
/*          LAPACK is used as the deflation tolerance in DSRTDF. */
/*          If (JOBACC.EQ.'A') then TAU2 is determined automatically from */
/*             TOL and the input value is ignored. */
/*          If (JOBACC.EQ.'M' .AND. TAU2.LT.20*EPS) then TAU2 is set to */
/*             the standard deflation tolerance from LAPACK. */

/*  EV      (output) DOUBLE PRECISION array, dimension (N) */
/*          If INFO = 0, then EV contains the computed eigenvalues of the */
/*          symmetric block tridiagonal matrix in ascending order. */

/*  Z       (output) DOUBLE PRECISION array, dimension (LDZ,N) */
/*          If (JOBZ.EQ.'D' .AND. INFO = 0) */
/*          then Z contains the orthonormal eigenvectors of the symmetric */
/*          block tridiagonal matrix computed by the routine DIBTDC */
/*          (accumulated in the divide-and-conquer process). */
/*          If (-199 < INFO < -99) then Z contains the orthonormal */
/*          eigenvectors of the symmetric block tridiagonal matrix, */
/*          computed without divide-and-conquer (quick returns). */
/*          Otherwise not referenced. */

/*  LDZ     (input) INTEGER */
/*          The leading dimension of the array Z.  LDZ >= max(1,N). */

/*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK) */

/*  LWORK   (input) INTEGER */
/*          The dimension of the array WORK. */
/*          If NBLKS.EQ.1, then LWORK has to be at least 2N^2+6N+1 */
/*          (for the call of DSYEVD). */
/*          If NBLKS.GE.2 and (JOBZ.EQ.'D') then the absolute minimum */
/*             required for DIBTDC is (N**2 + 3*N). This will not always */
/*             suffice, though, the routine will return a corresponding */
/*             error code and report how much work space was missing (see */
/*             INFO). */
/*          In order to guarantee correct results in all cases where */
/*          NBLKS.GE.2, LWORK must be at least (2*N**2 + 3*N). */

/*  IWORK   (workspace/output) INTEGER array, dimension (LIWORK) */

/*  LIWORK  (input) INTEGER */
/*          The dimension of the array IWORK. */
/*          LIWORK must be at least (5*N + 5*NBLKS - 1) (for DIBTDC) */
/*          Note that this should also suffice for the call of DSYEVD on a */
/*          diagonal block which requires (5*KMAX + 3). */

/*  MINGAP  (output) DOUBLE PRECISION */
/*          The minimum "gap" between the approximate eigenvalues */
/*          computed, i.e., MIN( ABS(EV(I+1)-EV(I)) for I=1,2,..., N-1 */
/*          IF (MINGAP.LE.TOL/10) THEN a warning flag is returned in INFO, */
/*          because the computed eigenvectors may be unreliable individually */
/*          (only the subspaces spanned are approximated reliably). */

/*  MINGAPI (output) INTEGER */
/*          Index I where the minimum gap in the spectrum occurred. */

/*  INFO    (output) INTEGER */
/*          = 0:  successful exit, no special cases occurred. */
/*          < -200: not enough workspace. Space for ABS(INFO + 200) */
/*                numbers is required in addition to the workspace provided, */
/*                otherwise some of the computed eigenvectors will be incorrect. */
/*          < -99, > -199: successful exit, but quick returns. */
/*                if INFO = -100, successful exit, but the input matrix */
/*                                was the zero matrix and no */
/*                                divide-and-conquer was performed */
/*                if INFO = -101, successful exit, but N was 1 and no */
/*                                divide-and-conquer was performed */
/*                if INFO = -102, successful exit, but only a single */
/*                                dense block. Standard dense solver */
/*                                was called, no divide-and-conquer was */
/*                                performed */
/*                if INFO = -103, successful exit, but warning that */
/*                                MINGAP.LE.TOL/10 and therefore the */
/*                                eigenvectors corresponding to close */
/*                                approximate eigenvalues may individually */
/*                                be unreliable (although taken together they */
/*                                do approximate the corresponding subspace to */
/*                                the desired accuracy) */
/*          = -99: error in the preprocessing in DIBTDC (when determining */
/*                 the merging order). */
/*          < 0, > -99: illegal arguments. */
/*                if INFO = -i, the i-th argument had an illegal value. */
/*          > 0:  The algorithm failed to compute an eigenvalue while */
/*                working on the submatrix lying in rows and columns */
/*                INFO/(N+1) through mod(INFO,N+1). */

/*  Further Details */
/*  =============== */

/*  Small modifications of code written by */
/*     Wilfried Gansterer and Bob Ward, */
/*     Department of Computer Science, University of Tennessee */
/*     see https://doi.org/10.1137/S1064827501399432 */

/*  Based on the design of the LAPACK code sstedc.f written by Jeff */
/*  Rutter, Computer Science Division, University of California at */
/*  Berkeley, and modified by Francoise Tisseur, University of Tennessee. */

/*  ===================================================================== */

/*     .. Parameters .. */

#define TOLMAX 0.1

/*        TOLMAX       .... upper bound for tolerances TOL, TAU1, TAU2 */
/*                          NOTE: in the routine DIBTDC, the value */
/*                                1.D-1 is hardcoded for TOLMAX ! */

  PetscBLASInt   i, j, k, i1, iwspc, lwmin, start;
  PetscBLASInt   ii, ip, nk, rk, np, iu, rp1, ldu;
  PetscBLASInt   ksk, ivt, iend, kchk=0, kmax=0, one=1, zero=0;
  PetscBLASInt   ldvt, ksum=0, kskp1, spneed, nrblks, liwmin, isvals;
  PetscReal      p, d2, eps, dmax, emax, done = 1.0;
  PetscReal      dnrm, tiny, anorm, exdnrm=0, dropsv, absdiff;

  PetscFunctionBegin;
  /* Determine machine epsilon. */
  eps = LAPACKlamch_("Epsilon");

  *oinfo = 0;

  if (*(unsigned char *)jobz != 'N' && *(unsigned char *)jobz != 'D') *oinfo = -1;
  else if (*(unsigned char *)jobacc != 'A' && *(unsigned char *)jobacc != 'M') *oinfo = -2;
  else if (n < 1) *oinfo = -3;
  else if (nblks < 1 || nblks > n) *oinfo = -4;
  if (*oinfo == 0) {
    for (k = 0; k < nblks; ++k) {
      ksk = ksizes[k];
      ksum += ksk;
      if (ksk > kmax) kmax = ksk;
      if (ksk < 1) kchk = 1;
    }
    if (nblks == 1) lwmin = 2*n*n + n*6 + 1;
    else lwmin = n*n + n*3;
    liwmin = n * 5 + nblks * 5 - 4;
    if (ksum != n || kchk == 1) *oinfo = -5;
    else if (l1d < PetscMax(3,kmax)) *oinfo = -7;
    else if (l2d < PetscMax(3,kmax)) *oinfo = -8;
    else if (l1e < PetscMax(3,2*kmax+1)) *oinfo = -10;
    else if (l2e < PetscMax(3,2*kmax+1)) *oinfo = -11;
    else if (*(unsigned char *)jobacc == 'A' && tol > TOLMAX) *oinfo = -12;
    else if (*(unsigned char *)jobacc == 'M' && tau1 > TOLMAX/2) *oinfo = -13;
    else if (*(unsigned char *)jobacc == 'M' && tau2 > TOLMAX/2) *oinfo = -14;
    else if (ldz < PetscMax(1,n)) *oinfo = -17;
    else if (lwork < lwmin) *oinfo = -19;
    else if (liwork < liwmin) *oinfo = -21;
  }

  PetscCheck(!*oinfo,PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"Wrong argument %" PetscBLASInt_FMT " in DSBTDC",-(*oinfo));

  /* Quick return if possible */

  if (n == 1) {
    ev[0] = d[0]; z[0] = 1.;
    *oinfo = -101;
    PetscFunctionReturn(PETSC_SUCCESS);
  }

  /* If NBLKS is equal to 1, then solve the problem with standard */
  /* dense solver (in this case KSIZES(1) = N). */

  if (nblks == 1) {
    for (i = 0; i < n; ++i) {
      for (j = 0; j <= i; ++j) {
        z[i + j*ldz] = d[i + j*l1d];
      }
    }
    PetscCallLAPACKInfo("LAPACKsyevd",LAPACKsyevd_("V", "L", &n, z, &ldz, ev, work, &lwork, iwork, &liwork, &info));
    *oinfo = -102;
    PetscFunctionReturn(PETSC_SUCCESS);
  }

  /* determine the accuracy parameters (if requested) */

  if (*(unsigned char *)jobacc == 'A') {
    tau1 = tol / 2;
    if (tau1 < eps * 20) tau1 = eps;
    tau2 = tol / 2;
  }

  /* Initialize Z as the identity matrix */

  if (*(unsigned char *)jobz == 'D') {
    for (j=0;j<n;j++) for (i=0;i<n;i++) z[i+j*ldz] = 0.0;
    for (i=0;i<n;i++) z[i+i*ldz] = 1.0;
  }

  /* Determine the off-diagonal ranks, form and store the lower rank */
  /* approximations based on the tolerance parameters, the */
  /* RANK(K) largest singular values and the associated singular */
  /* vectors of each subdiagonal block. Also find the maximum norm of */
  /* the subdiagonal blocks (in EMAX). */

  /* Compute SVDs of the subdiagonal blocks.... */

  /* EMAX .... maximum norm of the off-diagonal blocks */

  emax = 0.;
  for (k = 0; k < nblks-1; ++k) {
    ksk = ksizes[k];
    kskp1 = ksizes[k+1];
    isvals = 0;

    /* Note that min(KSKP1,KSK).LE.N/2 (equality possible for */
    /* NBLKS=2), and therefore storing the singular values requires */
    /* at most N/2 entries of the *        array WORK. */

    iu = isvals + n / 2;
    ivt = isvals + n / 2;

    /* Call of DGESVD: The space for U is not referenced, since */
    /* JOBU='O' and therefore this portion of the array WORK */
    /* is not referenced for U. */

    ldu = kskp1;
    ldvt = PetscMin(kskp1,ksk);
    iwspc = ivt + n * n / 2;

    /* Note that the minimum workspace required for this call */
    /* of DGESVD is: N/2 for storing the singular values + N**2/2 for */
    /* storing V^T + 5*N/2 workspace =  N**2/2 + 3*N. */

    i1 = lwork - iwspc;
    PetscCallLAPACKInfo("LAPACKgesvd",LAPACKgesvd_("O", "S", &kskp1, &ksk,
            &e[k*l1e*l2e], &l1e, &work[isvals],
            &work[iu], &ldu, &work[ivt], &ldvt, &work[iwspc], &i1, &info));

    /* Note that after the return from DGESVD U is stored in */
    /* E(*,*,K), and V^{\top} is stored in WORK(IVT, IVT+1, ....) */

    /* determine the ranks RANK() for the approximations */

    rk = PetscMin(ksk,kskp1);
L8:
    dropsv = work[isvals - 1 + rk];

    if (dropsv * 2. <= tau1) {

      /* the error caused by dropping singular value RK is */
      /* small enough, try to reduce the rank by one more */

      if (--rk > 0) goto L8;
      else iwork[k] = 0;
    } else {

      /* the error caused by dropping singular value RK is */
      /* too large already, RK is the rank required to achieve the */
      /* desired accuracy */

      iwork[k] = rk;
    }

/* ************************************************************************** */

    /* Store the first RANK(K) terms of the SVD of the current */
    /* off-diagonal block. */
    /* NOTE that here it is required that L1E, L2E >= 2*KMAX+1 in order */
    /* to have enough space for storing singular vectors and values up */
    /* to the full SVD of an off-diagonal block !!!! */

    /* u1-u_RANK(K) is already contained in E(:,1:RANK(K),K) (as a */
    /* result of the call of DGESVD !), the sigma1-sigmaK are to be */
    /* stored in E(1:RANK(K),RANK(K)+1,K),  and v1-v_RANK(K) are to be */
    /* stored in E(:,RANK(K)+2:2*RANK(K)+1,K) */

    rp1 = iwork[k];
    for (j = 0; j < iwork[k]; ++j) {

      /* store sigma_J in E(J,RANK(K)+1,K) */

      e[j + (rp1 + k*l2e)* l1e] = work[isvals + j];

      /* update maximum norm of subdiagonal blocks */

      if (e[j + (rp1 + k*l2e)*l1e] > emax) {
        emax = e[j + (rp1 + k*l2e)*l1e];
      }

      /* store v_J in E(:,RANK(K)+1+J,K) */
      /* (note that WORK contains V^{\top} and therefore */
      /* we need to read rowwise !) */

      for (i = 1; i <= ksk; ++i) {
        e[i-1 + (rp1+j+1 + k*l2e)*l1e] = work[ivt+j + (i-1)*ldvt];
      }
    }

  }

  /* Compute the maximum norm of diagonal blocks and store the norm */
  /* of each diagonal block in E(RP1,RP1,K) (after the singular values); */
  /* store the norm of the last diagonal block in EXDNRM. */

  /* DMAX .... maximum one-norm of the diagonal blocks */

  dmax = 0.;
  for (k = 0; k < nblks; ++k) {
    rp1 = iwork[k];

    /* compute the one-norm of diagonal block K */

    dnrm = LAPACKlansy_("1", "L", &ksizes[k], &d[k*l1d*l2d], &l1d, work);
    if (k+1 == nblks) exdnrm = dnrm;
    else e[rp1 + (rp1 + k*l2e)*l1e] = dnrm;
    if (dnrm > dmax) dmax = dnrm;
  }

  /* Check for zero matrix. */

  if (emax == 0. && dmax == 0.) {
    for (i = 0; i < n; ++i) ev[i] = 0.;
    *oinfo = -100;
    PetscFunctionReturn(PETSC_SUCCESS);
  }

/* **************************************************************** */

  /* ....Identify irreducible parts of the block tridiagonal matrix */
  /* [while (START <= NBLKS)].... */

  start = 0;
  np = 0;
L10:
  if (start < nblks) {

    /* Let IEND be the number of the next subdiagonal block such that */
    /* its RANK is 0 or IEND = NBLKS if no such subdiagonal exists. */
    /* The matrix identified by the elements between the diagonal block START */
    /* and the diagonal block IEND constitutes an independent (irreducible) */
    /* sub-problem. */

    iend = start;

L20:
    if (iend < nblks) {
      rk = iwork[iend];

      /* NOTE: if RANK(IEND).EQ.0 then decoupling happens due to */
      /*       reduced accuracy requirements ! (because in this case */
      /*       we would not merge the corresponding two diagonal blocks) */

      /* NOTE: seems like any combination may potentially happen: */
      /*       (i) RANK = 0 but no decoupling due to small norm of */
      /*           off-diagonal block (corresponding diagonal blocks */
      /*           also have small norm) as well as */
      /*       (ii) RANK > 0 but decoupling due to small norm of */
      /*           off-diagonal block (corresponding diagonal blocks */
      /*           have very large norm) */
      /*       case (i) is ruled out by checking for RANK = 0 above */
      /*       (we decide to decouple all the time when the rank */
      /*       of an off-diagonal block is zero, independently of */
      /*       the norms of the corresponding diagonal blocks. */

      if (rk > 0) {

        /* check for decoupling due to small norm of off-diagonal block */
        /* (relative to the norms of the corresponding diagonal blocks) */

        if (iend == nblks-2) {
          d2 = PetscSqrtReal(exdnrm);
        } else {
          d2 = PetscSqrtReal(e[iwork[iend+1] + (iwork[iend+1] + (iend+1)*l2e)*l1e]);
        }

        /* this definition of TINY is analogous to the definition */
        /* in the tridiagonal divide&conquer (dstedc) */

        tiny = eps * PetscSqrtReal(e[iwork[iend] + (iwork[iend] + iend*l2e)*l1e])*d2;
        if (e[(iwork[iend] + iend*l2e)*l1e] > tiny) {

          /* no decoupling due to small norm of off-diagonal block */

          ++iend;
          goto L20;
        }
      }
    }

    /* ....(Sub) Problem determined: between diagonal blocks */
    /*     START and IEND. Compute its size and solve it.... */

    nrblks = iend - start + 1;
    if (nrblks == 1) {

      /* Isolated problem is a single diagonal block */

      nk = ksizes[start];

      /* copy this isolated block into Z */

      for (i = 0; i < nk; ++i) {
        ip = np + i + 1;
        for (j = 0; j <= i; ++j) z[ip + (np+j+1)*ldz] = d[i + (j + start*l2d)*l1d];
      }

      /* check whether there is enough workspace */

      spneed = 2*nk*nk + nk * 6 + 1;
      PetscCheck(spneed<=lwork,PETSC_COMM_SELF,PETSC_ERR_MEM,"dsbtdc: not enough workspace for DSYEVD, info = %" PetscBLASInt_FMT,lwork - 200 - spneed);

      PetscCallLAPACKInfo("LAPACKsyevd",LAPACKsyevd_("V", "L", &nk,
                    &z[np + np*ldz], &ldz, &ev[np],
                    work, &lwork, &iwork[nblks-1], &liwork, &info));
      start = iend + 1;
      np += nk;

      /* go to the next irreducible subproblem */

      goto L10;
    }

    /* ....Isolated problem consists of more than one diagonal block. */
    /*     Start the divide and conquer algorithm.... */

    /* Scale: Divide by the maximum of all norms of diagonal blocks */
    /*        and singular values of the subdiagonal blocks */

    /* ....determine maximum of the norms of all diagonal and subdiagonal */
    /*     blocks.... */

    if (iend == nblks-1) anorm = exdnrm;
    else anorm = e[iwork[iend] + (iwork[iend] + iend*l2e)*l1e];
    for (k = start; k < iend; ++k) {
      rp1 = iwork[k];

      /* norm of diagonal block */
      anorm = PetscMax(anorm,e[rp1 + (rp1 + k*l2e)*l1e]);

      /* singular value of subdiagonal block */
      anorm = PetscMax(anorm,e[(rp1 + k*l2e)*l1e]);
    }

    nk = 0;
    for (k = start; k < iend+1; ++k) {
      ksk = ksizes[k];
      nk += ksk;

      /* scale the diagonal block */
      PetscCallLAPACKInfo("LAPACKlascl",LAPACKlascl_("L", &zero, &zero,
                    &anorm, &done, &ksk, &ksk, &d[k*l2d*l1d], &l1d, &info));

      /* scale the (approximated) off-diagonal block by dividing its */
      /* singular values */

      if (k != iend) {

        /* the last subdiagonal block has index IEND-1 !!!! */
        for (i = 0; i < iwork[k]; ++i) {
          e[i + (iwork[k] + k*l2e)*l1e] /= anorm;
        }
      }
    }

    /* call the block-tridiagonal divide-and-conquer on the */
    /* irreducible subproblem which has been identified */

    PetscCall(BDC_dibtdc_(jobz, nk, nrblks, &ksizes[start], &d[start*l1d*l2d], l1d, l2d,
                &e[start*l2e*l1e], &iwork[start], l1e, l2e, tau2, &ev[np],
                &z[np + np*ldz], ldz, work, lwork, &iwork[nblks-1], liwork, oinfo, 1));
    PetscCheck(!*oinfo,PETSC_COMM_SELF,PETSC_ERR_LIB,"dsbtdc: Error in DIBTDC, oinfo = %" PetscBLASInt_FMT,*oinfo);

/* ************************************************************************** */

    /* Scale back the computed eigenvalues. */

    PetscCallLAPACKInfo("LAPACKlascl",LAPACKlascl_("G", &zero, &zero, &done,
            &anorm, &nk, &one, &ev[np], &nk, &info));

    start = iend + 1;
    np += nk;

    /* Go to the next irreducible subproblem. */

    goto L10;
  }

  /* ....If the problem split any number of times, then the eigenvalues */
  /* will not be properly ordered. Here we permute the eigenvalues */
  /* (and the associated eigenvectors) across the irreducible parts */
  /* into ascending order.... */

  /*  IF(NRBLKS.LT.NBLKS)THEN */

  /*    Use Selection Sort to minimize swaps of eigenvectors */

  for (ii = 1; ii < n; ++ii) {
    i = ii;
    k = i;
    p = ev[i];
    for (j = ii; j < n; ++j) {
      if (ev[j] < p) {
        k = j;
        p = ev[j];
      }
    }
    if (k != i) {
      ev[k] = ev[i];
      ev[i] = p;
      PetscCallBLAS("BLASswap",BLASswap_(&n, &z[i*ldz], &one, &z[k*ldz], &one));
    }
  }

  /* ...Compute MINGAP (minimum difference between neighboring eigenvalue */
  /*    approximations).............................................. */

  *mingap = ev[1] - ev[0];
  PetscCheck(*mingap>=0.,PETSC_COMM_SELF,PETSC_ERR_LIB,"dsbtdc: Eigenvalue approximations are not ordered properly. Approximation 1 is larger than approximation 2.");
  *mingapi = 1;
  for (i = 2; i < n; ++i) {
    absdiff = ev[i] - ev[i-1];
    PetscCheck(absdiff>=0.,PETSC_COMM_SELF,PETSC_ERR_LIB,"dsbtdc: Eigenvalue approximations are not ordered properly. Approximation %" PetscBLASInt_FMT " is larger than approximation %" PetscBLASInt_FMT ".",i,i+1);
    if (absdiff < *mingap) {
      *mingap = absdiff;
      *mingapi = i;
    }
  }

  /* check whether the minimum gap between eigenvalue approximations */
  /* may indicate severe inaccuracies in the eigenvector approximations */

  if (*mingap <= tol / 10) *oinfo = -103;
  PetscFunctionReturn(PETSC_SUCCESS);
}
