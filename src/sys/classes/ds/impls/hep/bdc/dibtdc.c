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

static PetscErrorCode cutlr_(PetscBLASInt start,PetscBLASInt n,PetscBLASInt blkct,
        PetscBLASInt *bsizes,PetscBLASInt *ranks,PetscBLASInt *cut,
        PetscBLASInt *lsum,PetscBLASInt *lblks,PetscBLASInt *info)
{
/*  -- Routine written in LAPACK Version 3.0 style -- */
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

/*  CUTLR computes the optimal cut in a sequence of BLKCT neighboring */
/*  blocks whose sizes are given by the array BSIZES. */
/*  The sum of all block sizes in the sequence considered is given by N. */
/*  The cut is optimal in the sense that the difference of the sizes of */
/*  the resulting two halves is minimum over all cuts with minimum ranks */
/*  between blocks of the sequence considered. */

/*  Arguments */
/*  ========= */

/*  START  (input) INTEGER */
/*         In the original array KSIZES of the calling routine DIBTDC, */
/*         the position where the sequence considered in this routine starts. */
/*         START >= 1. */

/*  N      (input) INTEGER */
/*         The sum of all the block sizes of the sequence to be cut = */
/*         = sum_{i=1}^{BLKCT} BSIZES(I). */
/*         N >= 3. */

/*  BLKCT  (input) INTEGER */
/*         The number of blocks in the sequence to be cut. */
/*         BLKCT >= 3. */

/*  BSIZES (input) INTEGER array, dimension (BLKCT) */
/*         The dimensions of the (quadratic) blocks of the sequence to be */
/*         cut. sum_{i=1}^{BLKCT} BSIZES(I) = N. */

/*  RANKS  (input) INTEGER array, dimension (BLKCT-1) */
/*         The ranks determining the approximations of the off-diagonal */
/*         blocks in the sequence considered. */

/*  CUT    (output) INTEGER */
/*         After the optimum cut has been determined, the position (in the */
/*         overall problem as worked on in DIBTDC !) of the last block in */
/*         the first half of the sequence to be cut. */
/*         START <= CUT <= START+BLKCT-2. */

/*  LSUM   (output) INTEGER */
/*         After the optimum cut has been determined, the sum of the */
/*         block sizes in the first half of the sequence to be cut. */
/*         LSUM < N. */

/*  LBLKS  (output) INTEGER */
/*         After the optimum cut has been determined, the number of the */
/*         blocks in the first half of the sequence to be cut. */
/*         1 <= LBLKS < BLKCT. */

/*  INFO   (output) INTEGER */
/*          = 0:  successful exit. */
/*          < 0:  illegal arguments. */
/*                if INFO = -i, the i-th (input) argument had an illegal */
/*                value. */
/*          > 0:  illegal results. */
/*                if INFO = i, the i-th (output) argument had an illegal */
/*                value. */

/*  Further Details */
/*  =============== */

/*  Based on code written by */
/*     Wilfried Gansterer and Bob Ward, */
/*     Department of Computer Science, University of Tennessee */

/*  ===================================================================== */

  PetscBLASInt i, ksk, kchk, ksum, nhalf, deviat, mindev, minrnk, tmpsum;

  PetscFunctionBegin;
  *info = 0;
  *lblks = 1;
  *lsum = 1;
  *cut = start;

  if (start < 1) *info = -1;
  else if (n < 3) *info = -2;
  else if (blkct < 3) *info = -3;
  if (*info == 0) {
    ksum = 0;
    kchk = 0;
    for (i = 0; i < blkct; ++i) {
      ksk = bsizes[i];
      ksum += ksk;
      if (ksk < 1) kchk = 1;
    }
    if (ksum != n || kchk == 1) *info = -4;
  }
  PetscCheck(!*info,PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"Wrong argument %" PetscBLASInt_FMT " in CUTLR",-(*info));

  /* determine smallest rank in the range considered */

  minrnk = n;
  for (i = 0; i < blkct-1; ++i) {
    if (ranks[i] < minrnk) minrnk = ranks[i];
  }

  /* determine best cut among those with smallest rank */

  nhalf = n / 2;
  tmpsum = 0;
  mindev = n;
  for (i = 0; i < blkct; ++i) {
    tmpsum += bsizes[i];
    if (ranks[i] == minrnk) {

      /* determine deviation from "optimal" cut NHALF */

      deviat = tmpsum - nhalf;
      if (deviat<0) deviat = -deviat;

      /* compare to best deviation so far */

      if (deviat < mindev) {
        mindev = deviat;
        *cut = start + i;
        *lblks = i + 1;
        *lsum = tmpsum;
      }
    }
  }

  if (*cut < start || *cut >= start + blkct - 1) *info = 6;
  else if (*lsum < 1 || *lsum >= n) *info = 7;
  else if (*lblks < 1 || *lblks >= blkct) *info = 8;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode BDC_dibtdc_(const char *jobz,PetscBLASInt n,PetscBLASInt nblks,
        PetscBLASInt *ksizes,PetscReal *d,PetscBLASInt l1d,PetscBLASInt l2d,
        PetscReal *e,PetscBLASInt *rank,PetscBLASInt l1e,PetscBLASInt l2e,
        PetscReal tol,PetscReal *ev,PetscReal *z,PetscBLASInt ldz,PetscReal *work,
        PetscBLASInt lwork,PetscBLASInt *iwork,PetscBLASInt liwork,
        PetscBLASInt *oinfo,PetscBLASInt jobz_len)
{
/*  -- Routine written in LAPACK Version 3.0 style -- */
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

/*  DIBTDC computes all eigenvalues and corresponding eigenvectors of a */
/*  symmetric irreducible block tridiagonal matrix with rank RANK matrices */
/*  as the subdiagonal blocks using a block divide and conquer method. */

/*  Arguments */
/*  ========= */

/*  JOBZ    (input) CHARACTER*1 */
/*          = 'N':  Compute eigenvalues only (not implemented); */
/*          = 'D':  Compute eigenvalues and eigenvectors. */
/*                  Eigenvectors are accumulated in the */
/*                  divide-and-conquer process. */

/*  N      (input) INTEGER */
/*         The dimension of the symmetric irreducible block tridiagonal */
/*         matrix.  N >= 2. */

/*  NBLKS  (input) INTEGER, 2 <= NBLKS <= N */
/*         The number of diagonal blocks in the matrix. */

/*  KSIZES (input) INTEGER array, dimension (NBLKS) */
/*         The dimension of the square diagonal blocks from top left */
/*         to bottom right.  KSIZES(I) >= 1 for all I, and the sum of */
/*         KSIZES(I) for I = 1 to NBLKS has to be equal to N. */

/*  D      (input) DOUBLE PRECISION array, dimension (L1D,L2D,NBLKS) */
/*         The lower triangular elements of the symmetric diagonal */
/*         blocks of the block tridiagonal matrix.  Elements of the top */
/*         left diagonal block, which is of dimension KSIZES(1), are */
/*         contained in D(*,*,1); the elements of the next diagonal */
/*         block, which is of dimension KSIZES(2), are contained in */
/*         D(*,*,2); etc. */

/*  L1D    (input) INTEGER */
/*         The leading dimension of the array D.  L1D >= max(3,KMAX), */
/*         where KMAX is the dimension of the largest diagonal block. */

/*  L2D    (input) INTEGER */
/*         The second dimension of the array D.  L2D >= max(3,KMAX), */
/*         where KMAX is as stated in L1D above. */

/*  E      (input) DOUBLE PRECISION array, dimension (L1E,L2E,NBLKS-1) */
/*         Contains the elements of the scalars (singular values) and */
/*         vectors (singular vectors) defining the rank RANK subdiagonal */
/*         blocks of the matrix. */
/*         E(1:RANK(K),RANK(K)+1,K) holds the RANK(K) scalars, */
/*         E(:,1:RANK(K),K) holds the RANK(K) column vectors, and */
/*         E(:,RANK(K)+2:2*RANK(K)+1,K) holds the row vectors for the K-th */
/*         subdiagonal block. */

/*  RANK   (input) INTEGER array, dimension (NBLKS-1). */
/*         The ranks of all the subdiagonal blocks contained in the array E. */
/*         RANK(K) <= MIN(KSIZES(K), KSIZES(K+1)) */

/*  L1E    (input) INTEGER */
/*         The leading dimension of the array E.  L1E >= max(3,2*KMAX+1), */
/*         where KMAX is as stated in L1D above. */

/*  L2E    (input) INTEGER */
/*         The second dimension of the array E.  L2E >= max(3,2*KMAX+1), */
/*         where KMAX is as stated in L1D above. */

/*  TOL    (input) DOUBLE PRECISION, TOL <= 1.0D-1 */
/*         User specified deflation tolerance for the routine DMERG2. */
/*         If (1.0D-1 >= TOL >= 20*EPS) then TOL is used as */
/*         the deflation tolerance in DSRTDF. */
/*         If (TOL < 20*EPS) then the standard deflation tolerance from */
/*         LAPACK is used as the deflation tolerance in DSRTDF. */

/*  EV     (output) DOUBLE PRECISION array, dimension (N) */
/*         If INFO = 0, then EV contains the eigenvalues of the */
/*         symmetric block tridiagonal matrix in ascending order. */

/*  Z      (input/output) DOUBLE PRECISION array, dimension (LDZ, N) */
/*         On entry, Z will be the identity matrix. */
/*         On exit, Z contains the eigenvectors of the block tridiagonal */
/*         matrix. */

/*  LDZ    (input) INTEGER */
/*         The leading dimension of the array Z.  LDZ >= max(1,N). */

/*  WORK   (workspace) DOUBLE PRECISION array, dimension (LWORK) */

/*  LWORK   (input) INTEGER */
/*          The dimension of the array WORK. */
/*          In order to guarantee correct results in all cases, */
/*          LWORK must be at least (2*N**2 + 3*N). In many cases, */
/*          less workspace is required. The absolute minimum required is */
/*          (N**2 + 3*N). */
/*          If the workspace provided is not sufficient, the routine will */
/*          return a corresponding error code and report how much workspace */
/*          was missing (see INFO). */

/*  IWORK  (workspace) INTEGER array, dimension (LIWORK) */

/*  LIWORK  (input) INTEGER */
/*          The dimension of the array IWORK. */
/*          LIWORK must be at least (5*N + 3 + 4*NBLKS - 4): */
/*                 5*KMAX+3 for DSYEVD, 5*N for ????, */
/*                 4*NBLKS-4 for the preprocessing (merging order) */
/*          Summarizing, the minimum integer workspace needed is */
/*          MAX(5*N, 5*KMAX + 3) + 4*NBLKS - 4 */

/*  INFO   (output) INTEGER */
/*          = 0:  successful exit. */
/*          < 0, > -99: illegal arguments. */
/*                if INFO = -i, the i-th argument had an illegal value. */
/*          = -99: error in the preprocessing (call of CUTLR). */
/*          < -200: not enough workspace. Space for ABS(INFO + 200) */
/*                numbers is required in addition to the workspace provided, */
/*                otherwise some eigenvectors will be incorrect. */
/*          > 0:  The algorithm failed to compute an eigenvalue while */
/*                working on the submatrix lying in rows and columns */
/*                INFO/(N+1) through mod(INFO,N+1). */

/*  Further Details */
/*  =============== */

/*  Based on code written by */
/*     Wilfried Gansterer and Bob Ward, */
/*     Department of Computer Science, University of Tennessee */

/*  This routine is comparable to Dlaed0.f from LAPACK. */

/*  ===================================================================== */

  PetscBLASInt   i, j, k, np, rp1, ksk, one=1;
  PetscBLASInt   cut, mat1, kchk, kbrk, blks, kmax, icut, size, ksum, lsum;
  PetscBLASInt   lblks, rblks, isize, lwmin, ilsum;
  PetscBLASInt   start, istck1, istck2, istck3, merged;
  PetscBLASInt   liwmin, matsiz, startp, istrtp;
  PetscReal      rho, done=1.0, dmone=-1.0;

  PetscFunctionBegin;
  *oinfo = 0;

  if (*(unsigned char *)jobz != 'N' && *(unsigned char *)jobz != 'D') *oinfo = -1;
  else if (n < 2) *oinfo = -2;
  else if (nblks < 2 || nblks > n) *oinfo = -3;
  if (*oinfo == 0) {
    ksum = 0;
    kmax = 0;
    kchk = 0;
    for (k = 0; k < nblks; ++k) {
      ksk = ksizes[k];
      ksum += ksk;
      if (ksk > kmax) kmax = ksk;
      if (ksk < 1) kchk = 1;
    }
    lwmin = n*n + n * 3;
    liwmin = PetscMax(n * 5,kmax * 5 + 3) + 4*nblks - 4;
    if (ksum != n || kchk == 1) *oinfo = -4;
    else if (l1d < PetscMax(3,kmax)) *oinfo = -6;
    else if (l2d < PetscMax(3,kmax)) *oinfo = -7;
    else if (l1e < PetscMax(3,2*kmax + 1)) *oinfo = -10;
    else if (l2e < PetscMax(3,2*kmax + 1)) *oinfo = -11;
    else if (tol > .1) *oinfo = -12;
    else if (ldz < PetscMax(1,n)) *oinfo = -15;
    else if (lwork < lwmin) *oinfo = -17;
    else if (liwork < liwmin) *oinfo = -19;
  }
  if (*oinfo == 0) {
    for (k = 0; k < nblks-1; ++k) {
      if (rank[k] > PetscMin(ksizes[k],ksizes[k+1]) || rank[k] < 1) *oinfo = -9;
    }
  }

  PetscCheck(!*oinfo,PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"Wrong argument %" PetscBLASInt_FMT " in DIBTDC",-(*oinfo));

/* **************************************************************************** */

  /* ...Preprocessing..................................................... */
  /*    Determine the optimal order for merging the subblocks and how much */
  /*    workspace will be needed for the merging (determined by the last */
  /*    merge). Cutpoints for the merging operations are determined and stored */
  /*    in reverse chronological order (starting with the final merging */
  /*    operation). */

  /*    integer workspace requirements for the preprocessing: */
  /*         4*(NBLKS-1) for merging history */
  /*         at most 3*(NBLKS-1) for stack */

  start = 1;
  size = n;
  blks = nblks;
  merged = 0;
  k = 0;

  /* integer workspace used for the stack is not needed any more after the */
  /* preprocessing and therefore can use part of the 5*N */
  /* integer workspace needed later on in the code */

  istck1 = 0;
  istck2 = istck1 + nblks;
  istck3 = istck2 + nblks;

  /* integer workspace used for storing the order of merges starts AFTER */
  /* the integer workspace 5*N+3 which is needed later on in the code */
  /* (5*KMAX+3 for DSYEVD, 4*N in DMERG2) */

  istrtp = n * 5 + 4;
  icut = istrtp + nblks - 1;
  isize = icut + nblks - 1;
  ilsum = isize + nblks - 1;

L200:

  if (nblks >= 3) {

    /* Determine the cut point. Note that in the routine CUTLR it is */
    /* chosen such that it yields the best balanced merging operation */
    /* among all the rank modifications with minimum rank. */

    PetscCall(cutlr_(start, size, blks, &ksizes[start-1], &rank[start-1], &cut, &lsum, &lblks, oinfo));
    PetscCheck(!*oinfo,PETSC_COMM_SELF,PETSC_ERR_PLIB,"dibtdc: Error in cutlr, oinfo = %" PetscBLASInt_FMT,*oinfo);

  } else {
    cut = 1;
    lsum = ksizes[0];
    lblks = 1;
  }

  ++merged;
  startp = 0;
  for (i = 0; i < start-1; ++i) startp += ksizes[i];
  iwork[istrtp + (nblks - 1) - merged-1] = startp + 1;
  iwork[icut + (nblks - 1) - merged-1] = cut;
  iwork[isize + (nblks - 1) - merged-1] = size;
  iwork[ilsum + (nblks - 1) - merged-1] = lsum;

  if (lblks == 2) {

    /* one merge in left branch, left branch done */
    ++merged;
    iwork[istrtp + (nblks - 1) - merged-1] = startp + 1;
    iwork[icut + (nblks - 1) - merged-1] = start;
    iwork[isize + (nblks - 1) - merged-1] = lsum;
    iwork[ilsum + (nblks - 1) - merged-1] = ksizes[start-1];
  }

  if (lblks == 1 || lblks == 2) {

    /* left branch done, continue on the right side */
    start += lblks;
    size -= lsum;
    blks -= lblks;

    PetscCheck(blks>0,PETSC_COMM_SELF,PETSC_ERR_PLIB,"dibtdc: Error in preprocessing, blks = %" PetscBLASInt_FMT,blks);

    if (blks == 2) {

      /* one merge in right branch, right branch done */
      ++merged;
      startp += lsum;
      iwork[istrtp + (nblks - 1) - merged-1] = startp + 1;
      iwork[icut + (nblks - 1) - merged-1] = start;
      iwork[isize + (nblks - 1) - merged-1] = size;
      iwork[ilsum + (nblks - 1) - merged-1] = ksizes[start-1];
    }

    if (blks == 1 || blks == 2) {

      /* get the next subproblem from the stack or finished */

      if (k >= 1) {

        /* something left on the stack */
        start = iwork[istck1 + k-1];
        size = iwork[istck2 + k-1];
        blks = iwork[istck3 + k-1];
        --k;
        goto L200;
      } else {

        /* nothing left on the stack */
        PetscCheck(merged==nblks-1,PETSC_COMM_SELF,PETSC_ERR_PLIB,"ERROR in preprocessing - not enough merges performed");

        /* exit preprocessing */

      }
    } else {

      /* BLKS.GE.3, and therefore analyze the right side */

      goto L200;
    }
  } else {

    /* LBLKS.GE.3, and therefore check the right side and */
    /* put it on the stack if required */

    rblks = blks - lblks;
    if (rblks >= 3) {
      ++k;
      iwork[istck1 + k-1] = cut + 1;
      iwork[istck2 + k-1] = size - lsum;
      iwork[istck3 + k-1] = rblks;
    } else if (rblks == 2) {

      /* one merge in right branch, right branch done */
      /* (note that nothing needs to be done if RBLKS.EQ.1 !) */

      ++merged;
      startp += lsum;
      iwork[istrtp + (nblks - 1) - merged-1] = startp + 1;
      iwork[icut + (nblks - 1) - merged-1] = start + lblks;
      iwork[isize + (nblks - 1) - merged-1] = size - lsum;
      iwork[ilsum + (nblks - 1) - merged-1] = ksizes[start + lblks-1];
    }
    PetscCheck(rblks>0,PETSC_COMM_SELF,PETSC_ERR_PLIB,"dibtdc: ERROR in preprocessing - rblks = %" PetscBLASInt_FMT,rblks);

    /* continue on the left side */

    size = lsum;
    blks = lblks;
    goto L200;
  }

  /*  SIZE = IWORK(ISIZE+NBLKS-2) */
  /*  MAT1 = IWORK(ILSUM+NBLKS-2) */

  /* Note: after the dimensions SIZE and MAT1 of the last merging */
  /* operation have been determined, an upper bound for the workspace */
  /* requirements which is independent of how much deflation occurs in */
  /* the last merging operation could be determined as follows */
  /* (based on (3.15) and (3.19) from UT-CS-00-447): */

  /*  IF(MAT1.LE.N/2) THEN */
  /*     WSPREQ = 3*N + 3/2*(SIZE-MAT1)**2 + N*N/2 + MAT1*MAT1 */
  /*  ELSE */
  /*     WSPREQ = 3*N + 3/2*MAT1*MAT1 + N*N/2 + (SIZE-MAT1)**2 */
  /*  END IF */

  /*  IF(LWORK-WSPREQ.LT.0)THEN */
  /*          not enough work space provided */
  /*     INFO = -200 - (WSPREQ-LWORK) */
  /*     RETURN */
  /*  END IF */
  /*  However, this is not really useful, since the actual check whether */
  /*  enough workspace is provided happens in DMERG2.f ! */

/* ************************************************************************* */

  /* ...Solve subproblems................................... */

  /* Divide the matrix into NBLKS submatrices using rank-r */
  /* modifications (cuts) and solve for their eigenvalues and */
  /* eigenvectors. Initialize index array to sort eigenvalues. */

  /* first block: ...................................... */

  /*    correction for block 1: D1 - V1 \Sigma1 V1^T */

  ksk = ksizes[0];
  rp1 = rank[0];

  /* initialize the proper part of Z with the diagonal block D1 */
  /* (the correction will be made in Z and then the call of DSYEVD will */
  /*  overwrite it with the eigenvectors) */

  for (j=0;j<ksk;j++) for (i=j;i<ksk;i++) z[i+j*ldz] = d[i+j*l1d];

  /* copy D1 into WORK (in order to be able to restore it afterwards) */

  for (j=0;j<ksk;j++) for (i=j;i<ksk;i++) work[i+j*ksk] = d[i+j*l1d];

  /* copy V1 into the first RANK(1) columns of D1 and then */
  /* multiply with \Sigma1 */

  for (i = 0; i < rank[0]; ++i) {
    PetscCallBLAS("BLAScopy",BLAScopy_(&ksk, &e[(rp1 + i+1)*l1e], &one, &d[i*l1d], &one));
    PetscCallBLAS("BLASscal",BLASscal_(&ksk, &e[i + rp1*l1e], &d[i*l1d], &one));
  }

  /* multiply the first RANK(1) columns of D1 with V1^T and */
  /* subtract the result from the proper part of Z (previously */
  /* initialized with D1) */

  PetscCallBLAS("BLASgemm",BLASgemm_("N", "T", &ksk, &ksk, rank, &dmone,
          d, &l1d, &e[(rank[0]+1)*l1e], &l1e, &done, z, &ldz));

  /* restore the original D1 from WORK */

  for (j=0;j<ksk;j++) for (i=j;i<ksk;i++) d[i+j*l1d] = work[i+j*ksk];

  /* eigenanalysis of block 1 (using DSYEVD) */

  PetscCallLAPACKInfo("LAPACKsyev",LAPACKsyev_("V", "L", &ksk, z, &ldz, ev, work, &lwork, &info));

  /* EV(1:) contains the eigenvalues in ascending order */
  /* (they are returned this way by DSYEVD) */

  for (i = 0; i < ksk; ++i) iwork[i] = i+1;

    /* intermediate blocks: .............................. */

    np = ksk;

    /* remaining number of blocks */

    if (nblks > 2) {
      for (k = 1; k < nblks-1; ++k) {

      /* correction for block K: */
      /* Dk - U(k-1) \Sigma(k-1) U(k-1)^T - Vk \Sigmak Vk^T */

      ksk = ksizes[k];
      rp1 = rank[k];

      /* initialize the proper part of Z with the diagonal block Dk */
      /* (the correction will be made in Z and then the call of DSYEVD will */
      /*  overwrite it with the eigenvectors) */

      for (j=0;j<ksk;j++) for (i=j;i<ksk;i++) z[np+i+(np+j)*ldz] = d[i+(j+k*l2d)*l1d];

      /* copy Dk into WORK (in order to be able to restore it afterwards) */

      for (j=0;j<ksk;j++) for (i=j;i<ksk;i++) work[i+j*ksk] = d[i+(j+k*l2d)*l1d];

      /* copy U(K-1) into the first RANK(K-1) columns of Dk and then */
      /* multiply with \Sigma(K-1) */

      for (i = 0; i < rank[k-1]; ++i) {
        PetscCallBLAS("BLAScopy",BLAScopy_(&ksk, &e[(i+(k-1)*l2e)*l1e], &one, &d[(i+k*l2d)*l1d], &one));
        PetscCallBLAS("BLASscal",BLASscal_(&ksk, &e[i+(rank[k-1]+(k-1)*l2e)*l1e], &d[(i+k*l2d)*l1d], &one));
      }

      /* multiply the first RANK(K-1) columns of Dk with U(k-1)^T and */
      /* subtract the result from the proper part of Z (previously */
      /* initialized with Dk) */

      PetscCallBLAS("BLASgemm",BLASgemm_("N", "T", &ksk, &ksk, &rank[k-1],
                    &dmone, &d[k*l1d*l2d],
                    &l1d, &e[(k-1)*l1e*l2e], &l1e, &done, &z[np+np*ldz], &ldz));

      /* copy Vk into the first RANK(K) columns of Dk and then */
      /* multiply with \Sigmak */

      for (i = 0; i < rank[k]; ++i) {
        PetscCallBLAS("BLAScopy",BLAScopy_(&ksk, &e[(rp1+i+1 + k*l2e)*l1e], &one, &d[(i + k*l2d)*l1d], &one));
        PetscCallBLAS("BLASscal",BLASscal_(&ksk, &e[i + (rp1 + k*l2e)*l1e], &d[(i + k*l2d)*l1d], &one));
      }

      /* multiply the first RANK(K) columns of Dk with Vk^T and */
      /* subtract the result from the proper part of Z (previously */
      /* updated with [- U(k-1) \Sigma(k-1) U(k-1)^T]) */

      PetscCallBLAS("BLASgemm",BLASgemm_("N", "T", &ksk, &ksk, &rank[k],
                    &dmone, &d[k*l1d*l2d], &l1d,
                    &e[(rank[k]+1 + k*l2e)*l1e], &l1e, &done, &z[np+np*ldz], &ldz));

      /* restore the original Dk from WORK */

      for (j=0;j<ksk;j++) for (i=j;i<ksk;i++) d[i+(j+k*l2d)*l1d] = work[i+j*ksk];

      /* eigenanalysis of block K (using dsyevd) */

      PetscCallLAPACKInfo("LAPACKsyev",LAPACKsyev_("V", "L", &ksk, &z[np+np*ldz],
                     &ldz, &ev[np], work, &lwork, &info));

      /* EV(NPP1:) contains the eigenvalues in ascending order */
      /* (they are returned this way by DSYEVD) */

      for (i = 0; i < ksk; ++i) iwork[np + i] = i+1;

      /* update NP */
      np += ksk;
    }
  }

  /* last block: ....................................... */

  /*    correction for block NBLKS: */
  /*    D(nblks) - U(nblks-1) \Sigma(nblks-1) U(nblks-1)^T */

  ksk = ksizes[nblks-1];

  /* initialize the proper part of Z with the diagonal block D(nblks) */
  /* (the correction will be made in Z and then the call of DSYEVD will */
  /* overwrite it with the eigenvectors) */

  for (j=0;j<ksk;j++) for (i=j;i<ksk;i++) z[np+i+(np+j)*ldz] = d[i+(j+(nblks-1)*l2d)*l1d];

  /* copy D(nblks) into WORK (in order to be able to restore it afterwards) */

  for (j=0;j<ksk;j++) for (i=j;i<ksk;i++) work[i+j*ksk] = d[i+(j+(nblks-1)*l2d)*l1d];

  /* copy U(nblks-1) into the first RANK(nblks-1) columns of D(nblks) and then */
  /* multiply with \Sigma(nblks-1) */

  for (i = 0; i < rank[nblks-2]; ++i) {
    PetscCallBLAS("BLAScopy",BLAScopy_(&ksk, &e[(i + (nblks-2)*l2e)*l1e],
              &one, &d[(i + (nblks-1)*l2d)*l1d], &one));
    PetscCallBLAS("BLASscal",BLASscal_(&ksk,
              &e[i + (rank[nblks-2] + (nblks-2)*l2e)*l1e],
              &d[(i + (nblks-1)*l2d)*l1d], &one));
  }

  /* multiply the first RANK(nblks-1) columns of D(nblks) with U(nblks-1)^T */
  /* and subtract the result from the proper part of Z (previously */
  /* initialized with D(nblks)) */

  PetscCallBLAS("BLASgemm",BLASgemm_("N", "T", &ksk, &ksk, &rank[nblks - 2],
          &dmone, &d[(nblks-1)*l1d*l2d], &l1d,
          &e[(nblks-2)*l1e*l2e], &l1e, &done, &z[np+np*ldz], &ldz));

  /* restore the original D(nblks) from WORK */

  for (j=0;j<ksk;j++) for (i=j;i<ksk;i++) d[i+(j+(nblks-1)*l2d)*l1d] = work[i+j*ksk];

  /* eigenanalysis of block NBLKS (using dsyevd) */

  PetscCallLAPACKInfo("LAPACKsyev",LAPACKsyev_("V", "L", &ksk, &z[np+np*ldz], &ldz, &ev[np], work, &lwork, &info));

  /* EV(NPP1:) contains the eigenvalues in ascending order */
  /* (they are returned this way by DSYEVD) */

  for (i = 0; i < ksk; ++i) iwork[np + i] = i+1;

  /* note that from here on the entire workspace is available again */

  /* Perform all the merging operations. */

  for (i = 0; i < nblks-1; ++i) {

    /* MATSIZ = total size of the current rank RANK modification problem */

    matsiz = iwork[isize + i - 1];
    np = iwork[istrtp + i - 1];
    kbrk = iwork[icut + i - 1];
    mat1 = iwork[ilsum + i - 1];

    for (j = 0; j < rank[kbrk-1]; ++j) {

      /* NOTE: The parameter RHO in DMERG2 is modified in DSRTDF */
      /*       (multiplied by 2) ! In order not to change the */
      /*       singular value stored in E(:, RANK(KBRK)+1, KBRK), */
      /*       we do not pass on this variable as an argument to DMERG2, */
      /*       but we assign a separate variable RHO here which is passed */
      /*       on to DMERG2. */
      /*       Alternative solution in F90: */
      /*       pass E(:,RANK(KBRK)+1,KBRK) to an INTENT(IN) parameter */
      /*       in DMERG2. */

      rho = e[j + (rank[kbrk-1] + (kbrk-1)*l2e)*l1e];

      /* eigenvectors are accumulated (JOBZ.EQ.'D') */

      PetscCall(BDC_dmerg2_(jobz, j+1, matsiz, &ev[np-1], &z[np-1+(np-1)*ldz],
                    ldz, &iwork[np-1], &rho, &e[(j + (kbrk-1)*l2e)*l1e],
                    ksizes[kbrk], &e[(rank[kbrk-1]+j+1 + (kbrk-1)*l2e)*l1e],
                    ksizes[kbrk-1], mat1, work, lwork, &iwork[n], tol, oinfo, 1));
      PetscCheck(!*oinfo,PETSC_COMM_SELF,PETSC_ERR_PLIB,"dibtdc: Error in dmerg2, oinfo = %" PetscBLASInt_FMT,*oinfo);
    }

    /* at this point all RANK(KBRK) rank-one modifications corresponding */
    /* to the current off-diagonal block are finished. */
    /* Move on to the next off-diagonal block. */

  }

  /* Re-merge the eigenvalues/vectors which were deflated at the final */
  /* merging step by sorting all eigenvalues and eigenvectors according */
  /* to the permutation stored in IWORK. */

  /* copy eigenvalues and eigenvectors in ordered form into WORK */
  /* (eigenvalues into WORK(1:N), eigenvectors into WORK(N+1:N+1+N^2)) */

  for (i = 0; i < n; ++i) {
    j = iwork[i];
    work[i] = ev[j-1];
    PetscCallBLAS("BLAScopy",BLAScopy_(&n, &z[(j-1)*ldz], &one, &work[n*(i+1)], &one));
  }

  /* copy ordered eigenvalues back from WORK(1:N) into EV */

  PetscCallBLAS("BLAScopy",BLAScopy_(&n, work, &one, ev, &one));

  /* copy ordered eigenvectors back from WORK(N+1:N+1+N^2) into Z */

  for (j=0;j<n;j++) for (i=0;i<n;i++) z[i+j*ldz] = work[i+(j+1)*n];
  PetscFunctionReturn(PETSC_SUCCESS);
}
