/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Logarithm function  log(x)
*/

#include <slepc/private/fnimpl.h>      /*I "slepcfn.h" I*/
#include <slepcblaslapack.h>

static PetscErrorCode FNEvaluateFunction_Log(FN fn,PetscScalar x,PetscScalar *y)
{
  PetscFunctionBegin;
#if !defined(PETSC_USE_COMPLEX)
  PetscCheck(x>=0.0,PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Function not defined in the requested value");
#endif
  *y = PetscLogScalar(x);
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode FNEvaluateDerivative_Log(FN fn,PetscScalar x,PetscScalar *y)
{
  PetscFunctionBegin;
  PetscCheck(x!=0.0,PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Derivative not defined in the requested value");
#if !defined(PETSC_USE_COMPLEX)
  PetscCheck(x>0.0,PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Derivative not defined in the requested value");
#endif
  *y = 1.0/x;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   Block structure of a quasitriangular matrix T. Returns a list of n-1 numbers, where
   structure(j) encodes the block type of the j:j+1,j:j+1 diagonal block as one of:
      0 - not the start of a block
      1 - start of a 2x2 triangular block
      2 - start of a 2x2 quasi-triangular (full) block
*/
static PetscErrorCode qtri_struct(PetscBLASInt n,PetscScalar *T,PetscBLASInt ld,PetscInt *structure)
{
  PetscBLASInt j;

  PetscFunctionBegin;
#if defined(PETSC_USE_COMPLEX)
  for (j=0;j<n-1;j++) structure[j] = 1;
#else
  if (n==1) PetscFunctionReturn(PETSC_SUCCESS);
  else if (n==2) {
    structure[0] = (T[1]==0.0)? 1: 2;
    PetscFunctionReturn(PETSC_SUCCESS);
  }
  j = 0;
  while (j<n-2) {
    if (T[j+1+j*ld] != 0.0) { /* Start of a 2x2 full block */
      structure[j++] = 2;
      structure[j++] = 0;
    } else if (T[j+1+j*ld] == 0.0 && T[j+2+(j+1)*ld] == 0.0) { /* Start of a 2x2 triangular block */
      structure[j++] = 1;
    } else { /* Next block must start a 2x2 full block */
      structure[j++] = 0;
    }
  }
  if (T[n-1+(n-2)*ld] != 0.0) { /* 2x2 full block at the end */
    structure[n-2] = 2;
  } else if (structure[n-3] == 0 || structure[n-3] == 1) {
    structure[n-2] = 1;
  }
#endif
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   Compute scaling parameter (s) and order of Pade approximant (m).
   wr,wi is overwritten. Required workspace is 3*n*n.
   On output, Troot contains the sth square root of T.
*/
static PetscErrorCode FNlogm_params(FN fn,PetscBLASInt n,PetscScalar *T,PetscBLASInt ld,PetscScalar *wr,PetscScalar *wi,PetscInt maxroots,PetscInt *s,PetscInt *m,PetscScalar *Troot,PetscScalar *work)
{
  PetscInt        i,j,k,p,s0;
  PetscReal       inrm,eta,a2,a3,a4,d2,d3,d4,d5;
  PetscScalar     *TrootmI=work+2*n*n;
  PetscBool       foundm=PETSC_FALSE,more;
  PetscRandom     rand;
  const PetscReal xvals[] = { 1.586970738772063e-005, 2.313807884242979e-003, 1.938179313533253e-002,
       6.209171588994762e-002, 1.276404810806775e-001, 2.060962623452836e-001, 2.879093714241194e-001 };
  const PetscInt  mmax=PETSC_STATIC_ARRAY_LENGTH(xvals);

  PetscFunctionBegin;
  PetscCall(PetscRandomCreate(PETSC_COMM_SELF,&rand));
  /* get initial s0 so that T^(1/2^s0) < xvals(mmax) */
  *s = 0;
  do {
    inrm = SlepcAbsEigenvalue(wr[0]-PetscRealConstant(1.0),wi[0]);
    for (i=1;i<n;i++) inrm = PetscMax(inrm,SlepcAbsEigenvalue(wr[i]-PetscRealConstant(1.0),wi[i]));
    if (inrm < xvals[mmax-1]) break;
    for (i=0;i<n;i++) {
#if defined(PETSC_USE_COMPLEX)
      wr[i] = PetscSqrtScalar(wr[i]);
#else
#if defined(PETSC_HAVE_COMPLEX)
      PetscComplex z = PetscSqrtComplex(PetscCMPLX(wr[i],wi[i]));
      wr[i] = PetscRealPartComplex(z);
      wi[i] = PetscImaginaryPartComplex(z);
#else
      SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"This operation requires a compiler with C99 or C++ complex support");
#endif
#endif
    }
    *s = *s + 1;
  } while (*s<maxroots);
  s0 = *s;
  if (*s == maxroots) PetscCall(PetscInfo(fn,"Too many matrix square roots\n"));

  /* Troot = T */
  for (j=0;j<n;j++) PetscCall(PetscArraycpy(Troot+j*ld,T+j*ld,PetscMin(j+2,n)));
  for (k=1;k<=PetscMin(*s,maxroots);k++) PetscCall(FNSqrtmSchur(fn,n,Troot,ld,PETSC_FALSE));
  /* Compute value of s and m needed */
  /* TrootmI = Troot - I */
  for (j=0;j<n;j++) {
    PetscCall(PetscArraycpy(TrootmI+j*ld,Troot+j*ld,PetscMin(j+2,n)));
    TrootmI[j+j*ld] -= 1.0;
  }
  PetscCall(SlepcNormAm(n,TrootmI,2,work,rand,&d2));
  d2 = PetscPowReal(d2,1.0/2.0);
  PetscCall(SlepcNormAm(n,TrootmI,3,work,rand,&d3));
  d3 = PetscPowReal(d3,1.0/3.0);
  a2 = PetscMax(d2,d3);
  if (a2 <= xvals[1]) {
    *m = (a2 <= xvals[0])? 1: 2;
    foundm = PETSC_TRUE;
  }
  p = 0;
  while (!foundm) {
    more = PETSC_FALSE;  /* More norm checks needed */
    if (*s > s0) {
      PetscCall(SlepcNormAm(n,TrootmI,3,work,rand,&d3));
      d3 = PetscPowReal(d3,1.0/3.0);
    }
    PetscCall(SlepcNormAm(n,TrootmI,4,work,rand,&d4));
    d4 = PetscPowReal(d4,1.0/4.0);
    a3 = PetscMax(d3,d4);
    if (a3 <= xvals[mmax-1]) {
      for (j=2;j<mmax;j++) if (a3 <= xvals[j]) break;
      if (j <= 5) {
        *m = j+1;
        break;
      } else if (a3/2.0 <= xvals[4] && p < 2) {
        more = PETSC_TRUE;
        p = p + 1;
      }
    }
    if (!more) {
      PetscCall(SlepcNormAm(n,TrootmI,5,work,rand,&d5));
      d5 = PetscPowReal(d5,1.0/5.0);
      a4 = PetscMax(d4,d5);
      eta = PetscMin(a3,a4);
      if (eta <= xvals[mmax-1]) {
        for (j=5;j<mmax;j++) if (eta <= xvals[j]) break;
        *m = j + 1;
        break;
      }
    }
    if (*s == maxroots) {
      PetscCall(PetscInfo(fn,"Too many matrix square roots\n"));
      *m = mmax;  /* No good value found so take largest */
      break;
    }
    PetscCall(FNSqrtmSchur(fn,n,Troot,ld,PETSC_FALSE));
    /* TrootmI = Troot - I */
    for (j=0;j<n;j++) {
      PetscCall(PetscArraycpy(TrootmI+j*ld,Troot+j*ld,PetscMin(j+2,n)));
      TrootmI[j+j*ld] -= 1.0;
    }
    *s = *s + 1;
  }
  PetscCall(PetscRandomDestroy(&rand));
  PetscFunctionReturn(PETSC_SUCCESS);
}

#if !defined(PETSC_USE_COMPLEX)
/*
   Computes a^(1/2^s) - 1 accurately, avoiding subtractive cancellation
*/
static PetscScalar sqrt_obo(PetscScalar a,PetscInt s)
{
  PetscScalar val,z0,r;
  PetscReal   angle;
  PetscInt    i,n0;

  PetscFunctionBegin;
  if (s == 0) val = a-1.0;
  else {
    n0 = s;
    angle = PetscAtan2Real(PetscImaginaryPart(a),PetscRealPart(a));
    if (angle >= PETSC_PI/2.0) {
      a = PetscSqrtScalar(a);
      n0 = s - 1;
    }
    z0 = a - 1.0;
    a = PetscSqrtScalar(a);
    r = 1.0 + a;
    for (i=0;i<n0-1;i++) {
      a = PetscSqrtScalar(a);
      r = r*(1.0+a);
    }
    val = z0/r;
  }
  PetscFunctionReturn(val);
}
#endif

/*
   Square root of 2x2 matrix T from block diagonal of Schur form. Overwrites T.
*/
static PetscErrorCode sqrtm_tbt(PetscScalar *T)
{
  PetscScalar t11,t12,t21,t22,r11,r22;
#if !defined(PETSC_USE_COMPLEX)
  PetscScalar mu;
#endif

  PetscFunctionBegin;
  t11 = T[0]; t21 = T[1]; t12 = T[2]; t22 = T[3];
  if (t21 != 0.0) {
    /* Compute square root of 2x2 quasitriangular block */
    /* The algorithm assumes the special structure of real Schur form */
#if defined(PETSC_USE_COMPLEX)
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_PLIB,"Should not reach this line in complex scalars");
#else
    mu = PetscSqrtReal(-t21*t12);
    if (t11 > 0.0) r11 = PetscSqrtReal((t11+SlepcAbsEigenvalue(t11,mu))/2.0);
    else r11 = mu / PetscSqrtReal(2.0*(-t11+SlepcAbsEigenvalue(t11,mu)));
    T[0] = r11;
    T[1] = t21/(2.0*r11);
    T[2] = t12/(2.0*r11);
    T[3] = r11;
#endif
  } else {
    /* Compute square root of 2x2 upper triangular block */
    r11 = PetscSqrtScalar(t11);
    r22 = PetscSqrtScalar(t22);
    T[0] = r11;
    T[2] = t12/(r11+r22);
    T[3] = r22;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

#if defined(PETSC_USE_COMPLEX)
/*
   Unwinding number of z
*/
static inline PetscReal unwinding(PetscScalar z)
{
  PetscReal u;

  PetscFunctionBegin;
  u = PetscCeilReal((PetscImaginaryPart(z)-PETSC_PI)/(2.0*PETSC_PI));
  PetscFunctionReturn(u);
}
#endif

/*
   Power of 2-by-2 upper triangular matrix A.
   Returns the (1,2) element of the pth power of A, where p is an arbitrary real number
*/
static PetscScalar powerm2by2(PetscScalar A11,PetscScalar A22,PetscScalar A12,PetscReal p)
{
  PetscScalar a1,a2,a1p,a2p,loga1,loga2,w,dd,x12;

  PetscFunctionBegin;
  a1 = A11;
  a2 = A22;
  if (a1 == a2) {
    x12 = p*A12*PetscPowScalarReal(a1,p-1);
  } else if (PetscAbsScalar(a1) < 0.5*PetscAbsScalar(a2) || PetscAbsScalar(a2) < 0.5*PetscAbsScalar(a1)) {
    a1p = PetscPowScalarReal(a1,p);
    a2p = PetscPowScalarReal(a2,p);
    x12 = A12*(a2p-a1p)/(a2-a1);
  } else {  /* Close eigenvalues */
    loga1 = PetscLogScalar(a1);
    loga2 = PetscLogScalar(a2);
    w = PetscAtanhScalar((a2-a1)/(a2+a1));
#if defined(PETSC_USE_COMPLEX)
    w += PETSC_i*PETSC_PI*unwinding(loga2-loga1);
#endif
    dd = 2.0*PetscExpScalar((loga1+loga2)*p/2.0)*PetscSinhScalar(p*w)/(a2-a1);
    x12 = A12*dd;
  }
  PetscFunctionReturn(x12);
}

/*
   Recomputes diagonal blocks of T = X^(1/2^s) - 1 more accurately
*/
static PetscErrorCode recompute_diag_blocks_sqrt(PetscBLASInt n,PetscScalar *Troot,PetscScalar *T,PetscBLASInt ld,PetscInt *blockStruct,PetscInt s)
{
  PetscScalar    A[4],P[4],M[4],Z0[4],det;
  PetscInt       i,j;
#if !defined(PETSC_USE_COMPLEX)
  PetscInt       last_block=0;
  PetscScalar    a;
#endif

  PetscFunctionBegin;
  for (j=0;j<n-1;j++) {
#if !defined(PETSC_USE_COMPLEX)
    switch (blockStruct[j]) {
      case 0: /* Not start of a block */
        if (last_block != 0) {
          last_block = 0;
        } else { /* In a 1x1 block */
          a = T[j+j*ld];
          Troot[j+j*ld] = sqrt_obo(a,s);
        }
        break;
      default: /* In a 2x2 block */
        last_block = blockStruct[j];
#endif
        if (s == 0) {
          Troot[j+j*ld]       = T[j+j*ld]-1.0;
          Troot[j+1+j*ld]     = T[j+1+j*ld];
          Troot[j+(j+1)*ld]   = T[j+(j+1)*ld];
          Troot[j+1+(j+1)*ld] = T[j+1+(j+1)*ld]-1.0;
          continue;
        }
        A[0] = T[j+j*ld]; A[1] = T[j+1+j*ld]; A[2] = T[j+(j+1)*ld]; A[3] = T[j+1+(j+1)*ld];
        PetscCall(sqrtm_tbt(A));
        /* Z0 = A - I */
        Z0[0] = A[0]-1.0; Z0[1] = A[1]; Z0[2] = A[2]; Z0[3] = A[3]-1.0;
        if (s == 1) {
          Troot[j+j*ld]       = Z0[0];
          Troot[j+1+j*ld]     = Z0[1];
          Troot[j+(j+1)*ld]   = Z0[2];
          Troot[j+1+(j+1)*ld] = Z0[3];
          continue;
        }
        PetscCall(sqrtm_tbt(A));
        /* P = A + I */
        P[0] = A[0]+1.0; P[1] = A[1]; P[2] = A[2]; P[3] = A[3]+1.0;
        for (i=0;i<s-2;i++) {
          PetscCall(sqrtm_tbt(A));
          /* P = P*(I + A) */
          M[0] = P[0]*(A[0]+1.0)+P[2]*A[1];
          M[1] = P[1]*(A[0]+1.0)+P[3]*A[1];
          M[2] = P[0]*A[2]+P[2]*(A[3]+1.0);
          M[3] = P[1]*A[2]+P[3]*(A[3]+1.0);
          P[0] = M[0]; P[1] = M[1]; P[2] = M[2]; P[3] = M[3];
        }
        /* Troot(j:j+1,j:j+1) = Z0 / P (via Cramer) */
        det = P[0]*P[3]-P[1]*P[2];
        Troot[j+j*ld]       = (Z0[0]*P[3]-P[1]*Z0[2])/det;
        Troot[j+(j+1)*ld]   = (P[0]*Z0[2]-Z0[0]*P[2])/det;
        Troot[j+1+j*ld]     = (Z0[1]*P[3]-P[1]*Z0[3])/det;
        Troot[j+1+(j+1)*ld] = (P[0]*Z0[3]-Z0[1]*P[2])/det;
        /* If block is upper triangular recompute the (1,2) element.
           Skip when T(j,j) or T(j+1,j+1) < 0 since the implementation of atanh is nonstandard */
        if (T[j+1+j*ld]==0.0 && PetscRealPart(T[j+j*ld])>=0.0 && PetscRealPart(T[j+1+(j+1)*ld])>=0.0) {
          Troot[j+(j+1)*ld] = powerm2by2(T[j+j*ld],T[j+1+(j+1)*ld],T[j+(j+1)*ld],1.0/PetscPowInt(2,s));
        }
#if !defined(PETSC_USE_COMPLEX)
    }
#endif
  }
#if !defined(PETSC_USE_COMPLEX)
  /* If last diagonal entry is not in a block it will have been missed */
  if (blockStruct[n-2] == 0) {
    a = T[n-1+(n-1)*ld];
    Troot[n-1+(n-1)*ld] = sqrt_obo(a,s);
  }
#endif
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   Nodes x and weights w for n-point Gauss-Legendre quadrature (Q is n*n workspace)

   G. H. Golub and J. H. Welsch, Calculation of Gauss quadrature
      rules, Math. Comp., 23(106):221-230, 1969.
*/
static PetscErrorCode gauss_legendre(PetscBLASInt n,PetscScalar *x,PetscScalar *w,PetscScalar *Q)
{
  PetscScalar    v,a,*work;
  PetscReal      *eig,dummy;
  PetscBLASInt   k,ld=n,lwork,info;
#if defined(PETSC_USE_COMPLEX)
  PetscReal      *rwork,rdummy;
#endif

  PetscFunctionBegin;
  PetscCall(PetscArrayzero(Q,n*n));
  for (k=1;k<n;k++) {
    v = k/PetscSqrtReal(4.0*k*k-1.0);
    Q[k+(k-1)*n] = v;
    Q[(k-1)+k*n] = v;
  }

  /* workspace query and memory allocation */
  lwork = -1;
#if defined(PETSC_USE_COMPLEX)
  PetscCallBLAS("LAPACKsyev",LAPACKsyev_("V","L",&n,Q,&ld,&dummy,&a,&lwork,&rdummy,&info));
  PetscCall(PetscBLASIntCast((PetscInt)PetscRealPart(a),&lwork));
  PetscCall(PetscMalloc3(n,&eig,lwork,&work,PetscMax(1,3*n-2),&rwork));
#else
  PetscCallBLAS("LAPACKsyev",LAPACKsyev_("V","L",&n,Q,&ld,&dummy,&a,&lwork,&info));
  PetscCall(PetscBLASIntCast((PetscInt)a,&lwork));
  PetscCall(PetscMalloc2(n,&eig,lwork,&work));
#endif

  /* compute eigendecomposition */
#if defined(PETSC_USE_COMPLEX)
  PetscCallBLAS("LAPACKsyev",LAPACKsyev_("V","L",&n,Q,&ld,eig,work,&lwork,rwork,&info));
#else
  PetscCallBLAS("LAPACKsyev",LAPACKsyev_("V","L",&n,Q,&ld,eig,work,&lwork,&info));
#endif
  SlepcCheckLapackInfo("syev",info);

  for (k=0;k<n;k++) {
    x[k] = eig[k];
    w[k] = 2.0*Q[k*n]*Q[k*n];
  }
#if defined(PETSC_USE_COMPLEX)
  PetscCall(PetscFree3(eig,work,rwork));
#else
  PetscCall(PetscFree2(eig,work));
#endif
  PetscCall(PetscLogFlops(9.0*n*n*n+2.0*n*n*n));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   Pade approximation to log(1 + T) via partial fractions
*/
static PetscErrorCode pade_approx(PetscBLASInt n,PetscScalar *T,PetscScalar *L,PetscBLASInt ld,PetscInt m,PetscScalar *work)
{
  PetscScalar    *K,*W,*nodes,*wts;
  PetscBLASInt   *ipiv,info,mm;
  PetscInt       i,j,k;

  PetscFunctionBegin;
  K     = work;
  W     = work+n*n;
  nodes = work+2*n*n;
  wts   = work+2*n*n+m;
  ipiv  = (PetscBLASInt*)(work+2*n*n+2*m);
  PetscCall(PetscBLASIntCast(m,&mm));
  PetscCall(gauss_legendre(mm,nodes,wts,L));
  /* Convert from [-1,1] to [0,1] */
  for (i=0;i<m;i++) {
    nodes[i] = (nodes[i]+1.0)/2.0;
    wts[i] = wts[i]/2.0;
  }
  PetscCall(PetscArrayzero(L,n*n));
  for (k=0;k<m;k++) {
    for (i=0;i<n;i++) for (j=0;j<n;j++) K[i+j*ld] = nodes[k]*T[i+j*ld];
    for (i=0;i<n;i++) K[i+i*ld] += 1.0;
    for (i=0;i<n;i++) for (j=0;j<n;j++) W[i+j*ld] = T[i+j*ld];
    PetscCallBLAS("LAPACKgesv",LAPACKgesv_(&n,&n,K,&n,ipiv,W,&n,&info));
    for (i=0;i<n;i++) for (j=0;j<n;j++) L[i+j*ld] += wts[k]*W[i+j*ld];
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   Recomputes diagonal blocks of T = X^(1/2^s) - 1 more accurately
*/
static PetscErrorCode recompute_diag_blocks_log(PetscBLASInt n,PetscScalar *L,PetscScalar *T,PetscBLASInt ld,PetscInt *blockStruct)
{
  PetscScalar a1,a2,a12,loga1,loga2,z,dd;
  PetscInt    j;
#if !defined(PETSC_USE_COMPLEX)
  PetscInt    last_block=0;
  PetscScalar f,t;
#endif

  PetscFunctionBegin;
  for (j=0;j<n-1;j++) {
#if !defined(PETSC_USE_COMPLEX)
    switch (blockStruct[j]) {
      case 0: /* Not start of a block */
        if (last_block != 0) {
          last_block = 0;
        } else { /* In a 1x1 block */
          L[j+j*ld] = PetscLogScalar(T[j+j*ld]);
        }
        break;
      case 1: /* Start of upper-tri block */
        last_block = 1;
#endif
        a1 = T[j+j*ld];
        a2 = T[j+1+(j+1)*ld];
        loga1 = PetscLogScalar(a1);
        loga2 = PetscLogScalar(a2);
        L[j+j*ld] = loga1;
        L[j+1+(j+1)*ld] = loga2;
        if ((PetscRealPart(a1)<0.0 && PetscImaginaryPart(a1)==0.0) || (PetscRealPart(a2)<0.0 && PetscImaginaryPart(a1)==0.0)) {
         /* Problems with 2 x 2 formula for (1,2) block
            since atanh is nonstandard, just redo diagonal part */
          continue;
        }
        if (a1 == a2) {
          a12 = T[j+(j+1)*ld]/a1;
        } else if (PetscAbsScalar(a1)<0.5*PetscAbsScalar(a2) || PetscAbsScalar(a2)<0.5*PetscAbsScalar(a1)) {
          a12 = T[j+(j+1)*ld]*(loga2-loga1)/(a2-a1);
        } else {  /* Close eigenvalues */
          z = (a2-a1)/(a2+a1);
          dd = 2.0*PetscAtanhScalar(z);
#if defined(PETSC_USE_COMPLEX)
          dd += 2.0*PETSC_i*PETSC_PI*unwinding(loga2-loga1);
#endif
          dd /= (a2-a1);
          a12 = T[j+(j+1)*ld]*dd;
        }
        L[j+(j+1)*ld] = a12;
#if !defined(PETSC_USE_COMPLEX)
        break;
      case 2: /* Start of quasi-tri block */
        last_block = 2;
        f = 0.5*PetscLogScalar(T[j+j*ld]*T[j+j*ld]-T[j+(j+1)*ld]*T[j+1+j*ld]);
        t = PetscAtan2Real(PetscSqrtScalar(-T[j+(j+1)*ld]*T[j+1+j*ld]),T[j+j*ld])/PetscSqrtScalar(-T[j+(j+1)*ld]*T[j+1+j*ld]);
        L[j+j*ld]       = f;
        L[j+1+j*ld]     = t*T[j+1+j*ld];
        L[j+(j+1)*ld]   = t*T[j+(j+1)*ld];
        L[j+1+(j+1)*ld] = f;
    }
#endif
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}
/*
 * Matrix logarithm implementation based on algorithm and matlab code by N. Higham and co-authors
 *
 *     H. Al-Mohy and N. J. Higham, "Improved inverse scaling and squaring
 *     algorithms for the matrix logarithm", SIAM J. Sci. Comput. 34(4):C153-C169, 2012.
 */
static PetscErrorCode FNLogmPade(FN fn,PetscBLASInt n,PetscScalar *T,PetscBLASInt ld,PetscBool firstonly)
{
  PetscBLASInt   k,sdim,lwork,info;
  PetscScalar    *wr,*wi=NULL,*W,*Q,*Troot,*L,*work,one=1.0,zero=0.0,alpha;
  PetscInt       i,j,s=0,m=0,*blockformat;
#if defined(PETSC_USE_COMPLEX)
  PetscReal      *rwork;
#endif

  PetscFunctionBegin;
  lwork = 3*n*n; /* gees needs only 5*n, but work is also passed to FNlogm_params */
  k     = firstonly? 1: n;

  /* compute Schur decomposition A*Q = Q*T */
  PetscCall(PetscCalloc7(n,&wr,n*k,&W,n*n,&Q,n*n,&Troot,n*n,&L,lwork,&work,n-1,&blockformat));
#if !defined(PETSC_USE_COMPLEX)
  PetscCall(PetscMalloc1(n,&wi));
  PetscCallBLAS("LAPACKgees",LAPACKgees_("V","N",NULL,&n,T,&ld,&sdim,wr,wi,Q,&ld,work,&lwork,NULL,&info));
#else
  PetscCall(PetscMalloc1(n,&rwork));
  PetscCallBLAS("LAPACKgees",LAPACKgees_("V","N",NULL,&n,T,&ld,&sdim,wr,Q,&ld,work,&lwork,rwork,NULL,&info));
#endif
  SlepcCheckLapackInfo("gees",info);

#if !defined(PETSC_USE_COMPLEX)
  /* check for negative real eigenvalues */
  for (i=0;i<n;i++) {
    PetscCheck(wr[i]>=0.0 || wi[i]!=0.0,PETSC_COMM_SELF,PETSC_ERR_SUP,"Matrix has negative real eigenvalue; rerun with complex scalars");
  }
#endif

  /* get block structure of Schur factor */
  PetscCall(qtri_struct(n,T,ld,blockformat));

  /* get parameters */
  PetscCall(FNlogm_params(fn,n,T,ld,wr,wi,100,&s,&m,Troot,work));

  /* compute Troot - I = T(1/2^s) - I more accurately */
  PetscCall(recompute_diag_blocks_sqrt(n,Troot,T,ld,blockformat,s));

  /* compute Pade approximant */
  PetscCall(pade_approx(n,Troot,L,ld,m,work));

  /* scale back up, L = 2^s * L */
  alpha = PetscPowInt(2,s);
  for (i=0;i<n;i++) for (j=0;j<n;j++) L[i+j*ld] *= alpha;

  /* recompute diagonal blocks */
  PetscCall(recompute_diag_blocks_log(n,L,T,ld,blockformat));

  /* backtransform B = Q*L*Q' */
  PetscCallBLAS("BLASgemm",BLASgemm_("N","C",&n,&k,&n,&one,L,&ld,Q,&ld,&zero,W,&ld));
  PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&n,&k,&n,&one,Q,&ld,W,&ld,&zero,T,&ld));

  /* flop count: Schur decomposition, and backtransform */
  PetscCall(PetscLogFlops(25.0*n*n*n+4.0*n*n*k));

  PetscCall(PetscFree7(wr,W,Q,Troot,L,work,blockformat));
#if !defined(PETSC_USE_COMPLEX)
  PetscCall(PetscFree(wi));
#else
  PetscCall(PetscFree(rwork));
#endif
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode FNEvaluateFunctionMat_Log_Higham(FN fn,Mat A,Mat B)
{
  PetscBLASInt   n = 0;
  PetscScalar    *T;
  PetscInt       m;

  PetscFunctionBegin;
  if (A!=B) PetscCall(MatCopy(A,B,SAME_NONZERO_PATTERN));
  PetscCall(MatDenseGetArray(B,&T));
  PetscCall(MatGetSize(A,&m,NULL));
  PetscCall(PetscBLASIntCast(m,&n));
  PetscCall(FNLogmPade(fn,n,T,n,PETSC_FALSE));
  PetscCall(MatDenseRestoreArray(B,&T));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode FNEvaluateFunctionMatVec_Log_Higham(FN fn,Mat A,Vec v)
{
  PetscBLASInt   n = 0;
  PetscScalar    *T;
  PetscInt       m;
  Mat            B;

  PetscFunctionBegin;
  PetscCall(FN_AllocateWorkMat(fn,A,&B));
  PetscCall(MatDenseGetArray(B,&T));
  PetscCall(MatGetSize(A,&m,NULL));
  PetscCall(PetscBLASIntCast(m,&n));
  PetscCall(FNLogmPade(fn,n,T,n,PETSC_TRUE));
  PetscCall(MatDenseRestoreArray(B,&T));
  PetscCall(MatGetColumnVector(B,v,0));
  PetscCall(FN_FreeWorkMat(fn,&B));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode FNView_Log(FN fn,PetscViewer viewer)
{
  PetscBool      isascii;
  char           str[50];
  const char     *methodname[] = {
                  "scaling & squaring, [m/m] Pade approximant (Higham)"
  };
  const int      nmeth=PETSC_STATIC_ARRAY_LENGTH(methodname);

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) {
    if (fn->beta==(PetscScalar)1.0) {
      if (fn->alpha==(PetscScalar)1.0) PetscCall(PetscViewerASCIIPrintf(viewer,"  logarithm: log(x)\n"));
      else {
        PetscCall(SlepcSNPrintfScalar(str,sizeof(str),fn->alpha,PETSC_TRUE));
        PetscCall(PetscViewerASCIIPrintf(viewer,"  logarithm: log(%s*x)\n",str));
      }
    } else {
      PetscCall(SlepcSNPrintfScalar(str,sizeof(str),fn->beta,PETSC_TRUE));
      if (fn->alpha==(PetscScalar)1.0) PetscCall(PetscViewerASCIIPrintf(viewer,"  logarithm: %s*log(x)\n",str));
      else {
        PetscCall(PetscViewerASCIIPrintf(viewer,"  logarithm: %s",str));
        PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_FALSE));
        PetscCall(SlepcSNPrintfScalar(str,sizeof(str),fn->alpha,PETSC_TRUE));
        PetscCall(PetscViewerASCIIPrintf(viewer,"*log(%s*x)\n",str));
        PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_TRUE));
      }
    }
    if (fn->method<nmeth) PetscCall(PetscViewerASCIIPrintf(viewer,"  computing matrix functions with: %s\n",methodname[fn->method]));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode FNCreate_Log(FN fn)
{
  PetscFunctionBegin;
  fn->ops->evaluatefunction          = FNEvaluateFunction_Log;
  fn->ops->evaluatederivative        = FNEvaluateDerivative_Log;
  fn->ops->evaluatefunctionmat[0]    = FNEvaluateFunctionMat_Log_Higham;
  fn->ops->evaluatefunctionmatvec[0] = FNEvaluateFunctionMatVec_Log_Higham;
  fn->ops->view                      = FNView_Log;
  PetscFunctionReturn(PETSC_SUCCESS);
}
