/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   HZ iteration for generalized symmetric-indefinite eigenproblem.
   Based on Matlab code from David Watkins.

   References:

       [1] D.S. Watkins, The Matrix Eigenvalue Problem: GR and Krylov Subspace
           Methods, SIAM, 2007.

       [2] M.A. Brebner, J. Grad, "Eigenvalues of Ax = lambda Bx for real
           symmetric matrices A and B computed by reduction to pseudosymmetric
           form and the HR process", Linear Alg. Appl. 43:99-118, 1982.
*/

#include <slepc/private/dsimpl.h>
#include <slepcblaslapack.h>

/*
   Sets up a 2-by-2 matrix to eliminate y in the vector [x y]'.
   Transformation is rotator if sygn = 1 and hyperbolic if sygn = -1.
*/
static PetscErrorCode UnifiedRotation(PetscReal x,PetscReal y,PetscReal sygn,PetscReal *rot,PetscReal *rcond,PetscBool *swap)
{
  PetscReal nrm,c,s;

  PetscFunctionBegin;
  *swap = PETSC_FALSE;
  if (y == 0) {
    rot[0] = 1.0; rot[1] = 0.0; rot[2] = 0.0; rot[3] = 1.0;
    *rcond = 1.0;
  } else {
    nrm = PetscMax(PetscAbs(x),PetscAbs(y));
    c = x/nrm; s = y/nrm;
    PetscCheck(sygn==1.0 || sygn==-1,PETSC_COMM_SELF,PETSC_ERR_PLIB,"Value of sygn sent to transetup must be 1 or -1");
    if (sygn == 1.0) {  /* set up a rotator */
      nrm = SlepcAbs(c,s);
      c = c/nrm; s = s/nrm;
      /* rot = [c s; -s c]; */
      rot[0] = c; rot[1] = -s; rot[2] = s; rot[3] = c;
      *rcond = 1.0;
    } else {  /* sygn == -1, set up a hyperbolic transformation */
      nrm = c*c-s*s;
      if (nrm > 0) nrm = PetscSqrtReal(nrm);
      else {
        PetscCheck(nrm<0,PETSC_COMM_SELF,PETSC_ERR_PLIB,"Breakdown in construction of hyperbolic transformation");
        nrm = PetscSqrtReal(-nrm);
        *swap = PETSC_TRUE;
      }
      c = c/nrm; s = s/nrm;
      /* rot = [c -s; -s c]; */
      rot[0] = c; rot[1] = -s; rot[2] = -s; rot[3] = c;
      *rcond = PetscAbs(PetscAbs(s)-PetscAbs(c))/(PetscAbs(s)+PetscAbs(c));
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode HZStep(PetscBLASInt ntop,PetscBLASInt nn,PetscReal tr,PetscReal dt,PetscReal *aa,PetscReal *bb,PetscReal *dd,PetscScalar *uu,PetscInt n,PetscInt ld,PetscBool *flag)
{
  PetscBLASInt   one=1;
  PetscInt       k,jj,ii;
  PetscBLASInt   n_;
  PetscReal      bulge10,bulge20,bulge30,bulge31,bulge41,bulge42;
  PetscReal      sygn,rcond=1.0,worstcond,rot[4],buf[2],t;
  PetscScalar    rtmp;
  PetscBool      swap;

  PetscFunctionBegin;
  *flag = PETSC_FALSE;
  worstcond = 1.0;
  PetscCall(PetscBLASIntCast(n,&n_));

  /* Build initial bulge that sets step in motion */
  bulge10 = dd[ntop+1]*(aa[ntop]*(aa[ntop] - dd[ntop]*tr) + dt*dd[ntop]*dd[ntop]) + dd[ntop]*bb[ntop]*bb[ntop];
  bulge20 = bb[ntop]*(dd[ntop+1]*aa[ntop] + dd[ntop]*aa[ntop+1] - dd[ntop]*dd[ntop+1]*tr);
  bulge30 = bb[ntop]*bb[ntop+1]*dd[ntop];
  bulge31 = 0.0;
  bulge41 = 0.0;
  bulge42 = 0.0;

  /* Chase the bulge */
  for (jj=ntop;jj<nn-1;jj++) {

    /* Check for trivial bulge */
    if (jj>ntop && PetscMax(PetscMax(PetscAbs(bulge10),PetscAbs(bulge20)),PetscAbs(bulge30))<PETSC_MACHINE_EPSILON*(PetscAbs(aa[jj]) + PetscAbs(aa[jj+1]))) {
      bb[jj-1] = 0.0;  /* deflate and move on */

    } else { /* carry out the step */

      /* Annihilate tip entry bulge30 */
      if (bulge30 != 0.0) {

        /* Make an interchange if necessary to ensure that the
           first transformation is othogonal, not hyperbolic.  */
        if (dd[jj+1] != dd[jj+2]) { /* make an interchange */
          if (dd[jj] != dd[jj+1]) {  /* interchange 1st and 2nd */
            buf[0] = bulge20; bulge20 = bulge10; bulge10 = buf[0];
            buf[0] = aa[jj]; aa[jj] = aa[jj+1]; aa[jj+1] = buf[0];
            buf[0] = bb[jj+1]; bb[jj+1] = bulge31; bulge31 = buf[0];
            buf[0] = dd[jj]; dd[jj] = dd[jj+1]; dd[jj+1] = buf[0];
            for (k=0;k<n;k++) {
              rtmp = uu[k+jj*ld]; uu[k+jj*ld] = uu[k+(jj+1)*ld]; uu[k+(jj+1)*ld] = rtmp;
            }
          } else {  /* interchange 1st and 3rd */
            buf[0] = bulge30; bulge30 = bulge10; bulge10 = buf[0];
            buf[0] = aa[jj]; aa[jj] = aa[jj+2]; aa[jj+2] = buf[0];
            buf[0] = bb[jj]; bb[jj] = bb[jj+1]; bb[jj+1] = buf[0];
            buf[0] = dd[jj]; dd[jj] = dd[jj+2]; dd[jj+2] = buf[0];
            if (jj + 2 < nn-1) {
              bulge41 = bb[jj+2];
              bb[jj+2] = 0;
            }
            for (k=0;k<n;k++) {
              rtmp = uu[k+jj*ld]; uu[k+jj*ld] = uu[k+(jj+2)*ld]; uu[k+(jj+2)*ld] = rtmp;
            }
          }
        }

        /* Set up transforming matrix rot. */
        PetscCall(UnifiedRotation(bulge20,bulge30,1,rot,&rcond,&swap));

        /* Apply transforming matrix rot to T. */
        bulge20 = rot[0]*bulge20 + rot[2]*bulge30;
        buf[0] = rot[0]*bb[jj] + rot[2]*bulge31;
        buf[1] = rot[1]*bb[jj] + rot[3]*bulge31;
        bb[jj] = buf[0];
        bulge31 = buf[1];
        buf[0] = rot[0]*rot[0]*aa[jj+1] + 2.0*rot[0]*rot[2]*bb[jj+1] + rot[2]*rot[2]*aa[jj+2];
        buf[1] = rot[1]*rot[1]*aa[jj+1] + 2.0*rot[3]*rot[1]*bb[jj+1] + rot[3]*rot[3]*aa[jj+2];
        bb[jj+1] = rot[1]*rot[0]*aa[jj+1] + rot[3]*rot[2]*aa[jj+2] + (rot[3]*rot[0] + rot[1]*rot[2])*bb[jj+1];
        aa[jj+1] = buf[0];
        aa[jj+2] = buf[1];
        if (jj + 2 < nn-1) {
          bulge42 = bb[jj+2]*rot[2];
          bb[jj+2] = bb[jj+2]*rot[3];
        }

        /* Accumulate transforming matrix */
        PetscCallBLAS("BLASrot",BLASMIXEDrot_(&n_,uu+(jj+1)*ld,&one,uu+(jj+2)*ld,&one,&rot[0],&rot[2]));
      }

      /* Annihilate inner entry bulge20 */
      if (bulge20 != 0.0) {

        /* Begin by setting up transforming matrix rot */
        sygn = dd[jj]*dd[jj+1];
        PetscCall(UnifiedRotation(bulge10,bulge20,sygn,rot,&rcond,&swap));
        if (rcond<PETSC_MACHINE_EPSILON) {
          *flag = PETSC_TRUE;
          PetscFunctionReturn(PETSC_SUCCESS);
        }
        if (rcond < worstcond) worstcond = rcond;

        /* Apply transforming matrix rot to T */
        if (jj > ntop) bb[jj-1] = rot[0]*bulge10 + rot[2]*bulge20;
        buf[0] = rot[0]*rot[0]*aa[jj] + 2*rot[0]*rot[2]*bb[jj] + rot[2]*rot[2]*aa[jj+1];
        buf[1] = rot[1]*rot[1]*aa[jj] + 2*rot[3]*rot[1]*bb[jj] + rot[3]*rot[3]*aa[jj+1];
        bb[jj] = rot[1]*rot[0]*aa[jj] + rot[3]*rot[2]*aa[jj+1] + (rot[3]*rot[0] + rot[1]*rot[2])*bb[jj];
        aa[jj] = buf[0];
        aa[jj+1] = buf[1];
        if (jj + 1 < nn-1) {
          /* buf = [ bulge31 bb(jj+1) ] * rot' */
          buf[0] = rot[0]*bulge31 + rot[2]*bb[jj+1];
          buf[1] = rot[1]*bulge31 + rot[3]*bb[jj+1];
          bulge31 = buf[0];
          bb[jj+1] = buf[1];
        }
        if (jj + 2 < nn-1) {
          /* buf = [bulge41 bulge42] * rot' */
          buf[0] = rot[0]*bulge41 + rot[2]*bulge42;
          buf[1] = rot[1]*bulge41 + rot[3]*bulge42;
          bulge41 = buf[0];
          bulge42 = buf[1];
        }

        /* Apply transforming matrix rot to D */
        if (swap == 1) {
          buf[0] = dd[jj]; dd[jj] = dd[jj+1]; dd[jj+1] = buf[0];
        }

        /* Accumulate transforming matrix, uu(jj:jj+1,:) = rot*uu(jj:jj+1,:) */
        if (sygn==1) {
          PetscCallBLAS("BLASrot",BLASMIXEDrot_(&n_,uu+jj*ld,&one,uu+(jj+1)*ld,&one,&rot[0],&rot[2]));
        } else {
          if (PetscAbsReal(rot[0])>PetscAbsReal(rot[1])) { /* Type I */
            t = rot[1]/rot[0];
            for (ii=0;ii<n;ii++) {
              uu[jj*ld+ii] = rot[0]*uu[jj*ld+ii] + rot[1]*uu[(jj+1)*ld+ii];
              uu[(jj+1)*ld+ii] = t*uu[jj*ld+ii] + uu[(jj+1)*ld+ii]/rot[0];
            }
          } else { /* Type II */
            t = rot[0]/rot[1];
            for (ii=0;ii<n;ii++) {
              rtmp = uu[jj*ld+ii];
              uu[jj*ld+ii] = rot[0]*uu[jj*ld+ii] + rot[1]*uu[(jj+1)*ld+ii];
              uu[(jj+1)*ld+ii] = t*uu[jj*ld+ii] + rtmp/rot[1];
            }
          }
        }
      }
    }

    /* Adjust bulge for next step */
    bulge10 = bb[jj];
    bulge20 = bulge31;
    bulge30 = bulge41;
    bulge31 = bulge42;
    bulge41 = 0.0;
    bulge42 = 0.0;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode HZIteration(PetscBLASInt nn,PetscBLASInt cgd,PetscReal *aa,PetscReal *bb,PetscReal *dd,PetscScalar *uu,PetscBLASInt ld)
{
  PetscBLASInt   j2,one=1,its,nits,nstop,jj,ntop,nbot,ntry;
  PetscReal      htr,det,dis,dif,tn,kt,c,s,tr,dt;
  PetscBool      flag=PETSC_FALSE;

  PetscFunctionBegin;
  its = 0;
  nbot = nn-1;
  nits = 0;
  nstop = 40*(nn - cgd);

  while (nbot >= cgd && nits < nstop) {

    /* Check for zeros on the subdiagonal */
    jj = nbot - 1;
    while (jj>=cgd && PetscAbs(bb[jj])>PETSC_MACHINE_EPSILON*(PetscAbs(aa[jj])+PetscAbs(aa[jj+1]))) jj = jj-1;
    if (jj>=cgd) bb[jj]=0;
    ntop = jj + 1;  /* starting point for step */
    if (ntop == nbot) {  /* isolate single eigenvalue */
      nbot = ntop - 1;
      its = 0;
    } else if (ntop+1 == nbot) {  /* isolate pair of eigenvalues */
      htr = 0.5*(aa[ntop]*dd[ntop] + aa[nbot]*dd[nbot]);
      det = dd[ntop]*dd[nbot]*(aa[ntop]*aa[nbot]-bb[ntop]*bb[ntop]);
      dis = htr*htr - det;
      if (dis > 0) {  /* distinct real eigenvalues */
        if (dd[ntop] == dd[nbot]) {  /* separate the eigenvalues by a Jacobi rotator */
          dif = aa[ntop]-aa[nbot];
          if (2.0*PetscAbs(bb[ntop])<=dif) {
            tn = 2*bb[ntop]/dif;
            tn = tn/(1.0 + PetscSqrtReal(1.0+tn*tn));
          } else {
            kt = dif/(2.0*bb[ntop]);
            tn = PetscSign(kt)/(PetscAbsReal(kt)+PetscSqrtReal(1.0+kt*kt));
          }
          c = 1.0/PetscSqrtReal(1.0 + tn*tn);
          s = c*tn;
          aa[ntop] = aa[ntop] + tn*bb[ntop];
          aa[nbot] = aa[nbot] - tn*bb[ntop];
          bb[ntop] = 0;
          j2 = nn-cgd;
          PetscCallBLAS("BLASrot",BLASMIXEDrot_(&j2,uu+ntop*ld+cgd,&one,uu+nbot*ld+cgd,&one,&c,&s));
        }
      }
      nbot = ntop - 1;
    } else {  /* Do an HZ iteration */
      its = its + 1;
      nits = nits + 1;
      tr = aa[nbot-1]*dd[nbot-1] + aa[nbot]*dd[nbot];
      dt = dd[nbot-1]*dd[nbot]*(aa[nbot-1]*aa[nbot]-bb[nbot-1]*bb[nbot-1]);
      for (ntry=1;ntry<=6;ntry++) {
        PetscCall(HZStep(ntop,nbot+1,tr,dt,aa,bb,dd,uu,nn,ld,&flag));
        if (!flag) break;
        PetscCheck(ntry<6,PETSC_COMM_SELF,PETSC_ERR_CONV_FAILED,"Unable to complete hz step after six tries");
        tr = 0.9*tr; dt = 0.81*dt;
      }
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode DSSolve_GHIEP_HZ(DS ds,PetscScalar *wr,PetscScalar *wi)
{
  PetscInt          i,off;
  PetscBLASInt      n,l,n1,ld = 0;
  const PetscScalar *A,*B;
  PetscScalar       *Q;
  PetscReal         *d,*e,*s;

  PetscFunctionBegin;
#if !PetscDefined(USE_COMPLEX)
  PetscAssertPointer(wi,3);
#endif
  PetscCall(PetscBLASIntCast(ds->ld,&ld));
  PetscCall(PetscBLASIntCast(ds->n-ds->l,&n1));
  off = ds->l + ds->l*ld;
  PetscCall(DSGetArrayReal(ds,DS_MAT_T,&d));
  PetscCall(DSGetArrayReal(ds,DS_MAT_D,&s));
  e = d + ld;
#if PetscDefined(USE_DEBUG)
  /* Check signature */
  PetscCall(MatDenseGetArrayRead(ds->omat[DS_MAT_B],&B));
  for (i=0;i<ds->n;i++) {
    PetscReal de = (ds->compact)?s[i]:PetscRealPart(B[i*ld+i]);
    PetscCheck(de==1.0 || de==-1.0,PETSC_COMM_SELF,PETSC_ERR_PLIB,"Diagonal elements of the signature matrix must be 1 or -1");
  }
  PetscCall(MatDenseRestoreArrayRead(ds->omat[DS_MAT_B],&B));
#endif
  /* Quick return */
  if (n1 == 1) {
    PetscCall(MatDenseGetArray(ds->omat[DS_MAT_Q],&Q));
    for (i=0;i<=ds->l;i++) Q[i+i*ld] = 1.0;
    PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_Q],&Q));
    PetscCall(DSGHIEPComplexEigs(ds,0,ds->l,wr,wi));
    if (ds->compact) {
      wr[ds->l] = d[ds->l]/s[ds->l]; wi[ds->l] = 0.0;
    } else {
      PetscCall(MatDenseGetArrayRead(ds->omat[DS_MAT_A],&A));
      PetscCall(MatDenseGetArrayRead(ds->omat[DS_MAT_B],&B));
      d[ds->l] = PetscRealPart(A[off]);
      s[ds->l] = PetscRealPart(B[off]);
      PetscCall(MatDenseRestoreArrayRead(ds->omat[DS_MAT_A],&A));
      PetscCall(MatDenseRestoreArrayRead(ds->omat[DS_MAT_B],&B));
      wr[ds->l] = d[ds->l]/s[ds->l]; wi[ds->l] = 0.0;
    }
    PetscCall(DSRestoreArrayReal(ds,DS_MAT_T,&d));
    PetscCall(DSRestoreArrayReal(ds,DS_MAT_D,&s));
    PetscFunctionReturn(PETSC_SUCCESS);
  }
  /* Reduce to pseudotriadiagonal form */
  PetscCall(DSIntermediate_GHIEP(ds));
  PetscCall(MatDenseGetArray(ds->omat[DS_MAT_Q],&Q));
  PetscCall(PetscBLASIntCast(ds->n,&n));
  PetscCall(PetscBLASIntCast(ds->l,&l));
  PetscCall(HZIteration(n,l,d,e,s,Q,ld));
  PetscCall(MatDenseRestoreArray(ds->omat[DS_MAT_Q],&Q));
  PetscCall(DSRestoreArrayReal(ds,DS_MAT_T,&d));
  PetscCall(DSRestoreArrayReal(ds,DS_MAT_D,&s));
  if (!ds->compact) PetscCall(DSSwitchFormat_GHIEP(ds,PETSC_FALSE));
  /* Undo from diagonal the blocks with real eigenvalues*/
  PetscCall(DSGHIEPRealBlocks(ds));

  /* Recover eigenvalues from diagonal */
  PetscCall(DSGHIEPComplexEigs(ds,0,ds->n,wr,wi));
#if PetscDefined(USE_COMPLEX)
  if (wi) {
    for (i=ds->l;i<ds->n;i++) wi[i] = 0.0;
  }
#endif
  ds->t = ds->n;
  PetscFunctionReturn(PETSC_SUCCESS);
}
