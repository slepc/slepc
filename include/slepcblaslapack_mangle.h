/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#pragma once

/* MANSEC = Sys */

/* LAPACK functions without string parameters */
#define BLASrot_     PETSCBLAS(rot,ROT)
#define BLASMIXEDrot_ PETSCBLASMIXED(rot,ROT)
#if !defined(SLEPC_MISSING_LAPACK_LAEV2)
#define LAPACKlaev2_ PETSCBLAS(laev2,LAEV2)
#endif
#if !defined(SLEPC_MISSING_LAPACK_GEHRD)
#define LAPACKgehrd_ PETSCBLAS(gehrd,GEHRD)
#endif
#if !defined(SLEPC_MISSING_LAPACK_LARFG)
#define LAPACKlarfg_ PETSCBLAS(larfg,LARFG)
#endif
#if !defined(SLEPC_MISSING_LAPACK_LAG2)
#define LAPACKlag2_  PETSCBLASREAL(lag2,LAG2)
#endif
#if !defined(SLEPC_MISSING_LAPACK_LASV2)
#define LAPACKlasv2_ PETSCBLASREAL(lasv2,LASV2)
#endif
#if !defined(SLEPC_MISSING_LAPACK_LARTG)
#define LAPACKlartg_ PETSCBLAS(lartg,LARTG)
#define LAPACKREALlartg_ PETSCBLASREAL(lartg,LARTG)
#endif
#if !defined(SLEPC_MISSING_LAPACK_LAED4)
#define LAPACKlaed4_ PETSCBLASREAL(laed4,LAED4)
#endif
#if !defined(SLEPC_MISSING_LAPACK_LAMRG)
#define LAPACKlamrg_ PETSCBLASREAL(lamrg,LAMRG)
#endif
#if !defined(SLEPC_MISSING_LAPACK_ORGHR)
#if !defined(PETSC_USE_COMPLEX)
#define LAPACKorghr_ PETSCBLAS(orghr,ORGHR)
#else
#define LAPACKorghr_ PETSCBLAS(unghr,UNGHR)
#endif
#endif
#if !defined(SLEPC_MISSING_LAPACK_TGEXC)
#define LAPACKtgexc_ PETSCBLAS(tgexc,TGEXC)
#endif
#define LAPACKgeqp3_ PETSCBLAS(geqp3,GEQP3)

/* LAPACK functions with string parameters */

/* same name for real and complex */
#define BLAStrmm_    PETSCBLAS(trmm,TRMM)
#if !defined(SLEPC_MISSING_LAPACK_LANHS)
#define LAPACKlanhs_ PETSCBLAS(lanhs,LANHS)
#endif
#define LAPACKlange_ PETSCBLAS(lange,LANGE)
#if !defined(SLEPC_MISSING_LAPACK_LARF)
#define LAPACKlarf_  PETSCBLAS(larf,LARF)
#endif
#define LAPACKlansy_ PETSCBLAS(lansy,LANSY)
#if !defined(SLEPC_MISSING_LAPACK_TRSYL)
#define LAPACKtrsyl_ PETSCBLAS(trsyl,TRSYL)
#endif

/* subroutines in which we use only the real version, do not care whether they have different name */
#if !defined(SLEPC_MISSING_LAPACK_STEVR)
#define LAPACKstevr_ PETSCBLASREAL(stevr,STEVR)
#endif
#if !defined(SLEPC_MISSING_LAPACK_BDSDC)
#define LAPACKbdsdc_ PETSCBLASREAL(bdsdc,BDSDC)
#endif
#define LAPACKlamch_ PETSCBLASREAL(lamch,LAMCH)
#define LAPACKlamc3_ PETSCBLASREAL(lamc3,LAMC3)

/* subroutines with different name in real/complex */
#if !defined(PETSC_USE_COMPLEX)
#if !defined(SLEPC_MISSING_LAPACK_ORGTR)
#define LAPACKorgtr_ PETSCBLAS(orgtr,ORGTR)
#endif
#if !defined(SLEPC_MISSING_LAPACK_ORMBR)
#define LAPACKormbr_ PETSCBLAS(ormbr,ORMBR)
#endif
#if !defined(SLEPC_MISSING_LAPACK_SYTRD)
#define LAPACKsytrd_ PETSCBLAS(sytrd,SYTRD)
#endif
#define LAPACKsyevd_ PETSCBLAS(syevd,SYEVD)
#define LAPACKsygvd_ PETSCBLAS(sygvd,SYGVD)
#else
#if !defined(SLEPC_MISSING_LAPACK_ORGTR)
#define LAPACKorgtr_ PETSCBLAS(ungtr,UNGTR)
#endif
#if !defined(SLEPC_MISSING_LAPACK_ORMBR)
#define LAPACKormbr_ PETSCBLAS(unmbr,UNMBR)
#endif
#if !defined(SLEPC_MISSING_LAPACK_SYTRD)
#define LAPACKsytrd_ PETSCBLAS(hetrd,HETRD)
#endif
#define LAPACKsyevd_ PETSCBLAS(heevd,HEEVD)
#define LAPACKsygvd_ PETSCBLAS(hegvd,HEGVD)
#endif

/* subroutines with different signature in real/complex */
#define LAPACKggev_  PETSCBLAS(ggev,GGEV)
#if !defined(SLEPC_MISSING_LAPACK_GGEV3)   /* ggev and ggev3 have the same signature */
#define LAPACKggev3_ PETSCBLAS(ggev3,GGEV3)
#define LAPACKggevalt_ LAPACKggev3_
#else
#define LAPACKggevalt_ LAPACKggev_
#endif
#if !defined(SLEPC_MISSING_LAPACK_GGES3)   /* gges and gges3 have the same signature */
#define LAPACKgges3_ PETSCBLAS(gges3,GGES3)
#define LAPACKggesalt_ LAPACKgges3_
#else
#define LAPACKggesalt_ LAPACKgges_
#endif
#if !defined(SLEPC_MISSING_LAPACK_GGSVD)
#define LAPACKggsvd_ PETSCBLAS(ggsvd,GGSVD)
#endif
#if !defined(SLEPC_MISSING_LAPACK_GGSVD3)
#define LAPACKggsvd3_ PETSCBLAS(ggsvd3,GGSVD3)
#endif
#if !defined(SLEPC_MISSING_LAPACK_TREVC)
#define LAPACKtrevc_ PETSCBLAS(trevc,TREVC)
#endif
#define LAPACKgeevx_ PETSCBLAS(geevx,GEEVX)
#define LAPACKgees_  PETSCBLAS(gees,GEES)
#if !defined(SLEPC_MISSING_LAPACK_TREXC)
#define LAPACKtrexc_ PETSCBLAS(trexc,TREXC)
#endif
#define LAPACKgesdd_ PETSCBLAS(gesdd,GESDD)
#if !defined(SLEPC_MISSING_LAPACK_TGEVC)
#define LAPACKtgevc_ PETSCBLAS(tgevc,TGEVC)
#endif
#if !defined(SLEPC_MISSING_LAPACK_HSEIN)
#define LAPACKhsein_ PETSCBLAS(hsein,HSEIN)
#endif
#if !defined(SLEPC_MISSING_LAPACK_STEDC)
#define LAPACKstedc_ PETSCBLAS(stedc,STEDC)
#endif
#if !defined(SLEPC_MISSING_LAPACK_LASCL)
#define LAPACKlascl_ PETSCBLAS(lascl,LASCL)
#endif
#if !defined(SLEPC_MISSING_LAPACK_BDSQR)
#define LAPACKbdsqr_ PETSCBLAS(bdsqr,BDSQR)
#endif
#if !defined(SLEPC_MISSING_LAPACK_GEBRD)
#define LAPACKgebrd_ PETSCBLAS(gebrd,GEBRD)
#endif

#if defined(PETSC_HAVE_COMPLEX)
/* complex subroutines to be called with scalar-type=real */
#define BLASCOMPLEXgemm_   PETSCBLASCOMPLEX(gemm,GEMM)
#define BLASCOMPLEXscal_   PETSCBLASCOMPLEX(scal,SCAL)
#define LAPACKCOMPLEXgesv_ PETSCBLASCOMPLEX(gesv,GESV)
#endif
