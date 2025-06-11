/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* @@@ BLOPEX (version 1.1) LGPL Version 2.1 or above.See www.gnu.org. */
/* @@@ Copyright 2010 BLOPEX team https://github.com/lobpcg/blopex     */
/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

#pragma once

#include <interpreter.h>

#if !defined(PETSC_USE_COMPLEX)
SLEPC_INTERN BlopexInt PETSC_dpotrf_interface(char*,BlopexInt*,double*,BlopexInt*,BlopexInt*);
SLEPC_INTERN BlopexInt PETSC_dsygv_interface(BlopexInt*,char*,char*,BlopexInt*,double*,BlopexInt*,double*,BlopexInt*,double*,double*,BlopexInt*,BlopexInt*);
#else
SLEPC_INTERN BlopexInt PETSC_zpotrf_interface(char*,BlopexInt*,komplex*,BlopexInt*,BlopexInt*);
SLEPC_INTERN BlopexInt PETSC_zsygv_interface(BlopexInt*,char*,char*,BlopexInt*,komplex*,BlopexInt*,komplex*,BlopexInt*,double*,komplex*,BlopexInt*,double*,BlopexInt*);
#endif

SLEPC_INTERN int PETSCSetupInterpreter(mv_InterfaceInterpreter*);
