/* ---------------------------------------------------------------- */

%include petsc4py/petsc4py.i

/* ---------------------------------------------------------------- */

%header %{#include "slepc4py/slepc4py.h"%}
%init   %{import_slepc4py();%}

%define SWIG_TYPECHECK_SLEPC_ST  650 %enddef
%define SWIG_TYPECHECK_SLEPC_BV  651 %enddef
%define SWIG_TYPECHECK_SLEPC_DS  652 %enddef
%define SWIG_TYPECHECK_SLEPC_FN  653 %enddef
%define SWIG_TYPECHECK_SLEPC_RG  654 %enddef
%define SWIG_TYPECHECK_SLEPC_EPS 655 %enddef
%define SWIG_TYPECHECK_SLEPC_SVD 656 %enddef
%define SWIG_TYPECHECK_SLEPC_PEP 657 %enddef
%define SWIG_TYPECHECK_SLEPC_NEP 658 %enddef
%define SWIG_TYPECHECK_SLEPC_MFN 659 %enddef
%define SWIG_TYPECHECK_SLEPC_LME 660 %enddef

%define %slepc4py_objt(Pkg, PyType, Type, CODE)
%petsc4py_objt(Pkg, PyType, Type, CODE)
%enddef /* %slepc4py_objt */

/* ---------------------------------------------------------------- */

%slepc4py_objt( Slepc , ST ,  ST ,  SLEPC_ST  )
%slepc4py_objt( Slepc , BV ,  BV ,  SLEPC_BV  )
%slepc4py_objt( Slepc , DS ,  DS ,  SLEPC_DS  )
%slepc4py_objt( Slepc , FN ,  FN ,  SLEPC_FN  )
%slepc4py_objt( Slepc , RG ,  RG ,  SLEPC_RG  )
%slepc4py_objt( Slepc , EPS , EPS , SLEPC_EPS )
%slepc4py_objt( Slepc , SVD , SVD , SLEPC_SVD )
%slepc4py_objt( Slepc , PEP , PEP , SLEPC_PEP )
%slepc4py_objt( Slepc , NEP , NEP , SLEPC_NEP )
%slepc4py_objt( Slepc , MFN , MFN , SLEPC_MFN )
%slepc4py_objt( Slepc , LME , LME , SLEPC_LME )

/* ---------------------------------------------------------------- */

/*
 * Local Variables:
 * mode: C
 * End:
 */
