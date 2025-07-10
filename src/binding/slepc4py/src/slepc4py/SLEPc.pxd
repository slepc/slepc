# -----------------------------------------------------------------------------

cdef extern from "<slepc.h>":

    struct _p_ST
    ctypedef _p_ST* SlepcST "ST"

    struct _p_BV
    ctypedef _p_BV* SlepcBV "BV"

    struct _p_DS
    ctypedef _p_DS* SlepcDS "DS"

    struct _p_FN
    ctypedef _p_FN* SlepcFN "FN"

    struct _p_RG
    ctypedef _p_RG* SlepcRG "RG"

    struct _p_EPS
    ctypedef _p_EPS* SlepcEPS "EPS"

    struct _p_SVD
    ctypedef _p_SVD* SlepcSVD "SVD"

    struct _p_PEP
    ctypedef _p_PEP* SlepcPEP "PEP"

    struct _p_NEP
    ctypedef _p_NEP* SlepcNEP "NEP"

    struct _p_MFN
    ctypedef _p_MFN* SlepcMFN "MFN"

    struct _p_LME
    ctypedef _p_LME* SlepcLME "LME"

# -----------------------------------------------------------------------------

from petsc4py.PETSc cimport Object

ctypedef public api class ST(Object) [
    type   PySlepcST_Type,
    object PySlepcSTObject,
    ]:
    cdef SlepcST st

ctypedef public api class BV(Object) [
    type   PySlepcBV_Type,
    object PySlepcBVObject,
    ]:
    cdef SlepcBV bv

ctypedef public api class DS(Object) [
    type   PySlepcDS_Type,
    object PySlepcDSObject,
    ]:
    cdef SlepcDS ds

ctypedef public api class FN(Object) [
    type   PySlepcFN_Type,
    object PySlepcFNObject,
    ]:
    cdef SlepcFN fn

ctypedef public api class RG(Object) [
    type   PySlepcRG_Type,
    object PySlepcRGObject,
    ]:
    cdef SlepcRG rg

ctypedef public api class EPS(Object) [
    type PySlepcEPS_Type,
    object PySlepcEPSObject,
    ]:
    cdef SlepcEPS eps

ctypedef public api class SVD(Object) [
    type   PySlepcSVD_Type,
    object PySlepcSVDObject,
    ]:
    cdef SlepcSVD svd

ctypedef public api class PEP(Object) [
    type   PySlepcPEP_Type,
    object PySlepcPEPObject,
    ]:
    cdef SlepcPEP pep

ctypedef public api class NEP(Object) [
    type   PySlepcNEP_Type,
    object PySlepcNEPObject,
    ]:
    cdef SlepcNEP nep

ctypedef public api class MFN(Object) [
    type   PySlepcMFN_Type,
    object PySlepcMFNObject,
    ]:
    cdef SlepcMFN mfn

ctypedef public api class LME(Object) [
    type   PySlepcLME_Type,
    object PySlepcLMEObject,
    ]:
    cdef SlepcLME lme

# -----------------------------------------------------------------------------
