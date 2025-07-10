# -----------------------------------------------------------------------------

cdef inline int setref(void *d, void *s) except -1:
    cdef PetscObject *dest  = <PetscObject*> d
    cdef PetscObject source = <PetscObject>  s
    CHKERR( PetscINCREF(&source) )
    dest[0] = source
    return 0

# -----------------------------------------------------------------------------

# -- ST --

cdef api object PySlepcST_New(SlepcST arg):
    cdef ST retv = ST()
    setref(&retv.st, arg)
    return retv

cdef api SlepcST PySlepcST_Get(object arg) except ? NULL:
    cdef SlepcST retv = NULL
    cdef ST ob = <ST?> arg
    retv = ob.st
    return retv

# -----------------------------------------------------------------------------

# -- BV --

cdef api object PySlepcBV_New(SlepcBV arg):
    cdef BV retv = BV()
    setref(&retv.bv, arg)
    return retv

cdef api SlepcBV PySlepcBV_Get(object arg) except ? NULL:
    cdef SlepcBV retv = NULL
    cdef BV ob = <BV?> arg
    retv = ob.bv
    return retv

# -----------------------------------------------------------------------------

# -- DS --

cdef api object PySlepcDS_New(SlepcDS arg):
    cdef DS retv = DS()
    setref(&retv.ds, arg)
    return retv

cdef api SlepcDS PySlepcDS_Get(object arg) except ? NULL:
    cdef SlepcDS retv = NULL
    cdef DS ob = <DS?> arg
    retv = ob.ds
    return retv

# -----------------------------------------------------------------------------

# -- FN --

cdef api object PySlepcFN_New(SlepcFN arg):
    cdef FN retv = FN()
    setref(&retv.fn, arg)
    return retv

cdef api SlepcFN PySlepcFN_Get(object arg) except ? NULL:
    cdef SlepcFN retv = NULL
    cdef FN ob = <FN?> arg
    retv = ob.fn
    return retv

# -----------------------------------------------------------------------------

# -- RG --

cdef api object PySlepcRG_New(SlepcRG arg):
    cdef RG retv = RG()
    setref(&retv.rg, arg)
    return retv

cdef api SlepcRG PySlepcRG_Get(object arg) except ? NULL:
    cdef SlepcRG retv = NULL
    cdef RG ob = <RG?> arg
    retv = ob.rg
    return retv

# -----------------------------------------------------------------------------

# -- EPS --

cdef api object PySlepcEPS_New(SlepcEPS arg):
    cdef EPS retv = EPS()
    setref(&retv.eps, arg)
    return retv

cdef api SlepcEPS PySlepcEPS_Get(object arg) except ? NULL:
    cdef SlepcEPS retv = NULL
    cdef EPS ob = <EPS?> arg
    retv = ob.eps
    return retv

# -----------------------------------------------------------------------------

# -- SVD --

cdef api object PySlepcSVD_New(SlepcSVD arg):
    cdef SVD retv = SVD()
    setref(&retv.svd, arg)
    return retv

cdef api SlepcSVD PySlepcSVD_Get(object arg) except ? NULL:
    cdef SlepcSVD retv = NULL
    cdef SVD ob = <SVD?> arg
    retv = ob.svd
    return retv

# -----------------------------------------------------------------------------

# -- PEP --

cdef api object PySlepcPEP_New(SlepcPEP arg):
    cdef PEP retv = PEP()
    setref(&retv.pep, arg)
    return retv

cdef api SlepcPEP PySlepcPEP_Get(object arg) except ? NULL:
    cdef SlepcPEP retv = NULL
    cdef PEP ob = <PEP?> arg
    retv = ob.pep
    return retv

# -----------------------------------------------------------------------------

# -- NEP --

cdef api object PySlepcNEP_New(SlepcNEP arg):
    cdef NEP retv = NEP()
    setref(&retv.nep, arg)
    return retv

cdef api SlepcNEP PySlepcNEP_Get(object arg) except ? NULL:
    cdef SlepcNEP retv = NULL
    cdef NEP ob = <NEP?> arg
    retv = ob.nep
    return retv

# -----------------------------------------------------------------------------

# -- MFN --

cdef api object PySlepcMFN_New(SlepcMFN arg):
    cdef MFN retv = MFN()
    setref(&retv.mfn, arg)
    return retv

cdef api SlepcMFN PySlepcMFN_Get(object arg) except ? NULL:
    cdef SlepcMFN retv = NULL
    cdef MFN ob = <MFN?> arg
    retv = ob.mfn
    return retv

# -----------------------------------------------------------------------------

# -- LME --

cdef api object PySlepcLME_New(SlepcLME arg):
    cdef LME retv = LME()
    setref(&retv.lme, arg)
    return retv

cdef api SlepcLME PySlepcLME_Get(object arg) except ? NULL:
    cdef SlepcLME retv = NULL
    cdef LME ob = <LME?> arg
    retv = ob.lme
    return retv

# -----------------------------------------------------------------------------
