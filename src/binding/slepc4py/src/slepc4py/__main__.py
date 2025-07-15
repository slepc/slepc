"""
Command line access to the SLEPc Options Database.

This module provides command line access to SLEPc Options
Database. It outputs a listing of the many SLEPc options
indicating option names, default values and descriptions.
Usage::

  $ python -m slepc4py [eps|svd|pep|nep|mfn|st|bv|rg|fn|ds] [<slepc-option-list>]

"""

def help(args=None):
    import sys
    # program name
    try:
        prog = sys.argv[0]
    except Exception:
        prog = getattr(sys, 'executable', 'python')
    # arguments
    if args is None:
        args = sys.argv[1:]
    elif isinstance(args, str):
        args = args.split()
    else:
        args = [str(a) for a in args]
    if '-help' not in args:
        args.append('-help')
    # initialization
    import slepc4py
    slepc4py.init([prog] + args)
    from slepc4py import SLEPc
    # and finally ...
    COMM = SLEPc.COMM_SELF
    if 'eps' in args:
        eps = SLEPc.EPS().create(comm=COMM)
        eps.setFromOptions()
        eps.destroy()
        del eps
    if 'svd' in args:
        svd = SLEPc.SVD().create(comm=COMM)
        svd.setFromOptions()
        svd.destroy()
        del svd
    if 'pep' in args:
        pep = SLEPc.PEP().create(comm=COMM)
        pep.setFromOptions()
        pep.destroy()
        del pep
    if 'nep' in args:
        nep = SLEPc.NEP().create(comm=COMM)
        nep.setFromOptions()
        nep.destroy()
        del nep
    if 'mfn' in args:
        mfn = SLEPc.MFN().create(comm=COMM)
        mfn.setFromOptions()
        mfn.destroy()
        del mfn
    if 'lme' in args:
        lme = SLEPc.LME().create(comm=COMM)
        lme.setFromOptions()
        lme.destroy()
        del lme
    if 'st' in args:
        st = SLEPc.ST().create(comm=COMM)
        st.setFromOptions()
        st.destroy()
        del st
    if 'bv' in args:
        bv = SLEPc.BV().create(comm=COMM)
        bv.setFromOptions()
        bv.destroy()
        del bv
    if 'rg' in args:
        rg = SLEPc.RG().create(comm=COMM)
        rg.setFromOptions()
        rg.destroy()
        del rg
    if 'fn' in args:
        fn = SLEPc.FN().create(comm=COMM)
        fn.setFromOptions()
        fn.destroy()
        del fn
    if 'ds' in args:
        ds = SLEPc.DS().create(comm=COMM)
        ds.setFromOptions()
        ds.destroy()
        del ds

if __name__ == '__main__':
    help()
