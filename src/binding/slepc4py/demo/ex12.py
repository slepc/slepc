# ex12.py: Illustrate the use of arbitrary selection
# ==================================================
#
# This example solves a simple tridiagonal eigenproblem. It illustrates
# how to set up the arbitrary selection of eigenvalues, where the
# decision of which is the preferred eigenvalue is made based not only
# on the value of the approximate eigenvalue but also on the approximate
# eigenvector.
#
# In this example, the selection criterion is based on the projection
# of the approximate eigenvector onto a precomputed eigenvector. That is
# why we solve the problem twice.
#
# The full source code for this demo can be `downloaded here
# <../_static/ex12.py>`__.

# Initialization is similar to previous examples.

import sys, slepc4py
slepc4py.init(sys.argv)

from petsc4py import PETSc
from slepc4py import SLEPc
import numpy

# The matrix size ``n`` can be specified at the command line.

opts = PETSc.Options()
n = opts.getInt('n', 30)

# Create the matrix ``tridiag([-1 0 -1])``.

A = PETSc.Mat(); A.create()
A.setSizes([n, n])
A.setFromOptions()
rstart, rend = A.getOwnershipRange()
for i in range(rstart, rend):
    if i>0: A[i, i-1] = -1
    if i<n-1: A[i, i+1] = -1
A.assemble()

# Configure the linear eigensolver initially to compute leftmost
# eigenvalues.

E = SLEPc.EPS(); E.create()
E.setOperators(A)
E.setProblemType(SLEPc.EPS.ProblemType.HEP)
E.setWhichEigenpairs(SLEPc.EPS.Which.SMALLEST_REAL)
E.setFromOptions()

# Solve the eigenproblem and store the first computed eigenvector in
# ``sx`` to be used later. For the second solve, we configure the
# solver to select largest magnitude values with an arbitrary
# selection callback function ``myArbitrarySel()``. It means that instead
# of sorting eigenvalues, the solver will sort the approximations
# according to the largest values of the result of ``myArbitrarySel()``
# evaluated on the approximate eigenvectors. In this way, the same
# eigenvalue should be computed again.

E.solve()
nconv = E.getConverged()
Print = PETSc.Sys.Print
vw = PETSc.Viewer.STDOUT()
if nconv>0:
    sx, _ = A.createVecs()
    E.getEigenpair(0, sx)
    vw.pushFormat(PETSc.Viewer.Format.ASCII_INFO_DETAIL)
    E.errorView(viewer=vw)
    def myArbitrarySel(evalue, xr, xi, sx):
        return abs(xr.dot(sx))
    E.setArbitrarySelection(myArbitrarySel,sx)
    E.setWhichEigenpairs(SLEPc.EPS.Which.LARGEST_MAGNITUDE)
    E.solve()
    E.errorView(viewer=vw)
    vw.popFormat()
else:
    Print( "No eigenpairs converged" )
