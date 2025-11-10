# ex11.py: 2-D Laplacian eigenproblem solved with contour integral
# ================================================================
#
# This example is similar to ``ex2.py``, but employs a contour integral
# solver. It illustrates how to define a region of the complex plane
# using an `RG` object.
#
# The full source code for this demo can be `downloaded here
# <../_static/ex11.py>`__.

# Initialization is similar to previous examples.

try: range = xrange
except: pass

import sys, slepc4py
slepc4py.init(sys.argv)

from petsc4py import PETSc
from slepc4py import SLEPc

Print = PETSc.Sys.Print

# Build the finite-difference 2-D Laplacian matrix.

def construct_operator(m, n):
    # Create matrix for 2D Laplacian operator
    A = PETSc.Mat().create()
    A.setSizes([m*n, m*n])
    A.setFromOptions()
    # Fill matrix
    hx = 1.0/(m-1) # x grid spacing
    hy = 1.0/(n-1) # y grid spacing
    diagv = 2.0*hy/hx + 2.0*hx/hy
    offdx = -1.0*hy/hx
    offdy = -1.0*hx/hy
    Istart, Iend = A.getOwnershipRange()
    for I in range(Istart, Iend) :
        A[I,I] = diagv
        i = I//n    # map row number to
        j = I - i*n # grid coordinates
        if i> 0  : J = I-n; A[I,J] = offdx
        if i< m-1: J = I+n; A[I,J] = offdx
        if j> 0  : J = I-1; A[I,J] = offdy
        if j< n-1: J = I+1; A[I,J] = offdy
    A.assemble()
    return A

# In the main function, first two command-line options are processed to
# set the grid dimensions. Then the matrix is built and passed to the
# solver object. In this case, the solver is configured to use the contour
# integral method. Next, the region of interest is defined, in this case
# an ellipse centered at the origin, with radius 0.2 and vertical scaling
# of 0.1. Finally, the solver is run. In this example, we illustrate how to
# print the solution using the solver method `errorView() <EPS.errorView()>`.

def main():
    opts = PETSc.Options()
    n = opts.getInt('n', 32)
    m = opts.getInt('m', 32)
    Print("2-D Laplacian Eigenproblem solved with contour integral, "
          "N=%d (%dx%d grid)\n" % (m*n, m, n))
    A = construct_operator(m,n)

    E = SLEPc.EPS().create()
    E.setOperators(A)
    E.setProblemType(SLEPc.EPS.ProblemType.HEP)
    E.setType(SLEPc.EPS.Type.CISS)

    R = E.getRG()
    R.setType(SLEPc.RG.Type.ELLIPSE)
    R.setEllipseParameters(0.0,0.2,0.1)
    E.setFromOptions()

    E.solve()

    vw = PETSc.Viewer.STDOUT()
    vw.pushFormat(PETSc.Viewer.Format.ASCII_INFO_DETAIL)
    E.errorView(viewer=vw)
    vw.popFormat()

if __name__ == '__main__':
    main()
