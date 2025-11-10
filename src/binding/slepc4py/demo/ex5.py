# ex5.py: Simple quadratic eigenvalue problem
# ===========================================
#
# This example solves a polynomial eigenvalue problem of degree 2,
# which is the most commonly found in applications. Larger degree
# polynomials are handled similarly. The coefficient matrices in
# this example do not come from an application, they are just simple
# matrices that are easy to build, just to illustrate how it works.
# The eigenvalues in this example are purely imaginary and come in
# conjugate pairs.
#
# The full source code for this demo can be `downloaded here
# <../_static/ex5.py>`__.

# Initialization is similar to previous examples.

import sys, slepc4py
slepc4py.init(sys.argv)

from petsc4py import PETSc
from slepc4py import SLEPc

Print = PETSc.Sys.Print

# A function to build the matrices. The lowest degree coefficient is
# the 2-D Laplacian, the highest degree one is the identity matrix,
# and the other matrix is set to zero (which means that this problem
# could have been solved as a linear eigenproblem).

def construct_operators(m,n):
    Print("Quadratic Eigenproblem, N=%d (%dx%d grid)"% (m*n, m, n))
    # K is the 2-D Laplacian
    K = PETSc.Mat().create()
    K.setSizes([n*m, n*m])
    K.setFromOptions()
    Istart, Iend = K.getOwnershipRange()
    for I in range(Istart,Iend):
        v = -1.0; i = I//n; j = I-i*n;
        if i>0:
            J=I-n; K[I,J] = v
        if i<m-1:
            J=I+n; K[I,J] = v
        if j>0:
            J=I-1; K[I,J] = v
        if j<n-1:
            J=I+1; K[I,J] = v
        v=4.0; K[I,I] = v
    K.assemble()
    # C is the zero matrix
    C = PETSc.Mat().create()
    C.setSizes([n*m, n*m])
    C.setFromOptions()
    C.assemble()
    # M is the identity matrix
    M = PETSc.Mat().createConstantDiagonal([n*m, n*m], 1.0)

    return M, C, K

# The polynomial eigenvalue solver is similar to the linear eigensolver
# used in previous examples. The main difference is that we must provide
# a list of matrices, from lowest to highest degree.

def solve_eigensystem(M, C, K):
    # Setup the eigensolver
    Q = SLEPc.PEP().create()
    Q.setOperators([K, C, M])
    Q.setDimensions(6)
    Q.setProblemType(SLEPc.PEP.ProblemType.GENERAL)
    Q.setFromOptions()
    # Solve the eigensystem
    Q.solve()
    # Create the result vectors
    xr, xi = K.createVecs()

    its = Q.getIterationNumber()
    Print("Number of iterations of the method: %i" % its)
    sol_type = Q.getType()
    Print("Solution method: %s" % sol_type)
    nev, ncv, mpd = Q.getDimensions()
    Print("")
    Print("Number of requested eigenvalues: %i" % nev)
    tol, maxit = Q.getTolerances()
    Print("Stopping condition: tol=%.4g, maxit=%d" % (tol, maxit))
    nconv = Q.getConverged()
    Print("Number of converged approximate eigenpairs: %d" % nconv)
    if nconv > 0:
        Print("")
        Print("          k           ||(k^2M+Ck+K)x||/||kx|| ")
        Print("-------------------- -------------------------")
        for i in range(nconv):
            k = Q.getEigenpair(i, xr, xi)
            error = Q.computeError(i)
            if k.imag != 0.0:
                Print("%9f%+9f j    %12g" % (k.real, k.imag, error))
            else:
                Print("%12f         %12g" % (k.real, error))
    Print("")

# The main program simply processes two user-defined command-line options
# (the dimensions of the mesh) and calls the other two functions.

if __name__ == '__main__':
    opts = PETSc.Options()
    m = opts.getInt('m', 32)
    n = opts.getInt('n', m)
    M, C, K = construct_operators(m,n)
    solve_eigensystem(M, C, K)
    M = C = K = None
