# ex2.py: Standard symmetric eigenproblem for the 2-D Laplacian
# =============================================================
#
# This example computes eigenvalues and eigenvectors of the discrete Laplacian
# on a two-dimensional domain with finite differences.
#
# The full source code for this demo can be `downloaded here
# <../_static/ex2.py>`__.

# Initialization is similar to previous examples.

try: range = xrange
except: pass

import sys, slepc4py
slepc4py.init(sys.argv)

from petsc4py import PETSc
from slepc4py import SLEPc

Print = PETSc.Sys.Print

# In this example we have organized the code in several functions. This
# one builds the finite-difference Laplacian matrix by computing the
# indices of each entry. An alternative would be to use the functionality
# offered by `DMDA <petsc4py.PETSc.DMDA>`.

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

# This function receives the matrix and the problem type, then solves the
# eigenvalue problem and prints information about the computed solution.
# Although we know that eigenvalues and eigenvectors are real in this
# example, the function is prepared to solve it as a non-symmetric problem,
# by passing `SLEPc.EPS.ProblemType.NHEP`, that is why the code handles
# possibly complex eigenvalues and eigenvectors.

def solve_eigensystem(A, problem_type=SLEPc.EPS.ProblemType.HEP):
    # Create the result vectors
    xr, xi = A.createVecs()

    # Setup the eigensolver
    E = SLEPc.EPS().create()
    E.setOperators(A,None)
    E.setDimensions(3,PETSc.DECIDE)
    E.setProblemType(problem_type)
    E.setFromOptions()

    # Solve the eigensystem
    E.solve()

    Print("")
    its = E.getIterationNumber()
    Print("Number of iterations of the method: %i" % its)
    sol_type = E.getType()
    Print("Solution method: %s" % sol_type)
    nev, ncv, mpd = E.getDimensions()
    Print("Number of requested eigenvalues: %i" % nev)
    tol, maxit = E.getTolerances()
    Print("Stopping condition: tol=%.4g, maxit=%d" % (tol, maxit))
    nconv = E.getConverged()
    Print("Number of converged eigenpairs: %d" % nconv)
    if nconv > 0:
        Print("")
        Print("        k          ||Ax-kx||/||kx|| ")
        Print("----------------- ------------------")
        for i in range(nconv):
            k = E.getEigenpair(i, xr, xi)
            error = E.computeError(i)
            if k.imag != 0.0:
              Print(" %9f%+9f j  %12g" % (k.real, k.imag, error))
            else:
              Print(" %12f       %12g" % (k.real, error))
        Print("")

# The main program simply processes three user-defined command-line options
# and calls the other two functions.

def main():
    opts = PETSc.Options()
    N = opts.getInt('N', 32)
    m = opts.getInt('m', N)
    n = opts.getInt('n', m)
    Print("Symmetric Eigenproblem (sparse matrix), "
          "N=%d (%dx%d grid)" % (m*n, m, n))
    A = construct_operator(m,n)
    solve_eigensystem(A)

if __name__ == '__main__':
    main()
