# ex9.py: Generalized symmetric-definite eigenproblem
# ===================================================
#
# This example computes eigenvalues and eigenvectors of a generalized
# symmetric-definite eigenvalue problem, where the first matrix is the
# discrete Laplacian in two dimensions and the second matrix is quasi
# diagonal.
#
# The full source code for this demo can be `downloaded here
# <../_static/ex9.py>`__.

# Initialization is similar to previous examples.

try: range = xrange
except: pass

import sys, slepc4py
slepc4py.init(sys.argv)

from petsc4py import PETSc
from slepc4py import SLEPc

Print = PETSc.Sys.Print

# This function builds the discretized Laplacian operator in 2 dimensions.

def Laplacian2D(m, n):
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
    for I in range(Istart, Iend):
        A[I,I] = diagv
        i = I//n    # map row number to
        j = I - i*n # grid coordinates
        if i> 0  : J = I-n; A[I,J] = offdx
        if i< m-1: J = I+n; A[I,J] = offdx
        if j> 0  : J = I-1; A[I,J] = offdy
        if j< n-1: J = I+1; A[I,J] = offdy
    A.assemble()
    return A

# This function builds a quasi-diagonal matrix. It is two times the identity
# matrix except for the 2x2 leading submatrix ``[6 -1; -1 1]``.

def QuasiDiagonal(N):
    # Create matrix
    B = PETSc.Mat().create()
    B.setSizes([N, N])
    B.setFromOptions()
    # Fill matrix
    Istart, Iend = B.getOwnershipRange()
    for I in range(Istart, Iend):
        B[I,I] = 2.0
    if Istart==0:
        B[0,0] = 6.0
        B[0,1] = -1.0
        B[1,0] = -1.0
        B[1,1] = 1.0
    B.assemble()
    return B

# The following function receives the two matrices and solves the
# eigenproblem. In this example we illustrate how to pass objects
# that have been created beforehand, instead of extracting the internal
# objects. We are using a spectral transformation of type `ST.Type.PRECOND`
# and a Block Jacobi preconditioner. We want to compute the leftmost
# eigenvalues. The selected eigensolver is LOBPCG, which is appropriate
# for this use case. After the solve, we print the computed solution.

def solve_eigensystem(A, B, problem_type=SLEPc.EPS.ProblemType.GHEP):
    # Create the results vectors
    xr, xi = A.createVecs()

    pc = PETSc.PC().create()
    # pc.setType(pc.Type.HYPRE)
    pc.setType(pc.Type.BJACOBI)

    ksp = PETSc.KSP().create()
    ksp.setType(ksp.Type.PREONLY)
    ksp.setPC( pc )

    F = SLEPc.ST().create()
    F.setType(F.Type.PRECOND)
    F.setKSP( ksp )
    F.setShift(0)

    # Setup the eigensolver
    E = SLEPc.EPS().create()
    E.setST(F)
    E.setOperators(A,B)
    E.setType(E.Type.LOBPCG)
    E.setDimensions(10,PETSc.DECIDE)
    E.setWhichEigenpairs(E.Which.SMALLEST_REAL)
    E.setProblemType( problem_type )
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
# and calls the other functions.

def main():
    opts = PETSc.Options()
    N = opts.getInt('N', 10)
    m = opts.getInt('m', N)
    n = opts.getInt('n', m)
    Print("Symmetric-definite Eigenproblem, N=%d (%dx%d grid)" % (m*n, m, n))
    A = Laplacian2D(m,n)
    B = QuasiDiagonal(m*n)
    solve_eigensystem(A,B)

if __name__ == '__main__':
    main()
