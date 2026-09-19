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

import sys
import slepc4py

slepc4py.init(sys.argv)

from petsc4py import PETSc
from slepc4py import SLEPc

Print = PETSc.Sys.Print

# This function builds the discretized Laplacian operator in 2 dimensions.


def Laplacian2D(m, n):
    # Create matrix for 2D Laplacian operator
    A = PETSc.Mat().create()
    A.setSizes([m * n, m * n])
    A.setFromOptions()
    # Fill matrix
    hx = 1.0 / (m - 1)  # x grid spacing
    hy = 1.0 / (n - 1)  # y grid spacing
    diagv = 2.0 * hy / hx + 2.0 * hx / hy
    offdx = -1.0 * hy / hx
    offdy = -1.0 * hx / hy
    Istart, Iend = A.getOwnershipRange()
    for i in range(Istart, Iend):
        A[i, i] = diagv
        gi = i // n  # map row number to
        gj = i - gi * n  # grid coordinates
        if gi > 0:
            j = i - n
            A[i, j] = offdx
        if gi < m - 1:
            j = i + n
            A[i, j] = offdx
        if gj > 0:
            j = i - 1
            A[i, j] = offdy
        if gj < n - 1:
            j = i + 1
            A[i, j] = offdy
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
    for i in range(Istart, Iend):
        B[i, i] = 2.0
    if Istart == 0:
        B[0, 0] = 6.0
        B[0, 1] = -1.0
        B[1, 0] = -1.0
        B[1, 1] = 1.0
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
    ksp.setPC(pc)

    F = SLEPc.ST().create()
    F.setType(F.Type.PRECOND)
    F.setKSP(ksp)
    F.setShift(0)

    # Setup the eigensolver
    E = SLEPc.EPS().create()
    E.setST(F)
    E.setOperators(A, B)
    E.setType(E.Type.LOBPCG)
    E.setDimensions(10, PETSc.DECIDE)
    E.setWhichEigenpairs(E.Which.SMALLEST_REAL)
    E.setProblemType(problem_type)
    E.setFromOptions()

    # Solve the eigensystem
    E.solve()

    Print('')
    its = E.getIterationNumber()
    Print(f'Number of iterations of the method: {its}')
    sol_type = E.getType()
    Print(f'Solution method: {sol_type}')
    nev, _ncv, _mpd = E.getDimensions()
    Print(f'Number of requested eigenvalues: {nev}')
    tol, maxit = E.getTolerances()
    Print(f'Stopping condition: tol={tol:.4g}, maxit={maxit}')
    nconv = E.getConverged()
    Print(f'Number of converged eigenpairs: {nconv}')
    if nconv > 0:
        Print('')
        Print('        k          ||Ax-kx||/||kx|| ')
        Print('----------------- ------------------')
        for i in range(nconv):
            k = E.getEigenpair(i, xr, xi)
            error = E.computeError(i)
            if k.imag != 0.0:
                Print(f' {k.real:9f}{k.imag:+9f} j  {error:12g}')
            else:
                Print(f' {k.real:12f}       {error:12g}')
        Print('')


# The main program simply processes three user-defined command-line options
# and calls the other functions.


def main():
    opts = PETSc.Options()
    N = opts.getInt('N', 10)
    m = opts.getInt('m', N)
    n = opts.getInt('n', m)
    Print(f'Symmetric-definite Eigenproblem, N={m * n} ({m}x{n} grid)')
    A = Laplacian2D(m, n)
    B = QuasiDiagonal(m * n)
    solve_eigensystem(A, B)


if __name__ == '__main__':
    main()
