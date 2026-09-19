# ex2.py: Standard symmetric eigenproblem for the 2-D Laplacian
# =============================================================
#
# This example computes eigenvalues and eigenvectors of the discrete Laplacian
# on a two-dimensional domain with finite differences.
#
# The full source code for this demo can be `downloaded here
# <../_static/ex2.py>`__.

# Initialization is similar to previous examples.

import sys
import slepc4py

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
    E.setOperators(A, None)
    E.setDimensions(3, PETSc.DECIDE)
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
# and calls the other two functions.


def main():
    opts = PETSc.Options()
    N = opts.getInt('N', 32)
    m = opts.getInt('m', N)
    n = opts.getInt('n', m)
    Print(f'Symmetric Eigenproblem (sparse matrix), N={m * n} ({m}x{n} grid)')
    A = construct_operator(m, n)
    solve_eigensystem(A)


if __name__ == '__main__':
    main()
