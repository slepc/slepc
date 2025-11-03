# ------------------------------------------------------------------------
#   Solves a Lypunov equation with the shifted 2-D Laplacian
# ------------------------------------------------------------------------

import sys, slepc4py
slepc4py.init(sys.argv)

from petsc4py import PETSc
from slepc4py import SLEPc

Print = PETSc.Sys.Print


def Laplacian2D(m, n, sigma):
    """
    Builds shifted negative discretized Laplacian operator in 2 dimensions.
    """
    # Create matrix for 2D Laplacian operator
    A = PETSc.Mat().create()
    A.setSizes([m*n, m*n])
    A.setFromOptions()
    # Fill matrix
    Istart, Iend = A.getOwnershipRange()
    for I in range(Istart, Iend):
        A[I,I] = -4.0-sigma
        i = I//n    # map row number to
        j = I - i*n # grid coordinates
        if i> 0  : J = I-n; A[I,J] = 1.0
        if i< m-1: J = I+n; A[I,J] = 1.0
        if j> 0  : J = I-1; A[I,J] = 1.0
        if j< n-1: J = I+1; A[I,J] = 1.0
    A.assemble()
    return A

def solve_lyap(A, C, rk):
    # Setup the solver
    L = SLEPc.LME().create()
    L.setProblemType(SLEPc.LME.ProblemType.LYAPUNOV)
    L.setCoefficients(A)
    L.setRHS(C)

    N = C.size[0]
    if rk>0:
        X1 = PETSc.Mat().createDense((N,rk))
        X1.assemble()
        X = PETSc.Mat().createLRC(None, X1, None, None)
        L.setSolution(X)

    L.setTolerances(1e-7)
    L.setFromOptions()

    # Solve the problem
    L.solve()

    if rk==0:
        X = L.getSolution()

    its = L.getIterationNumber()
    Print(f'Number of iterations of the method: {its}')
    sol_type = L.getType()
    Print(f'Solution method: {sol_type}')
    tol, maxit = L.getTolerances()
    Print(f'Stopping condition: tol={tol}, maxit={maxit}')
    Print(f'Error estimate reported by the solver: {L.getErrorEstimate()}')
    if N<500:
        Print(f'Computed residual norm: {L.computeError()}')
    Print('')

    return X

if __name__ == '__main__':
    opts = PETSc.Options()
    n = opts.getInt('n', 15)
    m = opts.getInt('m', n)
    N = m*n
    sigma = opts.getScalar('sigma', 0.0)
    rk = opts.getInt('rank', 0)

    # Coefficient matrix A
    A = Laplacian2D(m, n, sigma)

    # Create a low-rank Mat to store the right-hand side C = C1*C1'
    C1 = PETSc.Mat().createDense((N,2))
    rstart, rend = C1.getOwnershipRange()
    v = C1.getDenseArray()
    for i in range(rstart,rend):
        if i<N//2: v[i-rstart,0] = 1.0
        if i==0: v[i-rstart,1] = -2.0
        if i==1 or i==2: v[i-rstart,1] = -1.0
    C1.assemble()
    C = PETSc.Mat().createLRC(None, C1, None, None)

    X = solve_lyap(A, C, rk)

