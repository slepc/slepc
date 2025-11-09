# ex3.py: Matrix-free eigenproblem for the 2-D Laplacian
# ======================================================
#
# This example solves the eigenproblem for the 2-D discrete Laplacian
# without building the matrix explicitly.
#
# The full source code for this demo can be `downloaded here
# <../_static/ex3.py>`__.

# Initialization is similar to previous examples.

import sys, slepc4py
slepc4py.init(sys.argv)

from petsc4py import PETSc
from slepc4py import SLEPc
import numpy as np

# In this case, the program cannot be run in parallel, so we check that
# the number of MPI processes is 1. In order to enable parallelism, we
# should implement a parallel matrix-vector operation ourselves, which
# is not done in this example.

assert PETSc.COMM_WORLD.getSize() == 1

Print = PETSc.Sys.Print

# This function computes the matrix-vector product f = L*x where the
# Laplacian L is not built explicitly, and the vectors x,f are viewed
# as two-dimensional arrays associated to grid points.

def laplace2d(U, x, f):
    U[:,:] = 0
    U[1:-1, 1:-1] = x
    # Grid spacing
    m, n = x.shape
    hx = 1.0/(m-1) # x grid spacing
    hy = 1.0/(n-1) # y grid spacing
    # Setup 5-points stencil
    u  = U[1:-1, 1:-1] # center
    uN = U[1:-1,  :-2] # north
    uS = U[1:-1, 2:  ] # south
    uW = U[ :-2, 1:-1] # west
    uE = U[2:,   1:-1] # east
    # Apply Laplacian
    f[:,:] = \
         (2*u - uE - uW) * (hy/hx) \
       + (2*u - uN - uS) * (hx/hy) \

# For a matrix-free solution in slepc4py we have to create a class that
# wraps the matrix-vector operation and optionally other operations of
# the matrix. In this case, we provide the constructor and the ``mult``
# operation, that simply calls the ``laplace2d`` function above.

class Laplacian2D(object):

    def __init__(self, m, n):
        self.m, self.n = m, n
        scalar = PETSc.ScalarType
        self.U = np.zeros([m+2, n+2], dtype=scalar)

    def mult(self, A, x, y):
        m, n = self.m, self.n
        xx = x.getArray(readonly=1).reshape(m,n)
        yy = y.getArray(readonly=0).reshape(m,n)
        laplace2d(self.U, xx, yy)

# In this example, building the matrix amounts to creating an object of
# the class defined above, and passing it to a special petsc4py matrix
# with `createPython() <petsc4py.PETSc.Mat.createPython>`.

def construct_operator(m, n):
    # Create shell matrix
    context = Laplacian2D(m,n)
    A = PETSc.Mat().createPython([m*n,m*n], context)
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
    Print("Symmetric Eigenproblem (matrix-free), "
          "N=%d (%dx%d grid)" % (m*n, m, n))
    A = construct_operator(m,n)
    solve_eigensystem(A)

if __name__ == '__main__':
    main()
