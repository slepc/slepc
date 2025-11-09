# ex6.py: Compute exp(t*A)*v for a matrix from a Markov model
# ===========================================================
#
# This example illustrates the functionality in slepc4py for computing
# matrix functions, or more precisely, the application of a matrix
# function on a given vector. The example works with the exponential
# function, which is most commonly found in applications.
#
# The main focus of slepc4py is eigenvalue and singular value problems,
# but it has some codes to deal with matrix functions, which sometimes
# are needed in the context of eigenproblems, but have interest on
# their own.
#
# The full source code for this demo can be `downloaded here
# <../_static/ex6.py>`__.

# Initialization is similar to previous examples.

import sys, slepc4py
slepc4py.init(sys.argv)

from petsc4py import PETSc
from slepc4py import SLEPc

Print = PETSc.Sys.Print

# This function builds a matrix that implements a Markov model of a random
# walk on a triangular grid. The entries of the matrix represent
# probabilities of moving to neighboring cells in the grid.

def build_matrix(m):
    N = m*(m+1)/2
    Print("Markov y=exp(t*A)*e_1, N=%d (m=%d)"% (N, m))
    A = PETSc.Mat().create()
    A.setSizes([N, N])
    A.setFromOptions()
    Istart, Iend = A.getOwnershipRange()
    ix = 0
    cst = 0.5/(m-1)
    for i in range(1,m+1):
       jmax = m-i+1
       for j in range(1,jmax+1):
           ix = ix + 1
           if ix-1<Istart or ix>Iend:
               continue  # compute only owned rows
           if j!=jmax:
               pd = cst*(i+j-1)
               # north
               if i==1:
                   A[ix-1,ix] = 2*pd
               else:
                   A[ix-1,ix] = pd
               # east
               if j==1:
                   A[ix-1,ix+jmax-1] = 2*pd
               else:
                   A[ix-1,ix+jmax-1] = pd
           # south
           pu = 0.5 - cst*(i+j-3)
           if j>1:
               A[ix-1,ix-2] = pu
           # west
           if i>1:
               A[ix-1,ix-jmax-2] = pu

    A.assemble()
    return A

# The following function solves the problem. This case is quite different
# from eigenproblems, and is more similar to solving a linear system of
# equations with `KSP <petsc4py.PETSc.KSP>`. To configure the problem we
# must provide the matrix and the function (the exponential in this case).
# Note how the internal `FN` object is extracted from the `MFN` solver.
# Also, it is often necessary to specify a scale factor, which in this
# case represents the time for which we want to obtain the evolved state.
# Once the solver is set up, we call `solve() <MFN.solve()>` passing the
# right-hand side vector ``b`` and the solution vector ``x``.

def solve_exp(t, A, b, x):
    # Setup the solver
    M = SLEPc.MFN().create()
    M.setOperator(A)
    f = M.getFN()
    f.setType(SLEPc.FN.Type.EXP)
    f.setScale(t)
    M.setTolerances(1e-7)
    M.setFromOptions()
    # Solve the problem
    M.solve(b,x)

    its = M.getIterationNumber()
    Print("Number of iterations of the method: %i" % its)
    sol_type = M.getType()
    Print("Solution method: %s" % sol_type)
    ncv = M.getDimensions()
    Print("")
    Print("Subspace dimension: %i" % ncv)
    tol, maxit = M.getTolerances()
    Print("Stopping condition: tol=%.4g, maxit=%d" % (tol, maxit))
    Print("Computed vector at time t=%.4g has norm %g" % (t.real, x.norm()))
    Print("")

# The main program processes the command-line option ``m`` (size of the
# grid), builds the matrix and calls the solver. Note how the vectors are
# created from the matrix. In this case, the right-hand side vector is
# the first element of the canonical basis.

if __name__ == '__main__':
    opts = PETSc.Options()
    m = opts.getInt('m', 15)
    A = build_matrix(m)   # transition probability matrix
    x, b = A.createVecs()
    x.set(0)
    b.set(0)
    b[0] = 1
    b.assemble()
    t = 2
    solve_exp(t, A, b, x)    # compute x=exp(t*A)*b
    A = None
    b = x = None

