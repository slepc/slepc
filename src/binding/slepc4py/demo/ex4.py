# ex4.py: Singular value decomposition of the Lauchli matrix
# ==========================================================
#
# This example illustrates the use of the SVD solver in slepc4py. It
# computes singular values and vectors of the Lauchli matrix, whose
# condition number depends on a parameter ``mu``.
#
# The full source code for this demo can be `downloaded here
# <../_static/ex4.py>`__.

# Initialization is similar to previous examples.

try: range = xrange
except: pass

import sys, slepc4py
slepc4py.init(sys.argv)

from petsc4py import PETSc
from slepc4py import SLEPc

# This example takes two command-line arguments, the matrix size ``n``
# and the ``mu`` parameter.

opts = PETSc.Options()
n  = opts.getInt('n', 30)
mu = opts.getReal('mu', 1e-6)

PETSc.Sys.Print( "Lauchli singular value decomposition, (%d x %d) mu=%g\n" % (n+1,n,mu) )

# Create the matrix and fill its nonzero entries. Every MPI process will
# insert its locally owned part only.

A = PETSc.Mat(); A.create()
A.setSizes([n+1, n])
A.setFromOptions()

rstart, rend = A.getOwnershipRange()

for i in range(rstart, rend):
  if i==0:
    for j in range(n):
      A[0,j] = 1.0
  else:
    A[i,i-1] = mu

A.assemble()

# The singular value solver is similar to the eigensolver used in previous
# examples. In this case, we select the thick-restart Lanczos
# bidiagonalization method.

S = SLEPc.SVD(); S.create()

S.setOperator(A)
S.setType(S.Type.TRLANCZOS)
S.setFromOptions()

S.solve()

# After solve, we print some informative data and extract the computed
# solution, showing the list of singular values and the corresponding
# residual errors.

Print = PETSc.Sys.Print

Print( "******************************" )
Print( "*** SLEPc Solution Results ***" )
Print( "******************************\n" )

svd_type = S.getType()
Print( "Solution method: %s" % svd_type )

its = S.getIterationNumber()
Print( "Number of iterations of the method: %d" % its )

nsv, ncv, mpd = S.getDimensions()
Print( "Number of requested singular values: %d" % nsv )

tol, maxit = S.getTolerances()
Print( "Stopping condition: tol=%.4g, maxit=%d" % (tol, maxit) )

nconv = S.getConverged()
Print( "Number of converged approximate singular triplets %d" % nconv )

if nconv > 0:
  v, u = A.createVecs()
  Print()
  Print("    sigma       residual norm ")
  Print("-------------  ---------------")
  for i in range(nconv):
    sigma = S.getSingularTriplet(i, u, v)
    error = S.computeError(i)
    Print( "   %6f     %12g" % (sigma, error) )
  Print()
