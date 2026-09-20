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

import sys
import slepc4py

slepc4py.init(sys.argv)

from petsc4py import PETSc
from slepc4py import SLEPc

Print = PETSc.Sys.Print

# This example takes two command-line arguments, the matrix size ``n``
# and the ``mu`` parameter.

opts = PETSc.Options()
n = opts.getInt('n', 30)
mu = opts.getReal('mu', 1e-6)

Print(f'Lauchli singular value decomposition, ({n + 1} x {n}) mu={mu}\n')

# Create the matrix and fill its nonzero entries. Every MPI process will
# insert its locally owned part only.

A = PETSc.Mat().create()
A.setSizes([n + 1, n])
A.setFromOptions()

rstart, rend = A.getOwnershipRange()

for i in range(rstart, rend):
    if i == 0:
        for j in range(n):
            A[0, j] = 1.0
    else:
        A[i, i - 1] = mu

A.assemble()

# The singular value solver is similar to the eigensolver used in previous
# examples. In this case, we select the thick-restart Lanczos
# bidiagonalization method.

S = SLEPc.SVD().create()

S.setOperator(A)
S.setType(S.Type.TRLANCZOS)
S.setFromOptions()

S.solve()

# After solve, we print some informative data and extract the computed
# solution, showing the list of singular values and the corresponding
# residual errors.

Print('******************************')
Print('*** SLEPc Solution Results ***')
Print('******************************\n')

svd_type = S.getType()
Print(f'Solution method: {svd_type}')

its = S.getIterationNumber()
Print(f'Number of iterations of the method: {its}')

nsv, _ncv, _mpd = S.getDimensions()
Print(f'Number of requested singular values: {nsv}')

tol, maxit = S.getTolerances()
Print(f'Stopping condition: tol={tol:.4g}, maxit={maxit}')

nconv = S.getConverged()
Print(f'Number of converged approximate singular triplets {nconv}')

if nconv > 0:
    v, u = A.createVecs()
    Print()
    Print('    sigma       residual norm ')
    Print('-------------  ---------------')
    for i in range(nconv):
        sigma = S.getSingularTriplet(i, u, v)
        error = S.computeError(i)
        Print(f'   {sigma:6f}     {error:12g}')
    Print()
