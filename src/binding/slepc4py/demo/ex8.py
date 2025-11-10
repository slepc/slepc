# ex8.py: Nonlinear eigenproblem with split form
# ==============================================
#
# This example solves a nonlinear eigenvalue problem where the
# nonlinear function is expressed in split form.
#
# We want to solve the following parabolic partial differential
# equation with time delay :math:`\tau`
#
# .. math::
#
#    u_t    &= u_{xx} + a u(t) + b u(t-\tau) \\
#    u(0,t) &= u(\pi,t) = 0
#
# with :math:`a = 20` and :math:`b(x) = -4.1+x (1-e^{x-\pi})`.
#
# Discretization leads to a DDE of dimension :math:`n`
#
# .. math::
#
#    -u' = A u(t) + B u(t-\tau)
#
# which results in the nonlinear eigenproblem
#
# .. math::
#
#    (-\lambda I + A + e^{-\tau\lambda}B)u = 0.
#
# The full source code for this demo can be `downloaded here
# <../_static/ex8.py>`__.

# Initialization is similar to previous examples. In this case we also
# need to import some math symbols.

import sys, slepc4py
slepc4py.init(sys.argv)

from petsc4py import PETSc
from slepc4py import SLEPc
from numpy import exp
from math import pi

Print = PETSc.Sys.Print

# This script has two command-line options: the discretization size ``n``
# and the time delay ``tau``.

opts = PETSc.Options()
n = opts.getInt('n', 128)
tau = opts.getReal('tau', 0.001)
a = 20
h = pi/(n+1)

# Next we have to set up the solver. In this case, we are going to
# represent the nonlinear problem in split form, i.e., as a sum of
# terms made of a constant matrix multiplied by a scalar nonlinear
# function.

nep = SLEPc.NEP().create()

# The first term involves the identity matrix.

Id = PETSc.Mat().createConstantDiagonal([n, n], 1.0)

# The second term has a tridiagonal matrix obtained from the
# discretization, :math:`A = \frac{1}{h^2}\operatorname{tridiag}(1,-2,1) + a I`.

A = PETSc.Mat().create()
A.setSizes([n, n])
A.setFromOptions()
rstart, rend = A.getOwnershipRange()
vd = -2.0/(h*h)+a
vo = 1.0/(h*h)
if rstart == 0:
  A[0, :2] = [vd, vo]
  rstart += 1
if rend == n:
  A[n-1, -2:] = [vo, vd]
  rend -= 1
for i in range(rstart, rend):
  A[i, i-1:i+2] = [vo, vd, vo]
A.assemble()

# The third term includes a diagonal matrix :math:`B = \operatorname{diag}(b(x_i))`.

B = PETSc.Mat().create()
B.setSizes([n, n])
B.setFromOptions()
rstart, rend = B.getOwnershipRange()
for i in range(rstart, rend):
  xi = (i+1)*h
  B[i, i] = -4.1+xi*(1.0-exp(xi-pi));
B.assemble()
B.setOption(PETSc.Mat.Option.HERMITIAN, True)

# Apart from the matrices, we have to create the functions, represented with
# `FN` objects: :math:`f_1=-\lambda, f_2=1, f_3=\exp(-\tau\lambda)`.

f1 = SLEPc.FN().create()
f1.setType(SLEPc.FN.Type.RATIONAL)
f1.setRationalNumerator([-1, 0])
f2 = SLEPc.FN().create()
f2.setType(SLEPc.FN.Type.RATIONAL)
f2.setRationalNumerator([1])
f3 = SLEPc.FN().create()
f3.setType(SLEPc.FN.Type.EXP)
f3.setScale(-tau)

# Put all the information together to define the split operator. Note that
# ``A`` is passed first so that `SUBSET <petsc4py.PETSc.Mat.Structure.SUBSET>`
# nonzero_pattern can be used.

nep.setSplitOperator([A, Id, B], [f2, f1, f3], PETSc.Mat.Structure.SUBSET)

# Now we can set some options and call the solver.

nep.setTolerances(tol=1e-9)
nep.setDimensions(1)
nep.setFromOptions()

nep.solve()

# Once the solver has finished, we print some information together with
# the computed solution. For each computed eigenpair, we print the
# eigenvalue and the residual norm.

its = nep.getIterationNumber()
Print("Number of iterations of the method: %i" % its)
sol_type = nep.getType()
Print("Solution method: %s" % sol_type)
nev, ncv, mpd = nep.getDimensions()
Print("")
Print("Subspace dimension: %i" % ncv)
tol, maxit = nep.getTolerances()
Print("Stopping condition: tol=%.4g" % tol)
Print("")

nconv = nep.getConverged()
Print( "Number of converged eigenpairs %d" % nconv )

if nconv > 0:
  x = Id.createVecs('right')
  x.set(1.0)
  Print()
  Print("        k              ||T(k)x||")
  Print("----------------- ------------------")
  for i in range(nconv):
    k = nep.getEigenpair(i, x)
    res = nep.computeError(i)
    if k.imag != 0.0:
      Print( " %9f%+9f j %12g" % (k.real, k.imag, res) )
    else:
      Print( " %12f       %12g" % (k.real, res) )
  Print()
