# ex1.py: Standard symmetric eigenproblem for the 1-D Laplacian
# =============================================================
#
# This example computes eigenvalues and eigenvectors of the discrete Laplacian
# on a one-dimensional domain with finite differences.
#
# The full source code for this demo can be `downloaded here
# <../_static/ex1.py>`__.

# The first thing to do is initialize the libraries. This is normally not
# required, as it is done automatically at import time. However, if you want to
# gain access to the facilities for setting command-line options, the
# following lines must be executed by the main script prior to any petsc4py or
# slepc4py calls:

import sys, slepc4py
slepc4py.init(sys.argv)

# Next, we have to import the relevant modules. Normally, both PETSc and SLEPc
# modules have to be imported in all slepc4py programs. It may be useful to
# import NumPy as well:

from petsc4py import PETSc
from slepc4py import SLEPc
import numpy

# At this point, we can use any petsc4py and slepc4py operations. For instance,
# the following lines allow the user to specify an integer command-line
# argument ``n`` with a default value of 30 (see below for example usage
# of command-line options):

opts = PETSc.Options()
n = opts.getInt('n', 30)

# It is necessary to build a matrix to define an eigenproblem (or two in the
# case of generalized eigenproblems). The following fragment of code creates
# the matrix object and then fills the non-zero elements one by one. The matrix
# of this particular example is tridiagonal, with value 2 in the diagonal, and
# -1 in off-diagonal positions. See petsc4py documentation for details about
# matrix objects:

A = PETSc.Mat(); A.create()
A.setSizes([n, n])
A.setFromOptions()

rstart, rend = A.getOwnershipRange()

# first row
if rstart == 0:
  A[0, :2] = [2, -1]
  rstart += 1
# last row
if rend == n:
  A[n-1, -2:] = [-1, 2]
  rend -= 1
# other rows
for i in range(rstart, rend):
  A[i, i-1:i+2] = [-1, 2, -1]

A.assemble()

# The solver object is created in a similar way as other objects in petsc4py:

E = SLEPc.EPS(); E.create()

# Once the object is created, the eigenvalue problem must be specified. At
# least one matrix must be provided. The problem type must be indicated as
# well, in this case it is HEP (Hermitian eigenvalue problem). Apart from
# these, other settings could be provided here (for instance, the tolerance for
# the computation). After all options have been set, the user should call the
# `setFromOptions() <EPS.setFromOptions()>` operation, so that any options
# specified at run time in the command line are passed to the solver object:

E.setOperators(A)
E.setProblemType(SLEPc.EPS.ProblemType.HEP)

history = []
def monitor(eps, its, nconv, eig, err):
    if nconv<len(err): history.append(err[nconv])
E.setMonitor(monitor)

E.setFromOptions()

# After that, the `solve() <EPS.solve()>` method will run the selected
# eigensolver, keeping the solution stored internally:

E.solve()

# Once the computation has finished, we are ready to print the results. First,
# some informative data can be retrieved from the solver object:

Print = PETSc.Sys.Print

Print()
Print("******************************")
Print("*** SLEPc Solution Results ***")
Print("******************************")
Print()

its = E.getIterationNumber()
Print( "Number of iterations of the method: %d" % its )

eps_type = E.getType()
Print( "Solution method: %s" % eps_type )

nev, ncv, mpd = E.getDimensions()
Print( "Number of requested eigenvalues: %d" % nev )

tol, maxit = E.getTolerances()
Print( "Stopping condition: tol=%.4g, maxit=%d" % (tol, maxit) )

# For retrieving the solution, it is necessary to find out how many eigenpairs
# have converged to the requested precision:

nconv = E.getConverged()
Print( "Number of converged eigenpairs %d" % nconv )

# For each of the ``nconv`` eigenpairs, we can retrieve the eigenvalue ``k``,
# and the eigenvector, which is represented by means of two petsc4py vectors
# ``vr`` and ``vi`` (the real and imaginary part of the eigenvector, since for
# real matrices the eigenvalue and eigenvector may be complex). In this
# example we know that both the eigenvalue and eigenvector are real, so only
# one vector ``v`` is needed. We also compute the corresponding relative errors
# in order to make sure that the computed solution is indeed correct:

if nconv > 0:
  # Create the results vectors
  v, _ = A.createVecs()
  #
  Print()
  Print("        k          ||Ax-kx||/||kx|| ")
  Print("----------------- ------------------")
  for i in range(nconv):
    k = E.getEigenpair(i, v)
    error = E.computeError(i)
    Print( " %12f       %12g" % (k, error) )
  Print()

#
# Example of command-line usage
# -----------------------------
#
# Now we illustrate how to specify command-line options in order to extract the
# full potential of slepc4py.
#
# A simple execution of the ``demo/ex1.py`` script will result in the following
# output:
#
# .. code-block:: console
#
#   $ python demo/ex1.py
#
#   ******************************
#   *** SLEPc Solution Results ***
#   ******************************
#
#   Number of iterations of the method: 4
#   Solution method: krylovschur
#   Number of requested eigenvalues: 1
#   Stopping condition: tol=1e-07, maxit=100
#   Number of converged eigenpairs 4
#
#       k          ||Ax-kx||/||kx||
#   ----------------- ------------------
#        3.989739        5.76012e-09
#        3.959060        1.41957e-08
#        3.908279        6.74118e-08
#        3.837916        8.34269e-08
#
# For specifying different setting for the solver parameters, we can use SLEPc
# command-line options with the ``-eps`` prefix. For instance, to change the number
# of requested eigenvalues and the tolerance:
#
# .. code-block:: console
#
#   $ python demo/ex1.py -eps_nev 10 -eps_tol 1e-11
#
# The method used by the solver object can also be set at run time:
#
# .. code-block:: console
#
#   $ python demo/ex1.py -eps_type subspace
#
# All the above settings can also be changed within the source code by making
# use of the appropriate slepc4py method. Since options can be set from within
# the code and the command-line, it is often useful to view the particular
# settings that are currently being used:
#
# .. code-block:: console
#
#   $ python demo/ex1.py -eps_view
#
#   EPS Object: 1 MPI process
#     type: krylovschur
#       50% of basis vectors kept after restart
#       using the locking variant
#     problem type: symmetric eigenvalue problem
#     selected portion of the spectrum: largest eigenvalues in magnitude
#     number of eigenvalues (nev): 1
#     number of column vectors (ncv): 16
#     maximum dimension of projected problem (mpd): 16
#     maximum number of iterations: 100
#     tolerance: 1e-08
#     convergence test: relative to the eigenvalue
#   BV Object: 1 MPI process
#     type: mat
#     17 columns of global length 30
#     orthogonalization method: classical Gram-Schmidt
#     orthogonalization refinement: if needed (eta: 0.7071)
#     block orthogonalization method: GS
#     doing matmult as a single matrix-matrix product
#   DS Object: 1 MPI process
#     type: hep
#     solving the problem with: Implicit QR method (_steqr)
#   ST Object: 1 MPI process
#     type: shift
#     shift: 0
#     number of matrices: 1
#
# Note that for computing eigenvalues of smallest magnitude we can use the
# option ``-eps_smallest_magnitude``, but for interior eigenvalues things are
# not so straightforward. One possibility is to try with harmonic extraction,
# for instance to get the eigenvalues closest to 0.6:
#
# .. code-block:: console
#
#   $ python demo/ex1.py -eps_harmonic -eps_target 0.6
#
# Depending on the problem, harmonic extraction may fail to converge. In those
# cases, it is necessary to specify a spectral transformation other than the
# default. In the command-line, this is indicated with the ``-st_`` prefix. For
# example, shift-and-invert with a value of the shift equal to 0.6 would be:
#
# .. code-block:: console
#
#   $ python demo/ex1.py -st_type sinvert -eps_target 0.6
