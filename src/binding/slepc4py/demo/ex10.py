# ex10.py: Laplace problem using the Proper Orthogonal Decomposition
# ==================================================================
#
# This example program solves the Laplace problem using the Proper Orthogonal
# Decomposition (POD) reduced-order modeling technique. For a full description
# of the technique the reader is referred to [1]_, [2]_.
#
# The method is split into an offline (computationally intensive) and an online
# (computationally cheap) phase. This has many applications including real-time
# simulation, uncertainty quantification and inverse problems, where similar
# models must be evaluated quickly and many times.
#
# Offline phase:
#
# 1. A set of solution snapshots of the 1D Laplace problem in the full
#    problem space are constructed and assembled into the columns of a dense
#    matrix :math:`S`.
# 2. A standard eigenvalue decomposition is performed on the
#    matrix :math:`S^T S`.
# 3. The eigenvectors and eigenvalues are projected back to the
#    original eigenvalue problem :math:`S`.
# 4. The leading eigenvectors then form the POD basis.
#
# Online phase:
#
# 1. The operator corresponding to the discrete Laplacian is
#    projected onto the POD basis.
# 2. The operator corresponding to the right-hand side is
#    projected onto the POD basis.
# 3. The reduced (dense) problem expressed in the POD basis is solved.
# 4. The reduced solution is projected back to the full
#    problem space.
#
# Contributed by: Elisa Schenone, Jack S. Hale.
#
# The full source code for this demo can be `downloaded here
# <../_static/ex10.py>`__.

# Initialization is similar to previous examples, but importing some
# additional modules.

try: range = xrange
except: pass

import sys, slepc4py
slepc4py.init(sys.argv)

from petsc4py import PETSc
from slepc4py import SLEPc
import numpy

import random
import math

# This function builds the discrete Laplacian operator in 1 dimension
# with homogeneous Dirichlet boundary conditions.

def construct_operator(m):
    # Create matrix for 1D Laplacian operator
    A = PETSc.Mat().create(PETSc.COMM_SELF)
    A.setSizes([m, m])
    A.setFromOptions()
    # Fill matrix
    hx = 1.0/(m-1) # x grid spacing
    diagv = 2.0/hx
    offdx = -1.0/hx
    Istart, Iend = A.getOwnershipRange()
    for i in range(Istart, Iend):
        if i != 0 and i != (m - 1):
            A[i, i] = diagv
            if i > 1: A[i, i - 1] = offdx
            if i < m - 2: A[i, i + 1] = offdx
        else:
            A[i, i] = 1.

    A.assemble()

    return A

# Set bell-shape function as the solution of the Laplacian problem in
# 1 dimension with homogeneous Dirichlet boundary conditions and
# compute the associated discrete RHS.

def set_problem_rhs(m):
    # Create 1D mass matrix operator
    M = PETSc.Mat().create(PETSc.COMM_SELF)
    M.setSizes([m, m])
    M.setFromOptions()
    # Fill matrix
    hx = 1.0/(m-1) # x grid spacing
    diagv = hx/3
    offdx = hx/6
    Istart, Iend = M.getOwnershipRange()
    for i in range(Istart, Iend):
        if i != 0 and i != (m - 1):
            M[i, i] = 2*diagv
        else:
            M[i, i] = diagv
        if i > 1: M[i, i - 1] = offdx
        if i < m - 2: M[i, i + 1] = offdx

    M.assemble()

    x_0 = 0.3
    x_f = 0.7
    mu = x_0 + (x_f - x_0)*random.random()
    sigma = 0.1**2
    uex, f = M.createVecs()
    for j in range(Istart, Iend):
        value = 2/sigma * math.exp(-(hx*j - mu)**2/sigma) * (1 - 2/sigma * (hx*j - mu)**2 )
        f.setValue(j, value)
        value = math.exp(-(hx*j - mu)**2/sigma)
        uex.setValue(j, value)
    f.assemble()
    uex.assemble()

    RHS = f.duplicate()
    M.mult(f, RHS)
    RHS.setValue(0, 0.)
    RHS.setValue(m-1, 0.)
    RHS.assemble()

    return RHS, uex

# Solve 1D Laplace problem with FEM.

def solve_laplace_problem(A, RHS):
    u = A.createVecs('right')
    r, c = A.getOrdering("natural")
    A.factorILU(r, c)
    A.solve(RHS, u)
    A.setUnfactored()
    return u

# Solve 1D Laplace problem with POD (dense matrix).

def solve_laplace_problem_pod(A, RHS, u):
    ksp = PETSc.KSP().create(PETSc.COMM_SELF)
    ksp.setOperators(A)
    ksp.setType('preonly')
    pc = ksp.getPC()
    pc.setType('none')
    ksp.setFromOptions()

    ksp.solve(RHS, u)

    return u

# Set ``N`` solution of the 1D Laplace problem as columns of a matrix
# (snapshot matrix).
#
# Note: For simplicity we do not perform a linear solve, but use
# some analytical solution:
# :math:`z(x) = \exp(-(x - \mu)^2 / \sigma)`.

def construct_snapshot_matrix(A, N, m):
    snapshots = PETSc.Mat().create(PETSc.COMM_SELF)
    snapshots.setSizes([m, N])
    snapshots.setType('seqdense')

    Istart, Iend = snapshots.getOwnershipRange()
    hx = 1.0/(m - 1)
    x_0 = 0.3
    x_f = 0.7
    sigma = 0.1**2
    for i in range(N):
        mu = x_0 + (x_f - x_0)*random.random()
        for j in range(Istart, Iend):
            value = math.exp(-(hx*j - mu)**2/sigma)
            snapshots.setValue(j, i, value)
    snapshots.assemble()

    return snapshots

# Solve the eigenvalue problem: the eigenvectors of this problem form the
# POD basis.

def solve_eigenproblem(snapshots, N):
    print('Solving POD basis eigenproblem using eigensolver...')

    Es = SLEPc.EPS()
    Es.create(PETSc.COMM_SELF)
    Es.setDimensions(N)
    Es.setProblemType(SLEPc.EPS.ProblemType.NHEP)
    Es.setTolerances(1.0e-8, 500);
    Es.setKrylovSchurRestart(0.6)
    Es.setWhichEigenpairs(SLEPc.EPS.Which.LARGEST_REAL)
    Es.setOperators(snapshots)
    Es.setFromOptions()

    Es.solve()
    print('Solved POD basis eigenproblem.')
    return Es

# Function to project :math:`S^TS` eigenvectors to :math:`S` eigenvectors.

def project_STS_eigenvectors_to_S_eigenvectors(bvEs, S):
    sizes = S.getSizes()[0]
    N = bvEs.getActiveColumns()[1]
    bv = SLEPc.BV().create(PETSc.COMM_SELF)
    bv.setSizes(sizes, N)
    bv.setActiveColumns(0, N)
    bv.setFromOptions()

    tmpvec2 = S.createVecs('left')
    for i in range(N):
        tmpvec = bvEs.getColumn(i)
        S.mult(tmpvec, tmpvec2)
        bv.insertVec(i, tmpvec2)
        bvEs.restoreColumn(i, tmpvec)

    return bv

# Function to project the reduced space to the full space.

def project_reduced_to_full_space(alpha, bv):
    uu = bv.getColumn(0)
    uPOD = uu.duplicate()
    bv.restoreColumn(0,uu)

    scatter, Wr = PETSc.Scatter.toAll(alpha)
    scatter.begin(alpha, Wr, PETSc.InsertMode.INSERT, PETSc.ScatterMode.FORWARD)
    scatter.end(alpha, Wr, PETSc.InsertMode.INSERT, PETSc.ScatterMode.FORWARD)
    PODcoeff = Wr.getArray(readonly=1)

    bv.multVec(1., 0., uPOD, PODcoeff)

    return uPOD

# The main function realizes the full procedure.

def main():
    problem_dim = 200
    num_snapshots = 30
    num_pod_basis_functions = 8

    assert(num_pod_basis_functions <= num_snapshots)

    A = construct_operator(problem_dim)
    S = construct_snapshot_matrix(A, num_snapshots, problem_dim)

    # Instead of solving the SVD of S, we solve the standard
    # eigenvalue problem on S.T*S
    STS = S.transposeMatMult(S)

    Es = solve_eigenproblem(STS, num_pod_basis_functions)
    nconv = Es.getConverged()
    print('Number of converged eigenvalues: %i' % nconv)
    Es.view()

    # get the EPS solution in a BV object
    bvEs = Es.getBV()
    bvEs.setActiveColumns(0, num_pod_basis_functions)

    # set the bv POD basis
    bv = project_STS_eigenvectors_to_S_eigenvectors(bvEs, S)
    # rescale the eigenvectors
    for i in range(num_pod_basis_functions):
        ll = Es.getEigenvalue(i)
        print('Eigenvalue '+str(i)+': '+str(ll.real))
        bv.scaleColumn(i,1.0/math.sqrt(ll.real))

    print('--------------------------------')
    # Verify that the active columns of bv form an orthonormal subspace, i.e. that X^H*X = Id
    print('Check that bv.dot(bv) is close to the identity matrix')
    XtX = bv.dot(bv)
    XtX.view()
    XtX_array = XtX.getDenseArray()
    n,m = XtX_array.shape
    assert numpy.allclose(XtX_array, numpy.eye(n, m))
    print('--------------------------------')
    print('Solve the problem with POD')

    # Project the linear operator A
    Ared = bv.matProject(A,bv)

    # Set the RHS and the exact solution
    RHS, uex = set_problem_rhs(problem_dim)

    # Project the RHS on the POD basis
    RHSred = PETSc.Vec().createWithArray(bv.dotVec(RHS))

    # Solve the problem with POD
    alpha = Ared.createVecs('right')
    alpha = solve_laplace_problem_pod(Ared,RHSred,alpha)

    # Project the POD solution back to the FE space
    uPOD = project_reduced_to_full_space(alpha, bv)

    # Compute the L2 and Linf norm of the error
    error = uex.copy()
    error.axpy(-1,uPOD)
    errorL2 = math.sqrt(error.dot(error).real)
    print('The L2-norm of the error is: '+str(errorL2))

    print("NORMAL END")

if __name__ == '__main__':
    main()

# .. [1] K. Kunisch and S. Volkwein. Galerkin proper orthogonal decomposition methods
#    for a general equation in fluid dynamics. SIAM Journal on Numerical Analysis,
#    40(2):492-515, 2003.
# .. [2] S. Volkwein, Optimal control of a phase-field model using the proper orthogonal
#    decomposition, Z. Angew. Math. Mech., 81:83-97, 2001.
