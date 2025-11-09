# ex13.py: Nonlinear eigenproblem with contour integral
# =====================================================
#
# This example solves a nonlinear eigenvalue problem with the `NEP
# <slepc4py.SLEPc.NEP>` module configured to apply a contour integral
# method.
#
# The problem arises from the PDE
#
# .. math::
#
#    u_{xx}(x) + n_c^2 \lambda^2 u(x) + g(\lambda) D_0 \lambda^2 u(x) = 0
#
# where
#
# .. math::
#
#    g(\lambda) &= g_t/(\lambda-k_a + i g_t), \\
#    k_a &= 8.0, \\
#    g_t &= 0.5, \\
#    D_0 &= 0.5, \\
#    n_c &= 1.2,
#
# and the boundary conditions are
#
# .. math::
#
#    u(0) &= 0 \\
#    u_x(1) &= i \lambda u(1).
#
# For the discretization, :math:`n` grid points are used, :math:`x_1=0,\dots,x_n=1`,
# with step size :math:`h = 1/(n-1)`. Hence,
#
# .. math::
#
#    u_{xx}(x_i) &= \frac{1}{h^2} (u_{i-1} - 2 u_i + u_{i+1}) \\
#                &= \frac{1}{h^2} [1, -2, 1] [u_{i-1}, u_i, u_{i+1}]^T.
#
# The boundary condition at :math:`x=0` is :math:`u_1=0`, and at :math:`x=1`:
#
# .. math::
#
#    u'(1) &= 1/2 ( u'(1+h/2) + u'(1-h/2) ) \\
#          &= 1/2 ( (u_{n+1}-u_n)/h  + (u_n-u_{n-1})/h ) \\
#          &= 1/(2h) (u_{n+1} - u_{n-1}) = i\lambda u_n,
#
# and therefore
#
# .. math::
#
#    u_{n+1} = 2i h \lambda u_n + u_{n-1}.
#
# The Laplace term for :math:`u_n` is
#
# .. math::
#
#    &= \frac{1}{h^2} (u_{n-1} - 2u_n + u_{n+1}) \\
#    &= \frac{1}{h^2} (u_{n-1} - 2u_n + 2ih\lambda u_n + u_{n-1}) \\
#    &= \frac{1}{h^2} (2 u_{n-1} + (2ih\lambda - 2) u_n) \\
#    &= \frac{1}{h^2} [2, 2ih\lambda -2] [u_{n-1}, u_n]^T.
#
# The above discretization allows us to write the nonlinear PDE in the
# following split-operator form
#
# .. math::
#
#    \{A + \lambda^2 n_c^2 I_d + g(\lambda) \lambda^2 D_0 I_d + 2i\lambda/h D\} u = 0
#
# so :math:`f_1 = 1`, :math:`f_2 = n_c^2 \lambda^2`, :math:`f_3 = g(\lambda)
# \lambda^2D_0`, :math:`f_4 = 2i\lambda/h`, with coefficient matrices
#
# .. math::
#
#    A = \begin{bmatrix}
#        1 & 0 & 0 & \dots & 0 \\
#        0 & * & * & \dots & 0 \\
#        0 & * & * & \dots & 0 \\
#        \vdots & \vdots & \vdots & \ddots & \vdots \\
#        0 & 0 & 0 & \dots & *
#        \end{bmatrix},
#    I_d = \begin{bmatrix}
#        0 & 0 & 0 & \dots & 0 \\
#        0 & 1 & 0 & \dots & 0 \\
#        0 & 0 & 1 & \dots & 0 \\
#        \vdots & \vdots & \vdots & \ddots & \vdots \\
#        0 & 0 & 0 & \dots & 0
#        \end{bmatrix},
#    D = \begin{bmatrix}
#        0 & 0 & 0 & \dots & 0 \\
#        0 & 0 & 0 & \dots & 0 \\
#        0 & 0 & 0 & \dots & 0 \\
#        \vdots & \vdots & \vdots & \ddots & \vdots \\
#        0 & 0 & 0 & \dots & 1
#        \end{bmatrix}.
#
# Contributed by: Thomas Hisch.
#
# The full source code for this demo can be `downloaded here
# <../_static/ex13.py>`__.

# Initialization by importing slepc4py, petsc4py, numpy and scipy.

import sys

import slepc4py

slepc4py.init(sys.argv)  # isort:skip

import numpy as np

try:
    import scipy
    import scipy.optimize
except ImportError:
    scipy = None

from petsc4py import PETSc
from slepc4py import SLEPc

Print = PETSc.Sys.Print

# Check that the selected PETSc/SLEPc are built with complex scalars.

if not np.issubdtype(PETSc.ScalarType, np.complexfloating):
    Print("Demo should only be executed with complex PETSc scalars")
    exit(0)

# Long function that defines the nonlinear eigenproblem and solves it.

def solve(n):
    L = 1.0
    h = L / (n - 1)
    nc = 1.2
    ka = 10.0
    gt = 4.0
    D0 = 0.5

    A = PETSc.Mat().create()
    A.setSizes([n, n])
    A.setFromOptions()
    A.setOption(PETSc.Mat.Option.HERMITIAN, False)

    rstart, rend = A.getOwnershipRange()
    d0, d1, d2 = (
        1 / h**2,
        -2 / h**2,
        1 / h**2,
    )
    Print(f"dterms={(d0, d1, d2)}")

    if rstart == 0:
        # dirichlet boundary condition at the left lead
        A[0, 0] = 1.0
        A[0, 1] = 0.0
        A[1, 0] = 0.0

        A[1, 1] = d1
        A[1, 2] = d2
        rstart += 2
    if rend == n:
        # at x=1.0 neumann boundary condition (not handled here but in a
        # different matrix (D))
        A[n - 1, n - 2] = 2.0 / h**2
        A[n - 1, n - 1] = (-2) / h**2  # + 2j*k*h / h**2 (neumann)
        rend -= 1

    for i in range(rstart, rend):
        A[i, i - 1 : i + 2] = [d0, d1, d2]

    A.assemble()

    Id = PETSc.Mat().create()
    Id.setSizes([n, n])
    Id.setFromOptions()
    Id.setOption(PETSc.Mat.Option.HERMITIAN, True)
    rstart, rend = Id.getOwnershipRange()
    if rstart == 0:
        # due to dirichlet BC
        rstart += 1
    for i in range(rstart, rend):
        Id[i, i] = 1.0
    Id.assemble()

    D = PETSc.Mat().create()
    D.setSizes([n, n])
    D.setFromOptions()
    D.setOption(PETSc.Mat.Option.HERMITIAN, True)
    _, rend = D.getOwnershipRange()
    if rend == n:
        D[n - 1, n - 1] = 1
    D.assemble()

    Print(f"DOF: {A.getInfo()['nz_used']}, MEM: {A.getInfo()['memory']}")

    f1 = SLEPc.FN().create()
    f1.setType(SLEPc.FN.Type.RATIONAL)
    f1.setRationalNumerator([1.0])
    f2 = SLEPc.FN().create()
    f2.setType(SLEPc.FN.Type.RATIONAL)
    f2.setRationalNumerator([nc**2, 0.0, 0.0])
    f3 = SLEPc.FN().create()
    f3.setType(SLEPc.FN.Type.RATIONAL)
    f3.setRationalNumerator([D0 * gt, 0.0, 0.0])
    f3.setRationalDenominator([1.0, -ka + 1j * gt])
    f4 = SLEPc.FN().create()
    f4.setType(SLEPc.FN.Type.RATIONAL)
    f4.setRationalNumerator([2j / h, 0])

    # Setup the solver
    nep = SLEPc.NEP().create()

    nep.setSplitOperator(
        [A, Id, Id, D],
        [f1, f2, f3, f4],
        PETSc.Mat.Structure.SUBSET,
    )

    # Customize options
    nep.setTolerances(tol=1e-7)

    nep.setDimensions(nev=24)
    nep.setType(SLEPc.NEP.Type.CISS)

    # the rg params are chosen s.t. the singularity at k = ka - 1j*gt is
    # outside of the contour.
    radius = 3 * gt
    vscale = 0.5 * gt / radius
    rg_params = (ka, 3 * gt, vscale)
    R = nep.getRG()
    R.setType(SLEPc.RG.Type.ELLIPSE)
    Print(f"RG params: {rg_params}")
    R.setEllipseParameters(*rg_params)

    nep.setFromOptions()

    # Solve the problem
    nep.solve()

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
    Print("Number of converged eigenpairs %d" % nconv)

    x = A.createVecs("right")

    evals = []
    modes = []
    if nconv > 0:
        Print()
        Print("        lam            ||T(lam)x||   |lam-lam_exact|/|lam_exact| ")
        Print("--------------------- ------------- -----------------------------")
        for i in range(nconv):
            lam = nep.getEigenpair(i, x)
            error = nep.computeError(i)

            def eigenvalue_error_term(k):
                gkmu = gt / (k - ka + 1j * gt)
                nceff = np.sqrt(nc**2 + gkmu * D0)
                return -1j / np.tan(nceff * k * L) - 1 / nceff

            # compute the expected_eigenvalue

            # we assume that the numerically calculated eigenvalue is close to
            # the exact one, which we can determine using a Newton-Raphson
            # method.
            if scipy:
                expected_lam = scipy.optimize.newton(
                    eigenvalue_error_term, np.complex128(lam), rtol=1e-11
                )
                rel_err = abs(lam - expected_lam) / abs(expected_lam)
                rel_err = "%6g" % rel_err
            else:
                rel_err = "scipy not installed"

            Print(" %9f%+9f j %12g   %s" % (lam.real, lam.imag, error, rel_err))

            evals.append(lam)
            modes.append(x.getArray().copy())
        Print()

    return np.asarray(evals), rg_params, ka, gt

# The main function reads the problem size ``n`` from the command line,
# solves the problem with the above function, and then plots the computed
# eigenvalues.

def main():
    opts = PETSc.Options()
    n = opts.getInt("n", 256)
    Print(f"n={n}")

    evals, rg_params, ka, gt = solve(n)

    if not opts.getBool("ploteigs", True) or PETSc.COMM_WORLD.getRank():
        return

    try:
        import matplotlib.pyplot as plt
        from matplotlib.patches import Ellipse
    except ImportError:
        print("plot is not shown, because matplotlib is not installed")
    else:
        fig, ax = plt.subplots()
        ax.plot(evals.real, evals.imag, "x")

        height = 2 * rg_params[1] * rg_params[2]
        ellipse = Ellipse(
            xy=(rg_params[0], 0.0),
            width=rg_params[1] * 2,
            height=height,
            edgecolor="r",
            fc="None",
            lw=2,
        )
        ax.add_patch(ellipse)

        ax.grid()
        ax.legend()
        plt.show()


if __name__ == "__main__":
    main()
