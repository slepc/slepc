"""Typing support."""

from __future__ import annotations  # novermin
from typing import (  # novermin
    Callable,
    Sequence,
    Literal,
)
from numpy.typing import (
    NDArray,
)
from petsc4py.PETSc import (
    Vec,
    Mat,
    KSP,
)
from .SLEPc import (
    BV,
    DS,
    FN,
    RG,
    ST,
    EPS,
    PEP,
    NEP,
    SVD,
    MFN,
    LME,
)

__all__ = [
    'Scalar',
    'ArrayInt',
    'ArrayReal',
    'ArrayComplex',
    'ArrayScalar',
    'LayoutSizeSpec',
    'EPSStoppingFunction',
    'EPSArbitraryFunction',
    'EPSEigenvalueComparison',
    'EPSMonitorFunction',
    'PEPStoppingFunction',
    'PEPMonitorFunction',
    'NEPStoppingFunction',
    'NEPMonitorFunction',
    'SVDStoppingFunction',
    'SVDMonitorFunction',
    'MFNMonitorFunction',
    'LMEMonitorFunction',
]

# --- PETSc Sys ---

Scalar = float | complex
"""Scalar type.

Scalars can be either `float` or `complex` (but not both) depending on how
PETSc was configured (``./configure --with-scalar-type=real|complex``).

"""

ArrayInt = NDArray[int]
"""Array of `int`."""

ArrayReal = NDArray[float]
"""Array of `float`."""

ArrayComplex = NDArray[complex]
"""Array of `complex`."""

ArrayScalar = NDArray[Scalar]
"""Array of `Scalar` numbers."""

LayoutSizeSpec = int | tuple[int, int]
"""`int` or 2-`tuple` of `int` describing the layout sizes.

   A single `int` indicates global size.
   A `tuple` of `int` indicates ``(local_size, global_size)``.
"""

# --- EPS ---

EPSStoppingFunction = Callable[[EPS, int, int, int, int], EPS.ConvergedReason]
"""`EPS` stopping test callback."""

EPSArbitraryFunction = Callable[[Scalar, Scalar, Vec, Vec, Scalar, Scalar], [Scalar, Scalar]]
"""`EPS` arbitrary selection callback."""

EPSEigenvalueComparison = Callable[[Scalar, Scalar, Scalar, Scalar], int]
"""`EPS` eigenvalue comparison callback."""

EPSMonitorFunction = Callable[[EPS, int, int, ArrayScalar, ArrayScalar, ArrayReal, int], None]
"""`EPS` monitor callback."""

# --- PEP ---

PEPStoppingFunction = Callable[[PEP, int, int, int, int], PEP.ConvergedReason]
""":py:class:`PEP <slepc4py.SLEPc.PEP>` stopping test callback."""

PEPMonitorFunction = Callable[[PEP, int, int, ArrayScalar, ArrayScalar, ArrayReal, int], None]
""":py:class:`PEP <slepc4py.SLEPc.PEP>` monitor callback."""

# --- NEP ---

NEPStoppingFunction = Callable[[NEP, int, int, int, int], NEP.ConvergedReason]
""":py:class:`NEP <slepc4py.SLEPc.NEP>` stopping test callback."""

NEPMonitorFunction = Callable[[NEP, int, int, ArrayScalar, ArrayScalar, ArrayReal, int], None]
""":py:class:`NEP <slepc4py.SLEPc.NEP>` monitor callback."""

NEPFunction = Callable[[NEP, Scalar, Mat, Mat], None]
""":py:class:`NEP <slepc4py.SLEPc.NEP>` Function callback."""

NEPJacobian = Callable[[NEP, Scalar, Mat], None]
""":py:class:`NEP <slepc4py.SLEPc.NEP>` Jacobian callback."""

# --- SVD ---

SVDStoppingFunction = Callable[[SVD, int, int, int, int], SVD.ConvergedReason]
""":py:class:`SVD <slepc4py.SLEPc.SVD>` stopping test callback."""

SVDMonitorFunction = Callable[[SVD, int, int, ArrayReal, ArrayReal, int], None]
""":py:class:`SVD <slepc4py.SLEPc.SVD>` monitor callback."""

MFNMonitorFunction = Callable[[MFN, int, float], None]
"""`MFN` monitor callback."""

LMEMonitorFunction = Callable[[LME, int, float], None]
"""`LME` monitor callback."""
