# Author:  Lisandro Dalcin
# Contact: dalcinl@gmail.com
"""Typing support."""

from __future__ import annotations  # novermin
from typing import (  # novermin
    Callable,
    Sequence,
    Literal,
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
]

# --- EPS ---

EPSStoppingFunction = Callable[[EPS, int, int, int, int], EPS.ConvergedReason]
"""`EPS` stopping test callback."""

EPSArbitraryFunction = Callable[[Scalar, Scalar, Vec, Vec, Scalar, Scalar], [Scalar, Scalar]]
"""`EPS` arbitrary selection callback."""

EPSEigenvalueComparison = Callable[[Scalar, Scalar, Scalar, Scalar], int]
"""`EPS` eigenvalue comparison callback."""

EPSMonitorFunction = Callable[[EPS, int, int, ScalarArray, ScalarArray, RealArray, int], None]
"""`EPS` monitor callback."""

# --- PEP ---

PEPStoppingFunction = Callable[[PEP, int, int, int, int], PEP.ConvergedReason]
"""`PEP` stopping test callback."""

PEPMonitorFunction = Callable[[PEP, int, int, ScalarArray, ScalarArray, RealArray, int], None]
"""`PEP` monitor callback."""

# --- NEP ---

NEPStoppingFunction = Callable[[NEP, int, int, int, int], NEP.ConvergedReason]
"""`NEP` stopping test callback."""

NEPMonitorFunction = Callable[[NEP, int, int, ScalarArray, ScalarArray, RealArray, int], None]
"""`NEP` monitor callback."""

NEPFunction = Callable[[NEP, Scalar, Mat, Mat], None]
"""`NEP` Function callback."""

NEPJacobian = Callable[[NEP, Scalar, Mat], None]
"""`NEP` Jacobian callback."""

# --- SVD ---

SVDStoppingFunction = Callable[[SVD, int, int, int, int], SVD.ConvergedReason]
"""`SVD` stopping test callback."""

SVDMonitorFunction = Callable[[SVD, int, int, RealArray, RealArray, int], None]
"""`SVD` monitor callback."""

