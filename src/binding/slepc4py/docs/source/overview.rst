Overview
========

*SLEPc for Python* (slepc4py) is a Python package that provides
convenient access to the functionality of SLEPc.

SLEPc [1]_, [2]_ implements algorithms and tools for the numerical
solution of large, sparse eigenvalue problems on parallel
computers. It can be used for linear eigenvalue problems in either
standard or generalized form, with real or complex arithmetic.
It can also be used for computing a partial SVD of a large, sparse,
rectangular matrix, and to solve nonlinear eigenvalue problems
(polynomial or general). Additionally, SLEPc provides solvers for
the computation of the action of a matrix function on a vector.

SLEPc is intended for computing a subset of the spectrum of a matrix
(or matrix pair). One can for instance approximate the largest
magnitude eigenvalues, or the smallest ones, or even those eigenvalues
located near a given region of the complex plane. Interior eigenvalues
are harder to compute, so SLEPc provides different methodologies. One
such method is to use a spectral transformation. Cheaper alternatives
are also available.

.. [1] J. E. Roman, C. Campos, L. Dalcin, E. Romero, A. Tomas.
   SLEPc Users Manual. DSIC-II/24/02 - Revision 3.23.
   D. Sistemas Informaticos y Computacion, Universitat Politecnica de
   Valencia. 2025.

.. [2] Vicente Hernandez, Jose E. Roman and Vicente Vidal.
   SLEPc: A Scalable and Flexible Toolkit for the Solution of
   Eigenvalue Problems, ACM Trans. Math. Softw. 31(3), pp. 351-362,
   2005.

.. include:: links.txt


Features
--------

Currently, the following types of eigenproblems can be addressed:

* Standard eigenvalue problem, *Ax=kx*, either for Hermitian or
  non-Hermitian matrices.

* Generalized eigenvalue problem, *Ax=kBx*, either Hermitian
  positive-definite or not.

* Partial singular value decomposition of a rectangular matrix,
  *Au=sv*.

* Polynomial eigenvalue problem, *P(k)x=0*.

* Nonlinear eigenvalue problem, *T(k)x=0*.

* Computing the action of a matrix function on a vector, *w=f(alpha A)v*.

For the linear eigenvalue problem, the following methods are available:

* Krylov eigensolvers, particularly Krylov-Schur, Arnoldi, and
  Lanczos.

* Davidson-type eigensolvers, including Generalized Davidson and
  Jacobi-Davidson.

* Subspace iteration and single vector iterations (inverse iteration,
  RQI).

* Conjugate gradient for the minimization of the Rayleigh quotient.

* A contour integral solver.

For singular value computations, the following alternatives can be
used:

* Use an eigensolver via the cross-product matrix *A'A* or the cyclic
  matrix *[0 A; A' 0]*.

* Explicitly restarted Lanczos bidiagonalization.

* Implicitly restarted Lanczos bidiagonalization (thick-restart
  Lanczos).

For polynomial eigenvalue problems, the following methods are available:

* Use an eigensolver to solve the generalized eigenvalue problem
  obtained after linearization.

* TOAR and Q-Arnoldi, memory efficient variants of Arnoldi for polynomial
  eigenproblems.

For general nonlinear eigenvalue problems, the following methods can be used:

* Solve a polynomial eigenproblem obtained via polynomial interpolation.

* Rational interpolation and linearization (NLEIGS).

* Newton-type methods such as SLP or RII.

Computation of interior eigenvalues is supported by means of the
following methodologies:

* Spectral transformations, such as shift-and-invert. This technique
  implicitly uses the inverse of the shifted matrix *(A-tI)* in order
  to compute eigenvalues closest to a given target value, *t*.

* Harmonic extraction, a cheap alternative to shift-and-invert that
  also tries to approximate eigenvalues closest to a target, *t*, but
  without requiring a matrix inversion.

Other remarkable features include:

* High computational efficiency, by using NumPy and SLEPc under the
  hood.

* Data-structure neutral implementation, by using efficient sparse
  matrix storage provided by PETSc. Implicit matrix representation is
  also available by providing basic operations such as matrix-vector
  products as user-defined Python functions.

* Run-time flexibility, by specifying numerous setting at the command
  line.

* Ability to do the computation in parallel.


Components
----------

SLEPc provides the following components, which are mirrored by slepc4py
for its use from Python. The first five components are solvers for
different classes of problems, while the rest can be considered
auxiliary object.

:EPS: The Eigenvalue Problem Solver is the component that provides all
      the functionality necessary to define and solve an
      eigenproblem. It provides mechanisms for completely specifying
      the problem: the problem type (e.g. standard symmetric), number
      of eigenvalues to compute, part of the spectrum of
      interest. Once the problem has been defined, a collection of
      solvers can be used to compute the required solutions.  The
      behaviour of the solvers can be tuned by means of a few
      parameters, such as the maximum dimension of the subspace to be
      used during the computation.

:SVD: This component is the analog of EPS for the case of Singular
      Value Decompositions. The user provides a rectangular matrix and
      specifies how many singular values and vectors are to be
      computed, whether the largest or smallest ones, as well as some
      other parameters for fine tuning the computation. Different
      solvers are available, as in the case of EPS.

:PEP: This component is the analog of EPS for the case of Polynomial
      Eigenvalue Problems. The user provides the coefficient matrices of
      the polynomial. Several parameters can be specified, as in
      the case of EPS. It is also possible to indicate whether the
      problem belongs to a special type, e.g., symmetric or gyroscopic.

:NEP: This component covers the case of general nonlinear eigenproblems,
      T(lambda)x=0. The user provides the parameter-dependent matrix T
      via the split form or by means of callback functions.

:MFN: This component provides the functionality for computing the action
      of a matrix function on a vector. Given a matrix A and a vector b,
      the call MFNSolve(mfn,b,x) computes x=f(A)b, where f is a function
      such as the exponential.

:ST:  The Spectral Transformation is a component that provides
      convenient implementations of common spectral
      transformations. These are simple transformations that map
      eigenvalues to different positions, in such a way that
      convergence to wanted eigenvalues is enhanced. The most common
      spectral transformation is shift-and-invert, that allows for the
      computation of eigenvalues closest to a given target value.

:BV:  This component encapsulates the concept of a set of Basis Vectors
      spanning a vector space. This component provides convenient access
      to common operations such as orthogonalization of vectors. The
      BV component is usually not required by end-users.

:DS:  The Dense System (or Direct Solver) component, used internally to
      solve dense eigenproblems of small size that appear in the course
      of iterative eigensolvers.

:FN:  A component used to define mathematical functions. This is required
      by the end-user for instance to define function T(.) when solving
      nonlinear eigenproblems with NEP in split form.

:RG:  A component used to define a region of the complex plane such as an
      ellipse or a rectangle. This is required by end-users in some cases
      such as contour-integral eigensolvers.

In addition to the above components, some extra functionality is provided
in the :Sys: and :Util: sections.
