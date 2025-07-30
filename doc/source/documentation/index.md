# Documentation

A good starting point for learning SLEPc is to read the Users Manual and then follow the tutorials. See also the video-tutorials. Application programmers can readily begin to use SLEPc from a high level (starting from one of the examples) and then gradually learn more details according to their needs.

## Video-Tutorials

The following video-tutorials have a 10 minute duration each, and include demonstration of actual code executions.

1. [SLEPc: What is it for?](https://media.upv.es/player/?id=440253a0-476d-11e7-9b33-83cdd974e088)
2. [SLEPc: Installation](https://media.upv.es/player/?id=a135c600-476e-11e7-9b33-83cdd974e088)
3. [SLEPc: Compiling and running programs](https://media.upv.es/player/?id=fa4af300-476e-11e7-9b33-83cdd974e088)
4. [SLEPc: Standard eigenvalue problem](https://media.upv.es/player/?id=7c46e660-4770-11e7-9b33-83cdd974e088)
5. [SLEPc: Viewing the solution](https://media.upv.es/player/?id=392e5510-4aa0-11e7-9b33-83cdd974e088)
6. [SLEPc: Spectral transformation](https://media.upv.es/player/?id=3ff6c850-4aa0-11e7-9b33-83cdd974e088)
7. [SLEPc: How to compute different parts of the spectrum](https://media.upv.es/player/?id=46b4db50-4aa0-11e7-9b33-83cdd974e088)
8. [SLEPc: Matrix-free eigenproblems](https://media.upv.es/player/?id=4de53820-4aa0-11e7-9b33-83cdd974e088)

## Printable Documentation

### SLEPc Users Manual

- [[PDF]](../_static/manual/slepc-manual.pdf)


{#str}
### SLEPc Technical Reports (STR)

- STR-1: Orthogonalization Routines in SLEPc - [[PDF]](../_static/reports/str1.pdf)
- STR-2: Single Vector Iteration Methods in SLEPc - [[PDF]](../_static/reports/str2.pdf)
- STR-3: Subspace Iteration in SLEPc - [[PDF]](../_static/reports/str3.pdf)
- STR-4: Arnoldi Methods in SLEPc - [[PDF]](../_static/reports/str4.pdf)
- STR-5: Lanczos Methods in SLEPc - [[PDF]](../_static/reports/str5.pdf)
- STR-6: A Survey of Software for Sparse Eigenvalue Problems - [[PDF]](../_static/reports/str6.pdf)
- STR-7: Krylov-Schur Methods in SLEPc - [[PDF]](../_static/reports/str7.pdf)
- STR-8: Restarted Lanczos Bidiagonalization for the SVD in SLEPc - [[PDF]](../_static/reports/str8.pdf)
- STR-9: Practical Implementation of Harmonic Krylov-Schur - [[PDF]](../_static/reports/str9.pdf)
- STR-10: Davidson Type Subspace Expansions for the Linear Eigenvalue Problem - [[PDF]](../_static/reports/str10.pdf)
- STR-11: Contour Integral Spectrum Slicing Method in SLEPc - [[PDF]](../_static/reports/str11.pdf)

## SLEPc Manual Pages

[Index](../manualpages/singleindex) of all manual pages

Main solver classes:

- [Eigenvalue Problem Solver (EPS)](../manualpages/EPS/index)
- [Singular Value Decomposition (SVD)](../manualpages/SVD/index)
- [Polynomial Eigenvalue Problem (PEP)](../manualpages/PEP/index)
- [Nonlinear Eigenvalue Problem (NEP)](../manualpages/NEP/index)
- [Matrix Function (MFN)](../manualpages/MFN/index)
- [Linear Matrix Equation (LME)](../manualpages/LME/index)

Auxiliary classes and system routines:

- [Spectral Transformation (ST)](../manualpages/ST/index)
- [Direct Solver (DS)](../manualpages/DS/index)
- [Basis Vectors (BV)](../manualpages/BV/index)
- [Mathematical Function (FN)](../manualpages/FN/index)
- [Region (RG)](../manualpages/RG/index)
- [System Routines](../manualpages/Sys/index)

:::{note}
The manual pages are organized in four categories:

* *Beginner* \- Basic usage
* *Intermediate* \- Setting options for algorithms and data structures
* *Advanced* \- Setting more advanced options and customization
* *Developer* \- Interfaces intended primarily for library developers, not for typical applications programmers

:::

## Additional Online Documentation

[Documentation of SLEPc's Python interface (slepc4py)](../slepc4py/index), also available at [readthedocs](https://slepc4py.readthedocs.io/en/stable/)

## PETSc Manual Pages

SLEPc is based on PETSc and therefore users are recommended to use the SLEPc documentation together with the one provided with PETSc.

PETSc Documentation: [[PETSc website]](https://petsc.org/release/docs/)

```{toctree}
:maxdepth: 2
:hidden:

manual/index
hands-on/index
presentations
faq
license
```
