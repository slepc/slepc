# Frequently Asked Questions

```{contents} Table Of Contents
:backlinks: top
:local: true
```

:::{seealso}
{{'[PETSc FAQ](https://petsc.org/{}/faq/)'.format(branch)}}.
:::

## 1.  Where should I send SLEPc bug reports and questions?

Send all maintenance requests to the SLEPc developers via <slepc-maint@upv.es>.

## 2.  Is there a SLEPc users mailing list?

No, but SLEPc-related queries can be posted in the {{'[petsc-users mailing list](https://petsc.org/{}/community/mailing/)'.format(branch)}}.

## 3.  How can I receive announcements of new SLEPc versions?

You can join the slepc-announce mailing list by following the instructions in the [Contact](../contact/mail_list) section. We will update users regarding new major releases through this mailing list.

An alternative is to subscribe to the RSS news feed. In addition to new releases, we will also notify the publication of patches.

## 4.  How should I cite SLEPc?

When writing a scientific paper that makes use of SLEPc, you should cite at least reference [1] in the [list of references](../material/index). In addition, if you use specific SLEPc features (such as computational intervals) that have papers on the list, we suggest citing them as well.

## 5.  Apart from PETSc, is it necessary to install other software to use SLEPc?

No, the only requirement to use SLEPc is to have PETSc installed on your system. Additionally, if you want to have access to eigensolvers not included in SLEPc, then you will have to install other libraries (e.g. ARPACK). See also the comment on linear solvers in FAQ #10 below.

## 6.  I do not see any speedup when using more than one process

:::{important}
Does not apply for version 3.6 or later, see FAQ #10 below.
:::

<span style="color:gray">
Probably you are dealing with a generalized eigenproblem (or a standard eigenproblem with shift-and-invert) and solving the linear systems with the default direct solver. By default, SLEPc uses a direct linear solver via PETSc's `redundant` mechanism, which allows the use of direct solvers in parallel executions but is not a really parallel factorization. In order to get speedup in parallel executions, you need to configure PETSc with a parallel direct linear solver such as MUMPS. For details, see the section "Solution of Linear Systems" in SLEPc's user manual.
</span>

## 7.  Which is the recommended way of learning SLEPc?

Possibly, the best way of learning to use SLEPc is to follow these steps:

  - First of all, get acquainted with PETSc if you are not already familiar with it (see the {{'[PETSc tutorials page](https://petsc.org/{}/tutorials/)'.format(branch)}}).
  - Read through the entire [](manual/index). In a first reading, one may skip the "advanced usage" sections.
  - Follow the steps provided by the [hands-on exercises](hands-on/index), trying the examples in an available SLEPc installation.
  - Use the example programs available in the SLEPc distribution as a basis for your own programs.
  - Use the on-line [manual pages](../manualpages/index) for reference for individual routines.

We also provide several [](#video-tutorials).

## 8.  From 3.0.0 to 3.1 the behavior of shift-and-invert has changed

The shift-and-invert spectral transformation (and Cayley as well) is intended for computing the eigenvalues closest to a given value {math}`\sigma` (the shift). Those eigenvalues closest to the shift become dominant in the transformed spectrum, so in SLEPc 3.0.0 one had to use `EPS_LARGEST_MAGNITUDE` (the default) for this situation. For example (the last option can be omitted because it is the default):

    $ ./ex1 -st_type sinvert -st_shift 3.5 -eps_largest_magnitude

In contrast, in SLEPc 3.1 the approach is to specify the target value directly in EPS (with `EPSSetTarget`) and indicate that we want to compute eigenvalues closest to the target, with `EPS_TARGET_MAGNITUDE`. For example (again, the last option can be omitted):

    $ ./ex1 -st_type sinvert -eps_target 3.5 -eps_target_magnitude

The value of the shift need not be provided because it is taken from the target value.

Note that another difference is that in 3.1 eigenvalues are returned in the correct order, that is, the first one is the closest to the target, and so on.

## 9.  I get an error when retrieving the eigenvector

After the solver has finished, the solution can be retrieved with `EPSGetEigenpair`.  In the `Vr` (and `Vi`) argument, one can pass `NULL` (if the eigenvector is not required), or a _valid_ {external:doc}`Vec` object. This means the vector must have been created, for example with {external:doc}`VecCreate`, {external:doc}`VecDuplicate`, or {external:doc}`MatCreateVecs`, see for instance {{'[ex7.c](https://slepc.upv.es/{}/src/eps/tutorials/ex7.c.html)'.format(branch)}}. The same occurs with analog functions in `SVD`, `PEP`, and `NEP`.

## 10. I get an error when running shift-and-invert in parallel

In 3.6 and later versions, the shift-and-invert spectral transformation defaults to using `preonly`+`lu` for solving linear systems. If you run with more than one MPI process this will fail, unless you use an external package for the parallel LU factorization. This is explained in section [](#sec:lin) in the Users Manual. In previous versions of SLEPc, this would not generate an error since it was using `redundant` rather than plain `lu`.

## 11. Building an application with CMake or pkg_config

SLEPc (and PETSc) provides `pkg_config` files that can be used from makefiles as well as from CMake files. See the discussion at <gl-issue:19>.

## 12. Why is it generally a bad idea to use EPS_SMALLEST_MAGNITUDE?

Krylov methods (and in particular the default SLEPc eigensolver, Krylov-Schur) are good for approximating eigenvalues in the periphery of the spectrum. Assuming an eigenproblem with real eigenvalues only, the use of `EPS_SMALLEST_MAGNITUDE` will be appropriate only if all eigenvalues are either positive or negative. Otherwise, the smallest magnitude eigenvalues lie in the interior of the spectrum, and therefore the convergence will likely be very slow. The usual approach for computing interior eigenvalues is the shift- and-invert spectral transformation (see [](ch:st) in the Users Manual). Hence, instead of `-eps_smallest_magnitude` one would generally prefer `-st_type sinvert -eps_target 0`.

## 13. Creating a sparse matrix gets terribly slow when I increase the matrix size

Matrix preallocation is extremely important, especially for large matrices.  See the {{'[chapter on matrices in the PETSc users manual](https://petsc.org/{}/manual/mat/#preallocation-of-memory-for-sequential-aij-sparse-matrices)'.format(branch)}}.

:::{note}
Since PETSc version 3.19 the {external:doc}`Mat` data structures have been changed so that the performance is reasonably good without preallocation.
:::

## 14. conda-forge: how to install slepc4py with complex scalars

The following instructions can be followed to install the conda-forge variant of slepc4py with scalar type complex.

> Download Miniforge
>
> ```{code} console
> $ wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
> ```
>
> Install miniforge in $HOME/miniforge
>
> ```{code} console
> $ bash Miniforge3-Linux-x86_64.sh -b -p ~/miniforge
> ```
>
> Activate Miniforge base environment
> ```{code} console
> $ source ~/miniforge/bin/activate
> ```
>
> Create a new test environment with a starting package specification set. Package specification is `<package>=<version>==<build_string>`. To request the latest version, we can use `*`
>
> For {petsc|slepc}[4py], the build string starts with either 'real' or 'complex', depending on the configured scalar type. The 'real' variants have a build number offset by 100, so they take precedence if a specific `<build_string>` is not requested.
>
> Long story short, we ask conda to create a new environment initally installing the complex variant of slepc4py latest version with the followng command
>
> ```{code} console
> $ conda create --name testenv slepc4py=*=*complex*
> ```
>
> Activate the test environment
>
> ```{code} console
> $ conda activate testenv
> ```
>
> Verify we are running the complex variant
> ```{code} console
> $ python -c '
>   from slepc4py import SLEPc
>   from petsc4py import PETSc
>   print(PETSc.ScalarType)
> '
> ```

## 15. spack: how to install slepc4py with complex scalars

Make sure you select the `petsc+complex` variant:

```{code} console
$ spack install py-slepc4py ^petsc+complex
```

This will install PETSc with complex scalars, together with SLEPc as well as petsc4py and slepc4py. Before that, you can also do `spack spec py-slepc4py ^petsc+complex` to check what it is going to install.

## 16. Eigenvectors have nonzero imaginary part

A real symmetric matrix has real eigenvectors, but when building SLEPc with complex scalars the computed eigenvectors have nonzero imaginary part. The rationale is the following. In real scalars, if $x$ is a unit-norm eigenvector then $-x$ is also a valid eigenvector. In complex scalars, if $x$ is a unit-norm eigenvector then $\alpha x$ is also a valid eigenvector, where $\alpha$ is a generalized sign, i.e., $\alpha=\exp(\theta j)$ for any $\theta$. So if one wants the imaginary part to be zero, the eigenvectors returned by SLEPc must be normalized a posteriori, as is done for example in {{'[ex20.c](https://slepc.upv.es/{}/src/eps/tutorials/ex20.c.html)'.format(branch)}} (or the equivalent python example `ex7.py`). SLEPc does not know if the input matrix is real or complex, so it cannot normalize the vectors internally.

Note that the simple scaling strategy shown in those examples will not be sufficient in case of degenerate eigenvalues, i.e., eigenvalues with multiplicity larger than one.
