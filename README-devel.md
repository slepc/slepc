
Notes for SLEPc developers and contributors
===========================================

This document is intended for pointing out the differences between SLEPc and PETSc that may be relevant for developers and contributors.

Additional information can be found in the [Developers section of PETSc's web site](https://petsc.org/release/developers).

Build system
------------

- The SLEPc `configure` script is independent of the PETSc counterpart, there is no overlap. It gets some PETSc configuration data from files `petscvariables` and `petscconf.h`.
- Makefiles in SLEPc follow a similar structure as in PETSc. When possible, relevant PETSc makefiles are included from SLEPc makefiles, in particular, files `rules`, `rules.util` and `variables` in `${PETSC_DIR}/lib/petsc/conf` are included.

Continuous integration
----------------------

- The SLEPc project uses its own gitlab-runners.
- Pipelines are generated automatically when pushing to a merge request, as in PETSc. In that case, a merge against the target branch is done (merged results pipeline).
- An alternative is to manually run a 'detached' pipeline (without merge) with the `New pipeline` button in the `Pipelines` item under the `Build` menu.
- It is possible to select a PETSc branch to be used in the pipeline, by using `CI_PETSC_BRANCH` as the variable key and the branch name as the variable value. In the case of a manual pipeline, this is done in the `New pipeline` screen. In the case of automatic pipelines, one has to enter the `pause-for-approval` job before pushing the play button, and insert the variable there.
- The test harness is run with `DIFF_NUMBERS` enabled by default, as opposed to PETSc. When adding a new test, all the output must match, including floating point numbers. Use filters to remove potentially problematic values such as small residual norms.
- Using filters in tests is preferred to adding `alt` output files.

Code
----

- In SLEPc, code style is not enforced via `clang-format`. Still, most coding conventions should be followed. They are the same as the ones from PETSc prior to the switch to `clang-format`, see the [PETSc Style and Usage Guide](https://petsc.org/release/developers/style/).

