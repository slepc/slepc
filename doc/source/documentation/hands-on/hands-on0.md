Tutorials

# Hello World

This exercise shows how to build and run a simple example program with SLEPc.

:::{note}
The description below related to directories and the use of the `PETSC_ARCH` variable will be different in the case of a prefix-based installation.
:::

## Compiling

SLEPc needs the following environment variables to be set:

`SLEPC_DIR` \- the location of SLEPc

`PETSC_DIR` \- the location of PETSc

`PETSC_ARCH` \- the architecture being used

Make sure that you have them correctly set.

Like in PETSc, a makefile is necessary to compile a SLEPc program. Paste this simple example into a file named `makefile` in your working directory:

```{code} make
hello: hello.o
	-${CLINKER} -o hello hello.o ${SLEPC_SYS_LIB}
	${RM} hello.o

include ${SLEPC_DIR}/lib/slepc/conf/slepc_common
```

:::{note}
In the above text, the blank space in the 2nd and 3rd lines is a tab.
:::

Also place the following source code into a file named "hello.c" in the same directory:

```{code} c
static char help[] = "Simple Hello World example program in SLEPc\n";

#include <slepcsys.h>

int main( int argc, char **argv )
{
  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,(char*)0,help));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Hello world\n"));
  PetscCall(SlepcFinalize());
  return 0;
}
```

Compile the program with the following command:

```{code} console
$ make hello
```

## Source Code Details

Examine the source code of the sample program. The following comments will help you understand the code thoroughly.

**Header File** : All SLEPc programs must include a header file with all the necessary definitions. In this simple example the file {{'[slepcsys.h](https://slepc.upv.es/{}/include/slepcsys.h.html)'.format(branch)}} (base SLEPc header) is enough because no solver components such as EPS are being used.

:::{note}
SLEPc header files automatically include some PETSc header files.
:::

**Library Initialization** : All SLEPc programs must begin with a call to [SlepcInitialize](../../manualpages/Sys/SlepcInitialize), which in turn initializes PETSc and MPI. Similarly, at the end of the program [SlepcFinalize](../../manualpages/Sys/SlepcFinalize) must be called for library cleanup.

**Input/Output** : In this example, we do input/output via a call to a PETSc function, [PetscPrintf](https://petsc.org/release/manualpages/Sys/PetscPrintf).  Remember that in parallel programs input/output cannot be done simply via C standard library functions. Note that in SLEPc programs we can freely use any PETSc function.

**Error Checking** : All SLEPc routines return an integer indicating whether an error has occurred during the call. The PETSc macro [PetscCall](https://petsc.org/release/manualpages/Sys/PetscCall) checks the return value and calls the PETSc error handler upon error detection. All function calls should be wrapped around `PetscCall` to enable a complete error traceback.

## Running the Program

SLEPc programs are executed as any other MPI program. Note that this typically differs from one system to another. To run the program with only one process, in some systems you can launch it as a normal program:

```{code} console
$ ./hello
```

but in other systems this would not work. Check the documentation of your system. Standard MPI implementations provide the `mpiexec` command to launch the applications

```{code} console
$ mpiexec -n 4 ./hello
```

In SLEPc (and PETSc) there are a lot of options (run-time parameters) to control program behavior. These options are usually equivalent to function calls, so the user can test its effect without changing the source code (this will be illustrated in the next exercises). To show which options are available in a program use:

```{code} console
$ ./hello -help
```

## Support for Debugging and Complex Numbers

The support for debugging capabilities, complex scalar arithmetic, and other features is managed by SLEPc and PETSc by means of different _architectures_ , represented by different values of the `PETSC_ARCH` variable. In a given system, you can typically find several versions of SLEPc and PETSc, each of them built with different configuration options. For instance, suppose the following values are available:

  * `arch-linux-gnu-c-opt`: built with compiler optimization
  * `arch-linux-gnu-c-debug`: built with debugging support
  * `arch-linux-gnu-c-opt-complex`: optimized with complex scalars
  * `arch-linux-gnu-c-debug-complex`: debug with complex scalars

:::{note}
In order to learn about the particular architectures available in your system, type `ls $SLEPC_DIR`. There should be a subdirectory for each allowed value of `PETSC_ARCH`.
:::

When using an architecture with support for complex scalars, all scalar values are complex instead of real. Try compiling the example program for complex numbers:

```{code} console
$ make PETSC_ARCH=arch-linux-gnu-c-opt-complex hello
```

When using the debug versions some options are available to support debugging.  For example

```{code} console
$ ./hello -start_in_debugger
```

opens the program in a debugger stopped at the [SlepcInitialize](../../manualpages/Sys/SlepcInitialize) function.

Other useful options are: `-info` to get informative messages about progress of the calculations, `-malloc_info` to print memory usage at end of run, `-log_trace [filename]` to get a full trace of the execution (in a file), `-malloc_dump` to list memory blocks not freed at the end of the program, and `-log_view` to get a summary including performance results.
