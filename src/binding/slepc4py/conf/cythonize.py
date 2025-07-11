#!/usr/bin/env python3
"""Run Cython with custom options."""
import os
import sys

appdir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, appdir)

import cyautodoc  # noqa: F401,E402
from Cython.Compiler.Main import main as cython_main  # noqa: E402


def cythonize(args=None):
    """Run `cython --3str --cleanup 3 <args>...`."""
    if args is None:
        argv = sys.argv[:]
    else:
        argv = [os.path.abspath(__file__)] + list(args)

    if '--cleanup' not in argv:
        argv[1:1] = ['--cleanup', '3']
    if '--3str' not in argv:
        argv[1:1] = ['--3str']

    cwd = os.getcwd()
    sys_argv = sys.argv[:]
    try:
        sys.argv[:] = argv
        cython_main(command_line=1)
        return 0
    except SystemExit as exc:
        return exc.code
    finally:
        os.chdir(cwd)
        sys.argv[:] = sys_argv


def main():
    """Entry-point to run Cython with custom options."""
    args = sys.argv[1:]
    if not args:
        topdir = os.path.dirname(appdir)
        srcdir = os.path.join(topdir, 'src')
        source = os.path.join('slepc4py', 'SLEPc.pyx')
        target = os.path.join('slepc4py', 'SLEPc.c')
        args += ['--working', srcdir]
        args += [source, '--output-file', target]
    sys.exit(cythonize(args))


if __name__ == "__main__":
    main()
