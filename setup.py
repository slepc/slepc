#!/usr/bin/env python3

"""
SLEPc: Scalable Library for Eigenvalue Problem Computations
===========================================================

SLEPc is a software library for the solution of large scale sparse
eigenvalue problems on parallel computers. It is an extension of PETSc
and can be used for either standard or generalized eigenproblems, with
real or complex arithmetic. It can also be used for computing a
partial SVD of a large, sparse, rectangular matrix, and to solve
nonlinear eigenvalue problems.

.. note::

   To install the ``PETSc``, ``SLEPc``, ``petsc4py``, and ``slepc4py`` packages
   (``mpi4py`` is optional but highly recommended) use::

     $ python -m pip install numpy mpi4py
     $ python -m pip install petsc petsc4py
     $ python -m pip install slepc slepc4py

.. tip::

  You can also install the in-development versions with::

    $ python -m pip install Cython numpy mpi4py
    $ python -m pip install --no-deps https://gitlab.com/petsc/petsc/-/archive/main/petsc-main.tar.gz
    $ python -m pip install --no-deps https://gitlab.com/slepc/slepc/-/archive/main/slepc-main.tar.gz

  Provide any ``SLEPc`` ./configure options using the environmental variable ``SLEPC_CONFIGURE_OPTIONS``.

"""

import re
import os
import sys
import shlex
import shutil
from setuptools import setup
from setuptools.command.install import install as _install
try:
    from setuptools.command.bdist_wheel import bdist_wheel as _bdist_wheel
except ImportError:
    from wheel.bdist_wheel import bdist_wheel as _bdist_wheel
from distutils import log

init_py = """\
# Author:  SLEPc Team
# Contact: slepc-maint@upv.es

def get_slepc_dir():
    import os
    return os.path.dirname(__file__)


def get_config():
    conf = {}
    conf['SLEPC_DIR'] = get_slepc_dir()
    return conf
"""

main_py = """\
# Author:  SLEPc Team
# Contact: slepc-maint@upv.es

if __name__ == "__main__":
    import sys
    if "--prefix" in sys.argv:
        from . import get_slepc_dir
        print(get_slepc_dir())
        del get_slepc_dir
    del sys
"""

metadata = {
    'provides' : ['slepc'],
    'zip_safe' : False,
}

CONFIGURE_OPTIONS = []


def bootstrap():
    # Set SLEPC_DIR
    SLEPC_DIR  = os.path.abspath(os.getcwd())
    os.environ['SLEPC_DIR']  = SLEPC_DIR
    # Check PETSC_DIR/PETSC_ARCH
    PETSC_DIR  = os.environ.get('PETSC_DIR',  "")
    PETSC_ARCH = os.environ.get('PETSC_ARCH', "")
    PETSC_ARCH_DIR = os.path.join(PETSC_DIR, PETSC_ARCH)
    if not (PETSC_DIR and os.path.isdir(PETSC_DIR)):
        PETSC_DIR = None
        os.environ.pop('PETSC_DIR', None)
        PETSC_ARCH = None
        os.environ.pop('PETSC_ARCH', None)
    elif not (PETSC_ARCH and os.path.isdir(PETSC_ARCH_DIR)):
        PETSC_ARCH = None
        os.environ.pop('PETSC_ARCH', None)

    # Generate package __init__.py and __main__.py files
    pkgdir = os.path.join('config', 'pypi')
    os.makedirs(pkgdir, exist_ok=True)
    for pyfile, contents in (
        ('__init__.py', init_py),
        ('__main__.py', main_py),
    ):
        with open(os.path.join(pkgdir, pyfile), 'w') as fh:
            fh.write(contents)

    # Configure options
    options = os.environ.get('SLEPC_CONFIGURE_OPTIONS', '')
    CONFIGURE_OPTIONS.extend(shlex.split(options))
    #
    if not PETSC_DIR:
        vstr = version()
        x, y = tuple(map(int, vstr.split('.')[:2]))
        dev = '.dev0' if '.dev' in vstr else ''
        reqs = ">=%s.%s%s,<%s.%s" % (x, y, dev, x, y+1)
        metadata['setup_requires'] = ['petsc'+reqs]
        metadata['install_requires'] = ['petsc'+reqs]


def get_petsc_dir():
    PETSC_DIR = os.environ.get('PETSC_DIR')
    if PETSC_DIR:
        return PETSC_DIR
    try:
        import petsc
        PETSC_DIR = petsc.get_petsc_dir()
    except ImportError:
        log.warn("PETSC_DIR not specified")
        PETSC_DIR = os.path.join(os.path.sep, 'usr', 'local', 'petsc')
    return PETSC_DIR


def get_petsc_arch():
    PETSC_ARCH = os.environ.get('PETSC_ARCH', "")
    return PETSC_ARCH


def config(prefix, dry_run=False):
    log.info('SLEPc: configure')
    options = [
        '--prefix=' + prefix,
        ]
    options.extend(CONFIGURE_OPTIONS)
    #
    log.info('configure options:')
    for opt in options:
        log.info(' '*4 + opt)
    # Run SLEPc configure
    if dry_run:
        return
    os.environ['PETSC_DIR'] = get_petsc_dir()
    os.environ['PETSC_ARCH'] = get_petsc_arch()
    python = sys.executable
    command = [python, './configure'] + options
    status = os.system(" ".join(command))
    if status != 0:
        raise RuntimeError(status)
    # Fix SLEPc configuration
    using_build_backend = any(
        os.environ.get(prefix + '_BUILD_BACKEND')
        for prefix in ('_PYPROJECT_HOOKS', 'PEP517')
    )
    if using_build_backend:
        pdir = os.environ['SLEPC_DIR']
        parch = os.environ['PETSC_ARCH']
        if not parch:
            makefile = os.path.join(pdir, 'lib', 'slepc', 'conf', 'slepcvariables')
            with open(makefile, 'r') as mfile:
                contents = mfile.readlines()
            for line in contents:
                if line.startswith('PETSC_ARCH'):
                    parch = line.split('=')[1].strip()
                    break
        include = os.path.join(pdir, parch, 'include')
        for filename in (
            'slepcconf.h',
        ):
            filename = os.path.join(include, filename)
            with open(filename, 'r') as old_fh:
                contents = old_fh.read()
            contents = contents.replace(prefix, '${SLEPC_DIR}')
            with open(filename, 'w') as new_fh:
                new_fh.write(contents)


def build(dry_run=False):
    log.info('SLEPc: build')
    # Run SLEPc build
    if dry_run:
        return
    PETSC_ARCH = get_petsc_arch()
    if PETSC_ARCH:
        PETSC_ARCH = 'PETSC_ARCH=' + PETSC_ARCH
    make = shutil.which('make')
    command = [make, 'all', 'PETSC_DIR='+get_petsc_dir(), PETSC_ARCH]
    status = os.system(" ".join(command))
    if status != 0:
        raise RuntimeError(status)


def install(dry_run=False):
    log.info('SLEPc: install')
    # Run SLEPc install
    if dry_run:
        return
    PETSC_ARCH = get_petsc_arch()
    if PETSC_ARCH:
        PETSC_ARCH = 'PETSC_ARCH=' + PETSC_ARCH
    make = shutil.which('make')
    command = [make, 'install', 'PETSC_DIR='+get_petsc_dir(), PETSC_ARCH]
    status = os.system(" ".join(command))
    if status != 0:
        raise RuntimeError(status)


class context(object):
    def __init__(self):
        self.sys_argv = sys.argv[:]
        self.wdir = os.getcwd()
    def enter(self):
        del sys.argv[1:]
        pdir = os.environ['SLEPC_DIR']
        os.chdir(pdir)
        return self
    def exit(self):
        sys.argv[:] = self.sys_argv
        os.chdir(self.wdir)


class cmd_install(_install):

    def initialize_options(self):
        _install.initialize_options(self)

    def finalize_options(self):
        _install.finalize_options(self)
        self.install_lib = self.install_platlib
        self.install_libbase = self.install_lib
        self.old_and_unmanageable = True

    def run(self):
        root_dir = os.path.abspath(self.install_lib)
        prefix = os.path.join(root_dir, 'slepc')
        #
        ctx = context().enter()
        try:
            config(prefix, self.dry_run)
            build(self.dry_run)
            install(self.dry_run)
        finally:
            ctx.exit()
        #
        self.outputs = []
        for dirpath, _, filenames in os.walk(prefix):
            for fn in filenames:
                self.outputs.append(os.path.join(dirpath, fn))
        #
        _install.run(self)

    def get_outputs(self):
        outputs = getattr(self, 'outputs', [])
        outputs += _install.get_outputs(self)
        return outputs


class cmd_bdist_wheel(_bdist_wheel):

    def finalize_options(self):
        super().finalize_options()
        self.root_is_pure = False
        self.build_number = None

    def get_tag(self):
        plat_tag = super().get_tag()[-1]
        return (self.python_tag, "none", plat_tag)


def version():
    version_re = {
        'major'  : re.compile(r"#define\s+SLEPC_VERSION_MAJOR\s+(\d+)"),
        'minor'  : re.compile(r"#define\s+SLEPC_VERSION_MINOR\s+(\d+)"),
        'micro'  : re.compile(r"#define\s+SLEPC_VERSION_SUBMINOR\s+(\d+)"),
        'release': re.compile(r"#define\s+SLEPC_VERSION_RELEASE\s+(\d+)"),
        }
    slepcversion_h = os.path.join('include','slepcversion.h')
    data = open(slepcversion_h, 'r').read()
    major = int(version_re['major'].search(data).groups()[0])
    minor = int(version_re['minor'].search(data).groups()[0])
    micro = int(version_re['micro'].search(data).groups()[0])
    release = int(version_re['release'].search(data).groups()[0])
    if release:
        v = "%d.%d.%d" % (major, minor, micro)
    else:
        v = "%d.%d.0.dev%d" % (major, minor+1, 0)
    return v


def tarball():
    VERSION = version()
    if '.dev' in VERSION:
        return None
    return ('https://slepc.upv.es/download/distrib/'
            'slepc-%s.tar.gz#egg=slepc-%s' % (VERSION, VERSION))


description = __doc__.split('\n')[1:-1]
del description[1:3]

classifiers = """
Development Status :: 5 - Production/Stable
Intended Audience :: Developers
Intended Audience :: Science/Research
License :: OSI Approved :: BSD License
Operating System :: POSIX
Programming Language :: C
Programming Language :: C++
Programming Language :: Fortran
Programming Language :: Python
Topic :: Scientific/Engineering
Topic :: Software Development :: Libraries
"""

bootstrap()
setup(
    name='slepc',
    version=version(),
    description=description.pop(0),
    long_description='\n'.join(description),
    long_description_content_type='text/x-rst',
    classifiers=classifiers.split('\n')[1:-1],
    keywords = ['SLEPc', 'PETSc', 'MPI'],
    platforms=['POSIX'],
    license='BSD-2-Clause',

    url='https://slepc.upv.es/',
    download_url=tarball(),

    author='SLEPc Team',
    author_email='slepc-maint@upv.es',
    maintainer='Lisandro Dalcin',
    maintainer_email='dalcinl@gmail.com',

    packages=['slepc'],
    package_dir={'slepc': 'config/pypi'},
    cmdclass={
        'install': cmd_install,
        'bdist_wheel': cmd_bdist_wheel,
    },
    **metadata
)
