# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

#from os import path
import sys
import os
import re
import subprocess
from datetime import datetime
import time
import shutil
import sphobjinv

#print('cwd:', os.getcwd())
sys.path.append(os.getcwd())
sys.path.append(os.path.abspath('./ext'))

import build_manpages_c2html
import update_htmlmap_links
import make_links_relative
import fix_man_page_edit_links
import add_man_page_redirects

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'SLEPc'
author = 'SLEPc Team'
version = 'main'
release = 'main'

with open(os.path.join('../..', 'include', 'slepcversion.h'),'r') as version_file:
    buf = version_file.read()
    slepc_release_flag = re.search(' SLEPC_VERSION_RELEASE[ ]*([0-9]*)',buf).group(1)
    major_version      = re.search(' SLEPC_VERSION_MAJOR[ ]*([0-9]*)',buf).group(1)
    minor_version      = re.search(' SLEPC_VERSION_MINOR[ ]*([0-9]*)',buf).group(1)
    subminor_version   = re.search(' SLEPC_VERSION_SUBMINOR[ ]*([0-9]*)',buf).group(1)

    version = '.'.join([major_version, minor_version])
    if slepc_release_flag == '0':
        release = '.'.join([major_version,minor_version]) + '-dev'
    else:
        release = '.'.join([major_version,minor_version,subminor_version])


release_date = subprocess.check_output(['git',
                                 'for-each-ref',
                                 '--format="%(creatordate:format:%B), %(creatordate:format:%Y)"',
                                 'refs/tags/v{}.0'.format(version)]
                                ).strip().decode('utf-8').strip('"')

release_year = subprocess.check_output(['git',
                                 'for-each-ref',
                                 '--format="%(creatordate:format:%Y)"',
                                 'refs/tags/v{}.0'.format(version)]
                                ).strip().decode('utf-8').strip('"')

if release_date == '':
    release_date = datetime.strftime(datetime.now(), '%B, %Y')
if release_year == '':
    release_year = datetime.strftime(datetime.now(), '%Y')

try:
    git_ref = subprocess.check_output(["git", "rev-parse", "HEAD"]).rstrip()
    git_ref_release = subprocess.check_output(["git", "rev-parse", "origin/release"]).rstrip()
    edit_branch = "release" if git_ref == git_ref_release else "main"
except subprocess.CalledProcessError:
    print("WARNING: determining branch for page edit links failed")
    edit_branch = "main"

copyright = '2002-{}, Universitat Politecnica de Valencia, Spain'.format(
        release_year)

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

needs_sphinx='7.3.7'

nitpicky = True # warn about all references where the target cannot be found

extensions = [
        'myst_parser',         # use structured markdown
        'sphinx.ext.duration',
        'sphinx_copybutton',   # Add a copy button to your code blocks
        'sphinx_design',       # tabs
        'sphinx_togglebutton', # dropdown box
        'sphinxcontrib.bibtex',
        'sphinxcontrib.rsvgconverter', # svg to pdf (manual)
        'sphinx.ext.intersphinx',
        'html5_petsc',
        ]

# intersphinx
intersphinx_mapping = {
    'petsc': ('https://petsc.org/release/', None),
}

def _mangle_petsc_intersphinx(branch):
    """Preprocess the keys in PETSc's intersphinx inventory.

    PETSc have intersphinx keys of the form:

        manualpages/Vec/VecShift

    This function downloads their object inventory and strips the leading path
    elements so that references to PETSc names actually resolve."""

    website = intersphinx_mapping['petsc'][0].partition('/release/')[0]
    doc_url = f'{website}/{branch}/'
    inventory_url = f'{doc_url}objects.inv'
    print('Using PETSC inventory from ' + inventory_url)
    inventory = sphobjinv.Inventory(url=inventory_url)
    print(inventory)

    for obj in inventory.objects:
        if obj.name.startswith('manualpages'):
            name = obj.name.split('/')[-1]
            if not name == 'index':
                obj.name = obj.name.split('/')[-1]

    new_inventory_filename = 'petsc_objects.inv'
    sphobjinv.writebytes(
        new_inventory_filename, sphobjinv.compress(inventory.data_file(contract=True))
    )
    intersphinx_mapping['petsc'] = (doc_url, new_inventory_filename)

_mangle_petsc_intersphinx(edit_branch)

myst_links_external_new_tab = True # open all external links in new tabs

bibtex_bibfiles = ['slepc.bib']
##########################################################################
import pybtex.plugin
from pybtex.style.formatting.alpha import Style as AlphaStyle
from pybtex.style.labels.alpha import LabelStyle as AlphaLabelStyle

class SLEPcLabelStyle(AlphaLabelStyle):
    def format_labels(self, sorted_entries):
        for entry in sorted_entries:
            yield entry.key

class SLEPcStyle(AlphaStyle):
    default_label_style = SLEPcLabelStyle

pybtex.plugin.register_plugin('pybtex.style.formatting', 'slepcstyle', SLEPcStyle)
##########################################################################
#bibtex_default_style = 'unsrtalpha' # alpha, plain, unsrt, unsrtalpha, and custom
bibtex_default_style = 'slepcstyle'
bibtex_reference_style = 'author_year' # label, author_year, super

# prevents incorrect WARNING: duplicate citation for key "xxxx" warnings
suppress_warnings = ['bibtex.duplicate_citation']

myst_enable_extensions = [
        'substitution', # enable substitutions
        'colon_fence',  # allow also the use of ::: delimiters to denote directives
        'smartquotes',  # convert quotations to their opening/closing variants
#        'replacements', # convert some common typographic texts (c) copyright
        'deflist', # definition lists
        'attrs_block',  # parsing of block attributes XXX user manual
        'attrs_inline', # parsing of inline attributes
        'dollarmath',   # use of $$ latex math
#        'html_image',
        'fieldlist',
        ]

myst_dmath_double_inline = True
myst_dmath_allow_labels = True # the default
myst_dmath_allow_space = True
myst_dmath_allow_digits=False

myst_substitutions = {
            'release_date': release_date,
            'release_year': release_year,
            'branch': edit_branch,
            }

myst_url_schemes = {
        "http": None,
        "https": None,
        "doi": {
            "url": "https://doi.org/{{path}}",
            "title": "[DOI]",
            },
        "gl-issue": {
            "url": "https://gitlab.com/slepc/slepc/-/issues/{{path}}",
            "title": "Issue #{{path}}",
            "classes": ["gitlab"],
            },
        }

copybutton_prompt_text = '$ ' # the prompt is not copied

templates_path = ['_templates']
exclude_patterns = ['_static/README.md']

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'pydata_sphinx_theme'

#html_theme_options = {}

html_static_path = ['_static']

html_favicon = '_static/images/favicon-slepc.ico'
html_logo = '_static/images/logo-slepc.gif' # width should not exceed 200 pixels
# html_logo_light = '_static/images/logo-slepc.gif' # width should not exceed 200 pixels
# html_logo_dark =
# html_logo = html_logo_light

#html_title = 'SLEPc' #defaults to '{} {} documentation'.format(project,release)

html_css_files = [ # relative to the html_static_path
            'css/slepc.css',
            ]

# Files here are included in all the pages
#html_js_files = [ # relative to the html_static_path
#        ]

html_theme_options = {
        "icon_links": [ # icon links at the top (right)
            {
                "name": "GitLab",
                "url": "https://gitlab.com/slepc/slepc",
                "icon": "fa-brands fa-square-gitlab",
                "type": "fontawesome",
                },
            {
                "name": "UPV",
                "url": "https://www.upv.es",
                "icon": "https://www.upv.es/favicon.ico",
                "type": "local"
                },
            {
                "name": "Feed",
                "url": "_static/rss/slepc-news.xml",
                "icon": "fa-solid fa-square-rss",
                "type": "fontawesome",
                }
            ],
        "use_edit_page_button": True,
        "header_links_before_dropdown": 10, # before "more"
        "logo": {
            "alt_text": "SLEPc Home",
            "image_light": '_static/images/logo-slepc.gif',
            "image_dark": '_static/images/logo-slepc.gif'
            },
        "navigation_with_keys":True
        }

html_context = {
        "display_gitlab": True,
        "gitlab_user": "slepc",
        "gitlab_repo": "slepc",
        "gitlab_version": edit_branch,
        "doc_path": "doc/source",
        }

# Remove the primary (left) sidebar from all pages
#html_sidebars = {
#        "**": []
#        }

# -- Options for LaTeX output --------------------------------------------------
latex_engine = 'xelatex'

# How to arrange the documents into LaTeX files, building only the manual.
latex_documents = [
        ('documentation/manual/index', 'slepc-manual.tex', 'SLEPc Users Manual',
        'J. E. Roman, C. Campos, L. Dalcin, E. Romero, A. Tomas', 'manual', False)
        ]

latex_additional_files = [
     'documentation/manual/latex/frontpage.tex.txt',
     '_static/images/manual/pdf/logo-upv.pdf',
     '_static/images/manual/pdf/logo-dsic-black.pdf',
]

latex_show_pagerefs = True
latex_show_urls = 'footnote'

latex_elements = {
        'papersize' : 'a4paper',
        'pointsize' : '10pt',
        'extrapackages': r'\usepackage{xspace}',
        'sphinxsetup' : 'TableRowColorHeader={white},'
                        + 'TableRowColorOdd={gray}{0.97},'
                        + 'TableRowColorEven={white},'
                        + 'VerbatimColor={white},'
                        + 'noteBgColor={white},',
        'maketitle': r'\newcommand{\slepcversion}{%s}' % version
                     + r'\newcommand{\releasedate}{%s}' % release_date
                     + r'\makeatletter\@ifundefined{bibfont}{\newcommand{\bibfont}{\small}}{\renewcommand{\bibfont}{\small}}\makeatother'
r'''
\input{frontpage.tex.txt}
''',
        'tableofcontents' : r'',
        'printindex': r'''
        \printindex
        ''',
}

def setup(app):

#    print('-----------------------------8<--------------------------------')
#    print(app.project)
#    print('app.srcdir: {}'.format(app.srcdir))
#    print('app.confdir: {}'.format(app.confdir))
#    print('app.doctreedir: {}'.format(app.doctreedir))
#    print('app.outdir: {}'.format(app.outdir))
#    print('app.fresh_env_used: {}'.format(app.fresh_env_used))
#    print('-----------------------------8<--------------------------------')

    if 'PETSC_DIR' not in os.environ:
        print('\nUnable to build the documentation, PETSC_DIR environment variable is not set')
        print('\nPlease configure PETSc and SLEPc before building the documentation')
        raise Exception('PETSC_DIR not set')
    if 'PETSC_ARCH' not in os.environ:
        print('\nUnable to build the documentation, PETSC_ARCH environment variable is not set')
        print('\nPlease configure PETSc and SLEPc before building the documentation')
        raise Exception('PETSC_ARCH not set')
    else:
        # We should know where we are
        app.slepc_dir = os.path.abspath('../')
        app.petsc_dir = os.path.abspath(os.environ['PETSC_DIR'])
        app.petsc_arch = os.environ['PETSC_ARCH']

    doctext = shutil.which('doctext')
    if not doctext:
        with open(os.path.join(app.petsc_dir,app.petsc_arch,'lib','petsc','conf','petscvariables')) as f:
            doctext = [line for line in f if line.find('DOCTEXT ') > -1]
        if not doctext:
            print('PETSc has been configured without DOCTEXT. Not building manpages.')
        else:
            doctext = re.sub('[ ]*DOCTEXT[ ]*=[ ]*','',doctext[0]).strip('\n').strip()
            app.doctext = doctext
            print('Using DOCTEXT:', doctext)

    c2html = shutil.which('c2html')
    if not c2html:
        with open(os.path.join(app.petsc_dir,app.petsc_arch,'lib','petsc','conf','petscvariables')) as f:
            c2html = [line for line in f if line.find('C2HTML ') > -1]
        if not c2html:
            print('PETSc has been configured without C2HTML. Not building source pages.')
        else:
            c2html = re.sub('[ ]*C2HTML[ ]*=[ ]*','',c2html[0]).strip('\n').strip()
            app.c2html = c2html
            print('Using C2HTML:', c2html)

    mapnames = shutil.which('mapnames')
    if not mapnames:
        with open(os.path.join(app.petsc_dir,app.petsc_arch,'lib','petsc','conf','petscvariables')) as f:
            mapnames = [line for line in f if line.find('MAPNAMES ') > -1]
        if not mapnames:
            print('PETSc has been configured without MAPNAMES. Not building source pages.')
        else:
            mapnames = re.sub('[ ]*MAPNAMES[ ]*=[ ]*','',mapnames[0]).strip('\n').strip()
            app.mapnames = mapnames
            print('Using MAPNAMES:', mapnames)


    with open(os.path.join(app.slepc_dir, app.petsc_arch, 'include', 'slepcconf.h'),'r') as slepcconf_file:
        buf = slepcconf_file.read()
        slepc4py = re.search(' SLEPC_HAVE_SLEPC4PY[ ]*(1)',buf)
        if slepc4py is not None:
            app.slepc4py = True
            print('slepc4py: ',app.slepc4py)
        else:
            print('SLEPc configured without slepc4py: '
                  'slepc4py documentation will not be created')
            app.slepc4py = False

    if doctext and c2html:
        app.connect('builder-inited', builder_init_handler)
        app.connect('build-finished', build_finished_handler)

def builder_init_handler(app):
    global xtime
    if app.builder.name.endswith('html'):
        build_manpages_c2html.pre(app.doctext,app.slepc_dir,app.srcdir,app.outdir)
        _update_htmlmap_links(app)
        ptype = 'html'
    else: ptype = 'pdf'
    print("============================================")
    print("    Running Sphinx on SLEPc " + ptype)
    xtime = time.clock_gettime(time.CLOCK_REALTIME)

def build_finished_handler(app, exception):
    global xtime
    print("Time: "+str(time.clock_gettime(time.CLOCK_REALTIME) - xtime))
    print("============================================")
    if app.builder.name.endswith('html'):
        build_manpages_c2html.post(app.c2html,app.mapnames,app.slepc_dir,app.srcdir,app.outdir)
        if app.slepc4py:
            build_slepc4py_docs(app)
        _fix_links(app, exception)
        _fix_man_page_edit_links(app, exception)
#        if app.builder.name == 'dirhtml':
#            _add_man_page_redirects(app, exception)
        _add_man_page_redirects(app, exception)
        # remove sources for manual pages since they are automatically generated and should not be looked at on the website
        if os.path.isdir(os.path.join(app.outdir,'_sources','manualpages')):
            shutil.rmtree(os.path.join(app.outdir,'_sources','manualpages'))
        if app.builder.name == 'html':
            print('==========================================================================')
            print(f'    open {app.outdir}/index.html in your browser to view the documentation')
            print('==========================================================================')

def _add_man_page_redirects(app, exception):
    if exception is None:
        import time
        print("============================================")
        print("    Adding man pages redirects")
        x = time.clock_gettime(time.CLOCK_REALTIME)
        add_man_page_redirects.add_man_page_redirects(app.outdir)
        print("Time: "+str(time.clock_gettime(time.CLOCK_REALTIME) - x))
        print("============================================")

def build_slepc4py_docs(app):
    '''Builds the slepc4py docs and puts the results into the same directory tree as the SLEPc docs'''
    import time
    command = ['make', 'docsclean']
    print('============================================')
    print('Cleaning slepc4py docs')
    subprocess.run(command, cwd=os.path.join(app.slepc_dir,'src','binding','slepc4py'), check=True)

    command = ['make', 'website',
               'SLEPC_DIR={}'.format(app.slepc_dir),
               'PETSC_DIR={}'.format(app.petsc_dir),
               'PETSC_ARCH={}'.format(app.petsc_arch),
               'LOC={}'.format(app.outdir)]
    print('============================================')

    print('Building slepc4py docs')
    print(command)
    x = time.clock_gettime(time.CLOCK_REALTIME)
    subprocess.run(command, cwd=os.path.join(app.slepc_dir,'src','binding','slepc4py'), check=True)
    print("End slepc4py docs Time: "+str(time.clock_gettime(time.CLOCK_REALTIME) - x))
    print('============================================')

def _fix_man_page_edit_links(app, exception):
    if exception is None:
        import time
        print("============================================")
        print("    Fixing manual page edit links")
        x = time.clock_gettime(time.CLOCK_REALTIME)
        fix_man_page_edit_links.fix_man_page_edit_links(app.outdir)
        print("Time: "+str(time.clock_gettime(time.CLOCK_REALTIME) - x))
        print("============================================")

#   The following two scripts are needed because the Sphinx html and dirhtml builds save the output html
#   files at different levels of the directory hierarchy. file.rst/md -> file.html with html but
#   file.rst/md -> file/index.html with dirhtml and we want both to work correctly using relative links.

def _fix_links(app, exception):
    """We need to manage our own relative paths in the User's Manual for the source code files which
       are auto-generated by c2html outside of Sphinx so Sphinx cannot directly handle those links for use.
       We use the string PETSC_DOC_OUT_ROOT_PLACEHOLDER in URLs in the Sphinx .rst files as a stand in
       for the root directory that needs to be constructed based on if the Sphinx build is html or dirhtml
    """
    if exception is None:
        import time
        print("============================================")
        print("    Fixing relative links")
        x = time.clock_gettime(time.CLOCK_REALTIME)
        make_links_relative.make_links_relative(app.outdir)
        print("Time: "+str(time.clock_gettime(time.CLOCK_REALTIME) - x))
        print("============================================")

def _update_htmlmap_links(app):
    """htmlmap maps from manualpage names to relative locations in the generated documentation directory
       hierarchy. The format of the directory location needs to be different for the Sphinx html and dirhtml
       builds
    """
    import time
    print("============================================")
    print("    Updating htmlmap")
    x = time.clock_gettime(time.CLOCK_REALTIME)
    update_htmlmap_links.update_htmlmap_links(app.builder,os.path.join('source','manualpages','htmlmap'))
    print("Time: "+str(time.clock_gettime(time.CLOCK_REALTIME) - x))
    print("============================================")
