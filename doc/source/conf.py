# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

from os import path
import re
import subprocess
from datetime import datetime


# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'SLEPc'
copyright = '2002-%Y, Universitat Politecnica de Valencia, Spain'
author = 'SLEPc Team'
version = 'development'
release = 'development'

with open(path.join('../..', 'include', 'slepcversion.h'),'r') as version_file:
    buf = version_file.read()
    petsc_release_flag = re.search(' SLEPC_VERSION_RELEASE[ ]*([0-9]*)',buf).group(1)
    major_version      = re.search(' SLEPC_VERSION_MAJOR[ ]*([0-9]*)',buf).group(1)
    minor_version      = re.search(' SLEPC_VERSION_MINOR[ ]*([0-9]*)',buf).group(1)
    subminor_version   = re.search(' SLEPC_VERSION_SUBMINOR[ ]*([0-9]*)',buf).group(1)

    git_describe_version = subprocess.check_output(['git', 'describe', '--always']).strip().decode('utf-8')
    if petsc_release_flag == '0':
        version = git_describe_version
        release = git_describe_version
    else:
        version = '.'.join([major_version, minor_version])
        release = '.'.join([major_version,minor_version,subminor_version])

release_date = subprocess.check_output(['git',
                                 'for-each-ref',
                                 '--format="%(creatordate:format:%B), %(creatordate:format:%Y)"',
                                 'refs/tags/v{}.0'.format(version)]
                                ).strip().decode('utf-8').strip('"')

if release_date == '':
    release_date = datetime.strftime(datetime.now(), '%B, %Y')

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

needs_sphinx='8.1'

#nitpicky = True # warn about all references where the target cannot be found

extensions = [
        'myst_parser',         # use structured markdown
        'sphinx.ext.duration',
        'sphinx_copybutton',   # Add a copy button to your code blocks
        'sphinx_design',       # tabs
        'sphinx_togglebutton', # dropdown box
        'sphinxcontrib.bibtex',
        ]

myst_links_external_new_tab = True # open all external links in new tabs

bibtex_bibfiles = ['slepc.bib']
##########################################################################
import pybtex.plugin
from pybtex.style.formatting.alpha import Style as AlphaStyle
from pybtex.style.labels.alpha import LabelStyle

class SLEPcLabelStyle(LabelStyle):
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
#        'attrs_inline', # parsing of inline attributes
        'dollarmath',   # use of $$ latex math
#        'html_image',
        'fieldlist',
        ]

myst_dmath_double_inline = True
myst_dmath_allow_labels = True # the default
# XXX check these two
myst_dmath_allow_space = True
myst_dmath_allow_digits=False

myst_substitutions = {
            'release_date': release_date
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
exclude_patterns = []



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

html_title = 'SLEPc' #defaults to '{} {} documentation'.format(project,release)

html_css_files = [ # relative to the html_static_path
            'css/slepc.css',
            ]

#html_js_files = [ # relative to the html_static_path
# XXX these are included in all the pages!
#        'js/slepc.js',
#        'js/news.js',
#        'js/apps.js',
#        ]


html_theme_options = {
#        "gitlab_url": "https://gitlab.com/slepc/slepc2.old",
#        "external_links": [
#            {"name": "link-one-name", "url": "https://<link-one>"},
#            {"name": "link-two-name", "url":
#             "https://<link-two>"}
#            ],

        "icon_links": [ # icon links at the top (right)
            {
                "name": "GitLab",
                "url": "https://gitlab.com/slepc/slepc",
                #"icon": "fa-brands fa-gitlab",
                "icon": "fa-brands fa-square-gitlab",
                "type": "fontawesome",
                },
            {
                "name": "UPV",
                "url": "https://www.upv.es",
                "icon": "_static/images/new/favicon-upv.ico",
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
#        "footer_end": ["theme-version", "last-updated"],
#        "secondary_sidebar_items" : ["edit-this-page"],
        "header_links_before_dropdown": 10, # before "more"
        "logo": {
#            "text": 'SLEPc-ng',
            "alt_text": "SLEPc Home",
            "image_light": '_static/images/logo-slepc.gif',
            "image_dark": '_static/images/logo-slepc.gif'
            },
        "navigation_with_keys":True
        }

try:
    git_ref = subprocess.check_output(["git", "rev-parse", "HEAD"]).rstrip()
    git_ref_release = subprocess.check_output(["git", "rev-parse", "origin/release"]).rstrip()
    edit_branch = "release" if git_ref == git_ref_release else "main"
except subprocess.CalledProcessError:
    print("WARNING: determining branch for page edit links failed")
    edit_branch = "main"

html_context = {
        "display_gitlab": True,
        "gitlab_user": "slepc",
        "gitlab_repo": "slepc",
        "gitlab_version": edit_branch,
        "doc_path": "doc",
        }

# remove the primary (left) sidebar from all pages
#html_sidebars = {
#        "**": []
#        }
