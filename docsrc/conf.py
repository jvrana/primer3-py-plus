#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Trident documentation build configuration file, created by
# sphinx-quickstart on Sun Nov 19 22:18:51 2017.
#
# This file is execfile()d with the current directory set to its
# containing dir.
#
# Note that not all possible configuration values are present in this
# autogenerated file.
#
# All configuration values have a default; values that are commented out
# serve to show the default.
# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
import os
import sys

sys.path.insert(0, os.path.abspath(".."))

import primer3plus


def boulderio_to_pandas(boulderio):
    import pandas as pd

    rows = []
    for k, v in boulderio._params.items():
        name = "`{} <{}>`_".format(v.name, v.help())
        if v.ptype.category == "extra":
            name = ":attr:`{x} <primer3plus.params.ExtraTypes.{x}>`".format(x=v.name)
        rows.append({"name": name, "type": str(v.ptype.type), "value": v.value})
    return pd.DataFrame(rows)


def boulderio_to_csv(boulderio):
    df = boulderio_to_pandas(boulderio)
    df.to_csv("./_static/boulderio.csv", index=False)


boulderio_to_csv(primer3plus.Design.DEFAULT_PARAMS)

# -- General configuration ------------------------------------------------

# AUTODOC
autoclass_content = "both"  # include both class docstring and __init__
autodoc_default_flags = [
    # Make sure that any autodoc declarations show the right members
    "members",
    "inherited-members",
    "private-members",
    "show-inheritance",
]
autosummary_generate = True  # Make _autosummary files and include them
napoleon_numpy_docstring = False  # Force consistency, leave only Google
napoleon_use_rtype = False  # More legible


# If your documentation needs a minimal Sphinx version, state it here.
#
# needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "autodocsumm",
    "sphinx_autodoc_typehints",
    "sphinx.ext.mathjax",
    "sphinx.ext.doctest",
    "sphinx.ext.coverage",
    "sphinx.ext.viewcode",
    "sphinx.ext.inheritance_diagram",
    "recommonmark",
]
autodoc_default_options = {"autosummary": True}

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
source_suffix = {".rst": "restructuredtext", ".txt": "markdown", ".md": "markdown"}
# source_suffix = '.rst'

# The master toctree document.
master_doc = "index"

# General information about the project.
project = primer3plus.__title__
copyright = "2017-2019, University of Washington"
author = ", ".join(primer3plus.__authors__)

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
# The short X.Y version.
version = primer3plus.__version__
homepage = primer3plus.__homepage__
# The full version, including alpha/beta/rc tags.
release = primer3plus.__version__

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = None

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This patterns also effect to html_static_path and html_extra_path
exclude_patterns = ["docs", "Thumbs.db", ".DS_Store", "_build", "**.ipynb_checkpoints"]

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = "sphinx"

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = False


html_context = {
    "version": version,
    "display_github": True,  # Integrate GitHub
    "github_user": "jvrana",  # Username
    "github_repo": "primer3plus-dna-design",  # Repo name
    "github_version": "master",  # Version
    "conf_py_path": "./",  # Path in the checkout to the docs root
}

# substitutations for the docsrc
substitutions = {"homepage": primer3plus.__homepage__, "repo": primer3plus.__repo__}

rst_epilog = "\n".join(
    ".. |{k}| replace:: {v}".format(k=k, v=v) for k, v in substitutions.items()
)

# -- Options for HTML output ----------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#

html_theme = "sphinx_rtd_theme"
html_theme_path = ["_themes"]

# Guzzle theme options (see theme.conf for more information)
html_theme_options = {
    "sticky_navigation": True,
    "display_version": True,
    "navigation_depth": 4,
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".

# html_sidebars = {
#     'index':    ['sidebar.html', 'globaltoc.html', 'relations.html', 'sourcelink.html', 'searchbox.html'],
#     '**':       ['sidebar.html', 'localtoc.html', 'relations.html', 'sourcelink.html', 'searchbox.html']
# }
# html_sidebars = {
#     'index':    ['sidebar.html'],
#     '**':       ['sidebar.html']
# }

# -- Options for HTMLHelp output ------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = "primer3plusdoc"


# -- Options for LaTeX output ---------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    # 'papersize': 'letterpaper',
    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',
    # Additional stuff for the LaTeX preamble.
    #
    # 'preamble': '',
    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (
        master_doc,
        "primer3plus.tex",
        "primer3plus Documentation",
        "Justin Vrana, Eric Klavins, Ben Keller",
        "manual",
    )
]


# -- Options for manual page output ---------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [(master_doc, "primer3plus", "primer3plus Documentation", [author], 1)]


# -- Options for Texinfo output -------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (
        master_doc,
        "primer3plus",
        "primer3plus Documentation",
        author,
        "primer3plus",
        "One line description of project.",
        "Miscellaneous",
    )
]

# def setup(app):
#     app.add_stylesheet('css/custom.css')
#

# Default language for syntax highlighting in reST and Markdown cells
highlight_language = "none"

# Don't add .txt suffix to source files (available for Sphinx >= 1.5):
html_sourcelink_suffix = ""

# Work-around until https://github.com/sphinx-doc/sphinx/issues/4229 is solved:
html_scaled_image_link = False

# List of arguments to be passed to the kernel that executes the notebooks:
nbsphinx_execute_arguments = [
    "--InlineBackend.figure_formats={'svg', 'pdf'}",
    "--InlineBackend.rc={'figure.dpi': 96}",
]

# This is processed by Jinja2 and inserted before each notebook
nbsphinx_prolog = r"""
{% set docname = env.doc2path(env.docname, base='doc') %}
.. only:: html
    .. role:: raw-html(raw)
        :format: html
    .. nbinfo::
        This page was generated from `{{ docname }}`__.
        Interactive online version:
        :raw-html:`<a href="https://mybinder.org/v2/gh/spatialaudio/nbsphinx/{{ env.config.release }}?filepath={{ docname }}"><img alt="Binder badge" src="https://mybinder.org/badge_logo.svg" style="vertical-align:text-bottom"></a>`
    __ https://github.com/spatialaudio/nbsphinx/blob/
        {{ env.config.release }}/{{ docname }}
.. raw:: latex
    \nbsphinxstartnotebook{\scriptsize\noindent\strut
    \textcolor{gray}{The following section was generated from
    \sphinxcode{\sphinxupquote{\strut {{ docname | escape_latex }}}} \dotfill}}
"""

# This is processed by Jinja2 and inserted after each notebook
nbsphinx_epilog = r"""
.. raw:: latex
    \nbsphinxstopnotebook{\scriptsize\noindent\strut
    \textcolor{gray}{\dotfill\ \sphinxcode{\sphinxupquote{\strut
    {{ env.doc2path(env.docname, base='doc') | escape_latex }}}} ends here.}}
"""

mathjax_path = (
    "https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"
)

mathjax_config = {
    "TeX": {"equationNumbers": {"autoNumber": "AMS", "useLabelIds": True}}
}
