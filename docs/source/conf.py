# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import re
from datetime import datetime

# import sys

# sys.path.insert(0, os.path.abspath('.'))

author_setup_str = "    author="
author_email_setup_str = "    author_email="
name_setup_str = "    name="
version_setup_str = "    version="


# Load the release info into a dict by explicit execution
info = {}
with open(os.path.abspath(os.path.join(
        os.path.dirname(__file__), "../..", "setup.py")), "r") as f:
    for line in f:
        if line.startswith(name_setup_str):
            project = (
                re.search(name_setup_str + "(.*),", line).group(1).strip("\'"))
        elif line.startswith(author_setup_str):
            _author = (
                re.search(
                    author_setup_str + "(.*),",
                    line).group(1).strip("\'").replace("\\", "")
            )
        elif line.startswith(author_email_setup_str):
            _email = (
                re.search(
                    author_email_setup_str + "(.*),",
                    line).group(1).strip("\'")
            )
        elif line.startswith(version_setup_str):
            _version = (
                re.search(
                    version_setup_str + "(.*),",
                    line).group(1).strip("\'")
            )

# -- Project information -----------------------------------------------------

copyright = f"2013-{datetime.now().year}, {_author} <{_email}>"
author = f"{_author}s"

# The full version, including alpha/beta/rc tags
release = _version


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
]

# Add any paths that contain templates here, relative to this directory.
templates_path = []

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["*tests*"]

# Sources
source_suffix = [".rst", ".md"]

# The master toctree document.
master_doc = "index"


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "furo"

html_favicon = "../_static/icon/favicon.ico"
html_logo = "../_static/wma_small_alpha.png"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = []
