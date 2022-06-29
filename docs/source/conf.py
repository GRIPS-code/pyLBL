# Sphinx documentation builder configuration file.
from os.path import abspath
import sys

sys.path.insert(0, abspath(".."))

# Project information.
project = "pyLBL"
copyright = "2021, GRIPS-code"
author = "GRIPS-code team"
release = "alpha"

# Sphinx extension module names.
extensions = ["sphinx.ext.autodoc",
              "sphinx.ext.autosummary",
              "sphinx.ext.napoleon",
              "sphinx.ext.viewcode",
              "sphinxcontrib.apidoc",
              "sphinx_autopackagesummary"]

# Templates.
templates = ["_templates"]

# Patterns, files, and directories to ignore.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# Apidoc settings
apidoc_module_dir = "../../pyLBL"

# Napoleon settings.
napoleon_google_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_numpy_docstring = False
napolean_use_keyword = True

# Autosummary settings.
autosummary_generate = True

# Theme for the HTML pages.
html_theme = "bizstyle"

# Paths to static files (like style sheets).
html_static_path = ["_static"]
