# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import sys
from pathlib import Path

# Add source to path for autodoc2
sys.path.insert(0, str(Path("../src").resolve()))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "airsspy"
copyright = "2024, Bonan Zhu"
author = "Bonan Zhu"
release = "0.1.4"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "myst_nb",                      # MyST Markdown + Jupyter notebooks
    "autodoc2",                     # Automatic API documentation
    "sphinx.ext.intersphinx",       # Link to other docs
    "sphinx.ext.viewcode",          # Add source code links
    "sphinx_design",                # Grid layouts, cards
    "sphinx_copybutton",            # Copy code button
    "sphinx_togglebutton",          # Collapsible sections
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "**.ipynb_checkpoints"]

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "pydata_sphinx_theme"
html_static_path = ["_static"]
html_css_files = ["custom.css"]

html_theme_options = {
    "github_url": "https://github.com/zhubonan/airsspy",
    "logo": {
        "text": "airsspy",
    },
    "navbar_end": ["navbar-icon-links", "theme-switcher"],
    "footer_start": ["copyright"],
    "footer_end": ["sphinx-version"],
}

# -- MyST configuration ------------------------------------------------------
# https://myst-parser.readthedocs.io/en/latest/configuration.html

myst_enable_extensions = [
    "colon_fence",      # ::: fences
    "deflist",          # Definition lists
    "dollarmath",       # $math$ syntax
    "fieldlist",        # Field lists
    "substitution",     # Variable substitution
]

# Execute notebooks
nb_execution_mode = "cache"
nb_execution_timeout = 300

# -- autodoc2 configuration --------------------------------------------------
# https://sphinx-autodoc2.readthedocs.io/en/latest/

autodoc2_packages = [
    {
        "path": "../src/airsspy",
        "exclude_files": ["_version.py"],
    }
]
autodoc2_render_plugin = "myst"
autodoc2_hidden_objects = ["private", "inherited"]

# -- Intersphinx configuration -----------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/extensions/intersphinx.html

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "ase": ("https://wiki.fysik.dtu.dk/ase/", None),
}

# -- Copybutton configuration ------------------------------------------------
# https://sphinx-copybutton.readthedocs.io/en/latest/

copybutton_prompt_text = r">>> |\.\.\. |\$ "
copybutton_prompt_is_regexp = True
