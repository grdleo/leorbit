from __future__ import annotations

import os
import sys
from pathlib import Path


ROOT_DIR = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT_DIR))

project = "LEOrbit"
author = "Leo Giroud"
release = "0a1"

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.mathjax",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

myst_enable_extensions = [
    "dollarmath",
    "amsmath",
]

suppress_warnings = ["myst.header"]

napoleon_google_docstring = False
napoleon_numpy_docstring = True

autodoc_member_order = "bysource"

html_theme = "alabaster"

os.chdir(ROOT_DIR)
