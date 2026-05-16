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
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

suppress_warnings = ["myst.header"]

napoleon_google_docstring = False
napoleon_numpy_docstring = True

autodoc_member_order = "bysource"

html_theme = "alabaster"

os.chdir(ROOT_DIR)
