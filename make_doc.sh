#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT_DIR"

if [[ -x "$ROOT_DIR/.venv/bin/python" ]]; then
  PYTHON_BIN="$ROOT_DIR/.venv/bin/python"
else
  PYTHON_BIN="python3"
fi

echo "Using Python: $PYTHON_BIN"

if ! "$PYTHON_BIN" -c "import sphinx, myst_parser" >/dev/null 2>&1; then
  echo "Sphinx dependencies are missing in this environment."
  echo "Installing Sphinx + MyST parser..."
  "$PYTHON_BIN" -m pip install "sphinx>=7" "myst-parser>=3"
fi

DOCS_SOURCE_DIR="$ROOT_DIR/docs/source"
DOCS_BUILD_DIR="$ROOT_DIR/docs/build"
DOCS_API_DIR="$DOCS_SOURCE_DIR/api"

mkdir -p "$DOCS_API_DIR"

"$PYTHON_BIN" -m sphinx.ext.apidoc \
  --force \
  --output-dir "$DOCS_API_DIR" \
  "$ROOT_DIR/leorbit" \
  "$ROOT_DIR/leorbit/tests"

"$PYTHON_BIN" -m sphinx \
  -M html \
  "$DOCS_SOURCE_DIR" \
  "$DOCS_BUILD_DIR"

echo "Documentation generated in: $DOCS_BUILD_DIR/html"
echo "Open: $DOCS_BUILD_DIR/html/index.html"
