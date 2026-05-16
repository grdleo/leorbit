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

if ! "$PYTHON_BIN" -c "import sphinx, myst_parser, nbconvert" >/dev/null 2>&1; then
  echo "Documentation dependencies are missing in this environment."
  echo "Installing Sphinx + MyST parser + nbconvert..."
  "$PYTHON_BIN" -m pip install "sphinx>=7" "myst-parser>=3" "nbconvert>=7"
fi

DOCS_SOURCE_DIR="$ROOT_DIR/docs/source"
DOCS_BUILD_DIR="$ROOT_DIR/docs/build"
DOCS_API_DIR="$DOCS_SOURCE_DIR/api"
DOCS_EXAMPLES_DIR="$DOCS_SOURCE_DIR/examples"
NOTEBOOKS_DIR="$ROOT_DIR/examples"

mkdir -p "$DOCS_SOURCE_DIR"
mkdir -p "$DOCS_API_DIR"
mkdir -p "$DOCS_EXAMPLES_DIR"

find "$DOCS_EXAMPLES_DIR" -maxdepth 1 -type f -name '*.md' ! -name 'index.md' -delete

for notebook in "$NOTEBOOKS_DIR"/*.ipynb; do
  [ -e "$notebook" ] || continue
  notebook_basename="$(basename "$notebook" .ipynb)"
  output_md="$DOCS_EXAMPLES_DIR/$notebook_basename.md"

  echo "Converting notebook: $(basename "$notebook")"
  "$PYTHON_BIN" -m jupyter nbconvert \
    --to markdown \
    --output-dir "$DOCS_EXAMPLES_DIR" \
    "$notebook"

  # Sphinx toctree requires each page to expose a title.
  if [[ -f "$output_md" ]]; then
    tmp_file="${output_md}.tmp"
    title="${notebook_basename//_/ }"
    {
      printf "# %s\n\n" "$title"
      cat "$output_md"
    } > "$tmp_file"
    mv "$tmp_file" "$output_md"
  fi
done

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
