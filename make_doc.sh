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

if ! "$PYTHON_BIN" -c "import pdoc" >/dev/null 2>&1; then
  echo "pdoc is not installed in this environment."
  echo "Installing pdoc..."
  "$PYTHON_BIN" -m pip install "pdoc>=14"
fi

DOCS_OUT_DIR="$ROOT_DIR/docs/api"
rm -rf "$DOCS_OUT_DIR"
mkdir -p "$DOCS_OUT_DIR"

"$PYTHON_BIN" -m pdoc \
  --docformat numpy \
  --output-directory "$DOCS_OUT_DIR" \
  leorbit

echo "Documentation generated in: $DOCS_OUT_DIR"
echo "Open: $DOCS_OUT_DIR/leorbit.html"
