import nbformat
from nbformat import NotebookNode
from nbconvert import RSTExporter
from glob import glob
from pathlib import Path
from typing import TypedDict

EXAMPLES_IPYNB_SOURCE_SCHEMA = r"examples/*.ipynb"
EXAMPLES_RST_TARGET_DIR = r"docs/source/examples"
EXAMPLE_RST_NAME = r"index.rst"

RST = RSTExporter()

class NotebookConversion(TypedDict):
    notebook: NotebookNode
    rst: str
    name: str

def convert_notebook(path: Path) -> str:
    """Returns: name"""
    name = path.stem
    nb: NotebookNode

    with open(path, "r") as f:
        nb = nbformat.reads(f.read(), as_version=4)
    body, res = RST.from_notebook_node(nb)

    outputs = res["outputs"]
    if outputs:
        raise NotImplementedError("Please handle outputs!")
    
    output_dir = Path(EXAMPLES_RST_TARGET_DIR) / name / EXAMPLE_RST_NAME
    output_dir.parent.mkdir(parents=True, exist_ok=True)

    with open(output_dir, "a") as f:
        f.write(body)
    
    return name

if __name__ == "__main__":
    nb_names = [convert_notebook(Path(p)) for p in glob(EXAMPLES_IPYNB_SOURCE_SCHEMA)]
    
    with open(Path(EXAMPLES_RST_TARGET_DIR) / "index.rst", "a") as f:
        f.writelines(
            f"   {name}/index" for name in nb_names
        )