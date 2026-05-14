import setuptools
from pathlib import Path

setuptools.setup(
    name="leorbit",
    license="MIT",
    author="Léo Giroud",
    author_email="leo@leog.dev",
    url="https://github.com/grdleo/leorbit",

    description="Python library for satellites in LEO (Low Earth Orbit): propagation, predictions & more",\
    long_description=Path("README.md").read_text(),
    long_description_content_type="text/markdown",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Development Status :: 3 - Alpha",
    ],
    
    version="0a1", 
    # Simple versioning scheme, with one component is the strategy. 
    # https://packaging.python.org/en/latest/discussions/versioning/
    
    packages=setuptools.find_packages(exclude=["tests"]),
    install_requires=[
        "numpy",
        "matplotlib",
        "pint",
        "requests",
        "pydantic"
    ],
    extras_require={
        "docs": [
            "pdoc>=14",
        ],
    },
    python_requires=">=3.14"
)

