import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="leorbit",
    licence="MIT",
    version="0.0.0",
    author="Léo Giroud",
    author_email="leo@leog.dev",
    description="Python library for satellites in LEO (Low Earth Orbit): propagation, predictions & more",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/grdleo",
    packages=setuptools.find_packages(exclude=["tests"]),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Development Status :: 3 - Alpha",
    ],
    install_requires=[
        "pint",
        "requests",
        "beautifulsoup4",
        "pdoc"
    ],
    python_requires=">=3.12"
)

