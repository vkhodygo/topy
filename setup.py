#!/usr/bin/env python
"""
ToPy install script.

Install ToPy through `pip install .`.
"""
import setuptools

# Get project description.
with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    version="0.4.0",
    name="topology-optimization",
    author="William Hunter",
    author_email="williamhunter@users.noreply.github.com",
    description="Topology optimization with Python",
    url="https://github.com/williamhunter/topy",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=["topy", "topy.data"],
    install_requires=[
        "matplotlib",
        "sympy",
        "numpy",
        "scipy",
        "scikit-sparse",
        "pyvtk",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Development Status :: 3 - Alpha",
    ],
    python_requires=">=3.8",
    keywords="topology optimization vtk",
)
