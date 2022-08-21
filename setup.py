#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ToPy install script.

Install ToPy through `pip install .`.
"""

import json

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

        "typing",
        "pathlib",
        "matplotlib",
        "sympy",
        "numpy<=1.14",
        "scipy",

        "pyvtk",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Development Status :: 3 - Alpha",
    ],
    python_requires=">=3",
    **metadata
)
