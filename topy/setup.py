#!/usr/bin/env python3

from distutils.core import setup

setup(name='ToPy',
      version='0.5',
      description='Topology optimization with Python',
      author='William Hunter. Ported to python3 by Francisco Sanchez',
      url='https://github.com/thebeachlab/topy',
      packages=['topy', 'topy.data'],
      package_dir={'topy': ''}
      )
