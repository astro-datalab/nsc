#!/usr/bin/env python

from distutils.core import setup

setup(name='noaosourcecatalog',
      version='1.0.0',
      description='NOIRLab Source Catalog Processing Software',
      author='David Nidever',
      author_email='dnidever@montana.edu',
      url='https://github.com/astro-datalab/noaosourcecatalog',
      packages=['nsc'],
      package_dir={'':'python'},
      scripts=['bin/nsc_instcal_calibrate','bin/nsc_instcal_combine'],
      #py_modules=['nsc_instcal',''],
      requires=['numpy','astropy','scipy','dlnpyutils','sep','healpy']
)
