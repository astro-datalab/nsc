#!/usr/bin/env python

from distutils.core import setup

setup(name='nsc',
      version='1.0.0',
      description='NOIRLab Source Catalog Processing Software',
      author='David Nidever',
      author_email='dnidever@montana.edu',
      url='https://github.com/astro-datalab/nsc',
      packages=['nsc'],
      package_dir={'':'python'},
      package_data={'nsc': ['data/*','data/params/*','data/params/*/*']},
      scripts=['bin/nscjobs','bin/nsc_instcal_measure','bin/nsc_instcal_calibrate',
               'bin/nsc_instcal_calibrate_healpix','bin/nsc_instcal_combine',
               'bin/nsc_archive_search','bin/decam_archive_search','bin/decam_parse_archive_search'],
      #py_modules=['nsc_instcal',''],
      requires=['numpy','astropy','scipy','dlnpyutils','sep','healpy','dustmaps','astroquery'],
      include_package_data=True
)
