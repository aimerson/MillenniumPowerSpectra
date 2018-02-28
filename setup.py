#! /usr/bin/env python

from setuptools import setup, find_packages

datafiles = ['datasets/millenniumWMAP1_matterPowerSpectra.hdf5']

setup(name='millenniumpk',
      version='0.1',
      description='Scripts for extracting power spectra for the Millennium Simulation and computing corresponding two-point correlation function.',
      url='http://bitbucket.org/aimerson/millenniumpk',
      author='Alex Merson',
      author_email='alex.i.merson@gmail.com',
      license='MIT',
      packages=find_packages(),
      package_dir={'millenniumpk':'millenniumpk'},
      package_data={'millenniumpk':datafiles},
      zip_safe=False)

