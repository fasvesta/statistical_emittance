# copyright ############################################## #
# This file is part of the statisticalEmittance Package.   #
# Copyright (c) CERN, 2022.                                #
# ######################################################## #

from setuptools import setup, find_packages, Extension

#######################################
# Prepare list of compiled extensions #
#######################################

extensions = []


#########
# Setup #
#########

setup(
    name='statisticalEmittance',
    version='0.1.0',
    description='calculation of statistical emittance',
    long_description=('Class to calculate the statistical emittance od an ensemble of particles'
                'Optics functions such as beta & dispersion functions are also statistically estimated.'),
    url='https://github.com/fasvesta/statistical_emittance',
    packages=find_packages(),
    ext_modules = extensions,
    include_package_data=True,
    install_requires=[
        'numpy>=1.0'
        ],
    author='F. Asvesta et al.',
    license='Apache 2.0',
    project_urls={
            "Source Code": "https://github.com/fasvesta/statistical_emittance",
        },
    )
