# coding=utf-8

"""
    vcfy.release
    ~~~~~~~~~~~~

    Module containing package release information.

    :copyright: (c) 2018 by Ali Ghaffaari.
    :license: MIT, see LICENSE for more details.
"""

# CONSTANTS ###################################################################
# Development statuses:
DS_PLANNING = 1
DS_PREALPHA = 2
DS_ALPHA = 3
DS_BETA = 4
DS_STABLE = 5
DS_MATURE = 6
DS_INACTIVE = 7
DS_STRING = {
    DS_PLANNING: 'Development Status :: 1 - Planning',
    DS_PREALPHA: 'Development Status :: 2 - Pre-Alpha',
    DS_ALPHA: 'Development Status :: 3 - Alpha',
    DS_BETA: 'Development Status :: 4 - Beta',
    DS_STABLE: 'Development Status :: 5 - Production/Stable',
    DS_MATURE: 'Development Status :: 6 - Mature',
    DS_INACTIVE: 'Development Status :: 7 - Inactive'
}
###############################################################################

# Package release information.
__title__ = 'vcfy'
__description__ = 'Generate VCF file with random variants from reference genome'
__author__ = 'Ali Ghaffaari'
__email__ = 'ali.ghaffaari@mpi-inf.mpg.de'
__license__ = 'MIT'

# Release
__version__ = '0.0.8'
__status__ = DS_PREALPHA

# Package data
__package_data__ = {__title__: ['resources/templates/template.vcf']}

# PyPI-related information
__keywords__ = 'vcf genome genomics variation DNA sequence ngs'
__classifiers__ = [
    # Development status
    DS_STRING[__status__],

    # Environment
    'Environment :: Console',

    # License
    'License :: OSI Approved :: MIT License',

    # Supported Python versions.
    'Programming Language :: Python :: 2',
    'Programming Language :: Python :: 2.7',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',

    # Intended Audience and Topic
    'Intended Audience :: Science/Research',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
]
__requires__ = [
    'PyVCF>=0.6.8',
    'biopython>=1.71',
    'numpy>=1.14.3',
    'BitVector>=3.4.8',
    'click>=6.7',
]
__tests_require__ = []
__extras_require__ = {
    'test': ['nose>=1.0', 'coverage'],
}
__setup_requires__ = ['nose>=1.0', 'coverage']
__entry_points__ = '''
[console_scripts]
vcfy=vcfy.cli:cli
ksnper=vcfy.ksnper:cli
'''
