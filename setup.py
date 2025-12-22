#!/usr/bin/env python3
"""
longRNA - RNA-seq Analysis Pipeline (Salmon & RSEM)
Setup script for installation
"""

import setuptools
from pathlib import Path

# Read README for long description
readme_file = Path(__file__).parent / "longRNA" / "README.md"
if readme_file.exists():
    with open(readme_file, 'r') as f:
        long_description = f.read()
else:
    long_description = "RNA-seq Analysis Pipeline with Salmon and RSEM workflows"

# Package metadata
NAME = 'longRNA'
DESCRIPTION = 'RNA-seq Analysis Pipeline - Salmon and RSEM workflows'
AUTHOR = 'shliu'
REQUIRES_PYTHON = '>=3.7.0'
VERSION = '1.0.0'

# Required packages
REQUIRED = [
    'pyyaml>=5.1',
]

setuptools.setup(
    name=NAME,
    version=VERSION,
    author=AUTHOR,
    description=DESCRIPTION,
    long_description=long_description,
    long_description_content_type='text/markdown',
    packages=setuptools.find_packages(),
    include_package_data=True,
    install_requires=REQUIRED,
    python_requires=REQUIRES_PYTHON,
    entry_points={
        'console_scripts': [
            'longRNA = longRNA.entry:main',
        ],
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'Operating System :: Unix',
    ],
    zip_safe=False,
)
