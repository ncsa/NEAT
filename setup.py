#!/usr/bin/env python3


__author__ = 'Zach Stephens'
__copyright__ = 'Copyright (c) 1999-2021'
__license__ = 'BSD 3-Clause License'
__version__ = '3.0.1'
__email__ = 'zstephe2@illinois.edu'
__status__ = 'prod'

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

with open('requirements.txt') as f:
    required = f.read().splitlines()

setuptools.setup(
    name="NEAT",
    version=__version__,
    author=__author__,
    author_email=__email__,
    description="Variant and read simulator for benchmarking NGS workflows",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ncsa/neat-genreads",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Development Status :: 5 - Production/Stable",
        "Programming Language :: Python :: Implementation :: PyPy",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    install_requires=required,
    scripts=[
        'gen_reads.py',
        'bacterial_genreads_wrapper.py',
        'utilities/compute_fraglen.py',
        'utilities/compute_gc.py',
        'utilities/genSeqErrorModel.py',
        'utilities/gen_mut_model.py',
        'utilities/generate_random_dna.py',
        'utilities/plotMutModel.py',
        'utilities/validateBam.py',
        'utilities/validateFQ.py',
    ],
    project_urls={  # Optional
        'Bug Reports': 'https://github.com/ncsa/NEAT/issues',
        'Source': 'https://github.com/ncsa/NEAT',
    }
)
