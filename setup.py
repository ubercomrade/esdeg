#!/usr/bin/env python3
from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

install_requires = ['setuptools>=18.0', 'pandas', 'numpy', 'scipy>=1.10.0', 'plotly', 'biopython', 'pyjaspar>=3.0', 'xlsxwriter', 'panel']
setup(
    name='ESDEG',
    version='0.0.1',
    description='DEG analisys',
    author='Anton Tsukanov',
    author_email='tsukanov@bionet.nsc.ru',
    url='http://github.com/ubercomrade/enrest',
    package_dir={'esdeg' : 'esdeg'},
    packages=['esdeg'],
    package_data={
        'esdeg': ['tomtom/*.tsv'],
        'esdeg': ['clusters/*.tsv'],
        'esdeg': ['logos/*.png'],
    },
    #scripts=['bin/ESDEG.py',],
    entry_points={
        "console_scripts": [
            "ESDEG = bin.ESDEG:main",
        ]
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Unix",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
    zip_safe=False,
    include_package_data=True,
    install_requires=install_requires,
    setup_requires=install_requires,
    python_requires='>=3.7')
