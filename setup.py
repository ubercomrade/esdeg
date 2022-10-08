#!/usr/bin/env python3
from setuptools import setup

try:
    from pythran.dist import PythranExtension, PythranBuildExt
    setup_args = {
        'cmdclass': {"build_ext": PythranBuildExt},
        'ext_modules': [PythranExtension('esdeg.speedup', sources = ['esdeg/speedup.py'],  extra_compile_args=["-Ofast", "-DUSE_XSIMD"])],
    }
except ImportError:
    print("Not building Pythran extension")
    setup_args = {}
    
with open("README.md", "r") as fh:
    long_description = fh.read()

install_requires = ['setuptools>=18.0', 'pandas', 'numpy', 'scipy', 'statsmodels', 'pythran', 'biopython', 'pyjaspar']
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
        'data': ['*.fa'],
    },
    scripts=['bin/ESDEG.py',],
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
    python_requires='>=3.7',
    **setup_args)
