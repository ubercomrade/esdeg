"""kek"""
from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

install_requires = ['setuptools>=18.0', 'numba', 'numpy', 'scipy', 'matplotlib']


setup(
    name='enrest',
    version='0.0.1',
    description='DEG analisys',
    author='Anton Tsukanov',
    author_email='tsukanov@bionet.nsc.ru',
    url='http://github.com/ubercomrade/enrest',
    package_dir={'enrest' : 'enrest'},
    packages=[
        'enrest',
    ],
    package_data={
        'data': ['*.fasta'],
    },
    scripts=['bin/enrest.py',],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Unix",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
    zip_safe=False,
    include_package_data=True,
    install_requires=install_requires,
    setup_requires=install_requires,
    python_requires='>=3.7',
)