#!/usr/bin/env python
# -*- encoding: utf-8 -*-
import io
import re
from glob import glob
from os.path import basename, dirname, join, splitext

from setuptools import find_packages
from setuptools import setup


def read(*names, **kwargs):
    with io.open(
        join(dirname(__file__), *names),
        encoding=kwargs.get('encoding', 'utf8')
    ) as fh:
        return fh.read()


setup(
    name='analyze_molcas',
    version='0.0.0',
    license='LGPLv2',
    description='Module for molcas utility functions.',
    author='Oskar Weser',
    author_email='oskar.weser@gmail.com',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    py_modules=[splitext(basename(path))[0] for path in glob('src/*.py')],
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        # complete classifier list: http://pypi.python.org/pypi?%3Aaction=list_classifiers
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: GNU Lesser General Public License v2 (LGPLv2)',
        'Operating System :: Unix',
        'Operating System :: POSIX',
        'Operating System :: Microsoft :: Windows',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: Implementation :: CPython',
        'Topic :: Utilities',
    ],
    keywords=[
        # eg: 'keyword1', 'keyword2', 'keyword3',
    ],
    python_requires='>=3.6',
    install_requires=[
        'click', 'attrs'
    ],
    #  extras_require={
    #      # eg:
    #      #   'rst': ['docutils>=0.11'],
    #      #   ':python_version=="2.6"': ['argparse'],
    #  },
    entry_points={
        'console_scripts': [
            'spat_to_spin = analyze_molcas.cli:spat_to_spin',
            'assimilate = analyze_molcas.cli:assimilate',
            'combine_orbs = analyze_molcas.cli:combine_orbs',
        ]
    },
)
