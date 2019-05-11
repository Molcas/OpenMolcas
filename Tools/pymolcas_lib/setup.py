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
    name='pymolcas-lib',
    version='0.0.0',
    license='MIT license',
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
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Operating System :: Unix',
        'Operating System :: POSIX',
        'Operating System :: Microsoft :: Windows',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: Implementation :: CPython',
        'Programming Language :: Python :: Implementation :: PyPy',
        'Topic :: Utilities',
    ],
    keywords=[
        # eg: 'keyword1', 'keyword2', 'keyword3',
    ],
    python_requires='>=3.4',
    install_requires=[
        'click', 'attrs'
    ],
    #  extras_require={
    #      # eg:
    #      #   'rst': ['docutils>=0.11'],
    #      #   ':python_version=="2.6"': ['argparse'],
    #  },
    #  entry_points={
    #      'console_scripts': [
    #          'pymolcas-lib = pymolcas_lib.cli:main',
    #      ]
    #  },
)
