#!/usr/bin/env python
# -*- coding: utf-8 -*-


from setuptools import setup

setup(name='MRBIGR',
    python_requires='>3.7',
    version='1.0',
    description='MRBIGR: Mendelian Randomization-Based Inference of Genetic Regulation',
    url='https://github.com/CrazyHsu/MRBIGR.git',
    author='Feng Xu',
    author_email='',
    license='GPLv3',
    packages=['mrbigr'],
    scripts=['MRBIGR.py'],
    install_requires=[
        'cython',
        'pandas_plink',
        'numpy',
        'pandas',
        'scipy',
        'scikit-learn',
        'matplotlib'
    ]
)
