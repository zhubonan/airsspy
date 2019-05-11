#!/usr/bin/env python

from setuptools import setup, find_packages

version = '0.1.0'
if __name__ == '__main__':
    import os
    install_folder = os.path.split(__file__)[0]
    with open(os.path.join(install_folder, 'README.md')) as fh:
        long_description = fh.read()

    setup(
        name='airss-ase',
        version=version,
        url='https://www.gitlab.com/bz1/airss-ase',
        packages=find_packages(),
        install_requires=[
            'ase',
            'castepinput == 0.1.4',
            'six',
        ],
        extras_requrie={'testing': ['pytest']},
        maintainer='Bonan Zhu',
        maintainer_email='zhubonan@outlook.com',
        long_description=long_description,
    )
