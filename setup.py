#!/usr/bin/env python
# -*- coding: utf-8 -*-
###########################################################################
# airsspy                                                               #
# Copyright (C) 2019  Bonan Zhu                                           #
#                                                                         #
# This program is free software; you can redistribute it and/or modify    #
# it under the terms of the GNU General Public License as published by    #
# the Free Software Foundation; either version 2 of the License, or       #
# (at your option) any later version.                                     #
#                                                                         #
# This program is distributed in the hope that it will be useful,         #
# but WITHOUT ANY WARRANTY; without even the implied warranty of          #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           #
# GNU General Public License for more details.                            #
#                                                                         #
# You should have received a copy of the GNU General Public License along #
# with this program; if not, write to the Free Software Foundation, Inc., #
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.             #
###########################################################################

from setuptools import setup, find_packages
version = '0.1.2'

if __name__ == '__main__':
    from os import path
    import os

    # Check if in a CI environment
    is_tagged = False
    if os.environ.get('CI_COMMIT_TAG'):
        ci_version = os.environ['CI_COMMIT_TAG']
        is_tagged = True
    elif os.environ.get('CI_JOB_ID'):
        ci_version = os.environ['CI_JOB_ID']
    else:
        # Note in CI
        ci_version = None

    # If in a CI environment, set the version accordingly
    if ci_version:
        # If this a release, check the consistency
        if is_tagged:
            assert ci_version == version, 'Inonsistency between versions'
        else:
            version = ci_version

    README_PATH = path.join(path.dirname(__file__), "README.md")
    with open(README_PATH) as fh:
        long_description = fh.read()

    setup(
        name='airsspy',
        version=version,
        url='https://github.com/zhubonan/airsspy',
        packages=find_packages(),
        install_requires=[
            'ase~=3.17',
            'castepinput~=0.1',
        ],
        description='A wrapper for using AIRSS with python and ase.',
        extras_require={
            'testing': ['pytest'],
            "pre-commit": [
                "pre-commit==1.11.0",
                "yapf==0.24.0",
            ]
        },
        maintainer='Bonan Zhu',
        maintainer_email='bon.zhu@protonmail.com',
        long_description=long_description,
        long_description_content_type='text/markdown',
    )
