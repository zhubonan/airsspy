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
version = '0.1.1'

if __name__ == '__main__':
    import os
    install_folder = os.path.split(__file__)[0]
    with open(os.path.join(install_folder, 'README.md')) as fh:
        long_description = fh.read()

    setup(
        name='airsspy',
        version=version,
        url='https://www.gitlab.com/bz1/airsspy',
        packages=find_packages(),
        install_requires=[
            'ase',
            'castepinput == 0.1.4',
        ],
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
