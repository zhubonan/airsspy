# -*- coding: utf-8 -*-
###########################################################################
# airss-ase                                                               #
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
"""
Tools for hanlding res files
"""
import os
import ase.io


def extract_res(fname):
    """
    Extract information from res file.
    Returns a dictionary.
    The structure of he res file is not extracted
    """
    rems = []
    with open(fname) as fh:
        for line in fh:
            if 'TITL' in line:
                title = line.strip()
            if 'REM' in line:
                rems.append(line.replace("REM", "").strip())
            # Break when data starts
            if 'cell' in line:
                break
    entries = title.split()
    res = {}
    res['rem'] = rems
    res['uid'] = entries[1]
    res['P'] = float(entries[2])
    res['V'] = float(entries[3])
    res['H'] = float(entries[4])
    res['nat'] = int(entries[7])
    res['sym'] = entries[8]
    res['fname'] = fname
    return res


def save_airss_res(atoms, info_dict, fname=None, force_write=False):
    """
    Save the relaxed structure in res format which is compatible with the
    ``cryan`` program.
    """

    # Prepare output file
    if fname is None:
        fname = info_dict['uid'] + '.res'
    if os.path.isfile(fname) and not force_write:
        raise FileExistsError(
            "Switch on force_write to overwrite existing files")

    P, V, H = info_dict['P'], info_dict['V'], info_dict['H']
    # Get number of atoms, spin
    nat, sg = info_dict['nat'], info_dict['sym']
    # Construct title line
    pvh = " {:.3f} {:.3f} {:.6f} ".format(P, V, H)
    title = 'TITL ' + info_dict['uid'] + ' ' + pvh + ' 0 ' + ' 0 ' + ' ' + str(
        nat) + ' ' + sg + ' n - 1\n'
    rems = info_dict.get('rem', [])

    # Write to the top of res file
    restmp = info_dict['uid'] + '.rtmp'
    ase.io.write(restmp, atoms, format="res")
    with open(restmp) as resin:
        with open(fname, 'w') as fout:
            # Write the title
            fout.write(title + '\n')
            for line in rems:
                fout.write('REM ' + line + '\n')
            # Write the data lines
            for n, line in enumerate(resin):
                if n > 0:
                    fout.write(line)

    os.remove(restmp)
    return
