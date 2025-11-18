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
Import common used stuff to the model namespace
"""

from typing import Final

from .build import Buildcell
from .restools import (
    RESFile,
    TitlInfo,
    extract_res,
    format_minsep,
    get_minsep,
    get_spacegroup_atoms,
    parse_titl,
    read_res_atoms,
    read_res_pmg,
    save_airss_res,
    unique,
)
from .seed import SeedAtoms
from .utils import (
    calc_kpt_tuple_recip,
    count_pattern_in_file,
    extract_number_from_string,
    filter_out_stream,
    find_pattern_in_file,
    format_time_elapsed,
    safe_cast_float,
    safe_cast_int,
    stream_to_list,
    trim_stream,
)
from .utils import (
    unique as utils_unique,
)

__version__: Final = "0.1.4"

__all__ = [
    # Core classes
    "SeedAtoms",
    "Buildcell",
    # RES file handling
    "RESFile",
    "TitlInfo",
    "extract_res",
    "save_airss_res",
    "parse_titl",
    "read_res_atoms",
    "read_res_pmg",
    "get_spacegroup_atoms",
    "get_minsep",
    "format_minsep",
    "unique",
    # Utilities
    "trim_stream",
    "filter_out_stream",
    "calc_kpt_tuple_recip",
    "utils_unique",
    "stream_to_list",
    "safe_cast_float",
    "safe_cast_int",
    "extract_number_from_string",
    "find_pattern_in_file",
    "count_pattern_in_file",
    "format_time_elapsed",
]
