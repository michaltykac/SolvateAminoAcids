# -*- coding: UTF-8 -*-
#   \file solvate_globals.py
#   \brief This file provides global variables to all other files.
#
#   Copyright by the Authors and individual contributors. All rights reserved.
#
#   Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
#       1) Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#       2) Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
#       3) Neither the name of Michal Tykac nor the names of this code's contributors may be used to endorse or promote products derived from this software without specific prior written permission.
#
#   This software is provided by the copyright holder and contributors "as is" and any express or implied warranties, including, but not limitted to, the implied warranties of merchantibility and fitness for a particular purpose are disclaimed. In no event shall the copyright owner or the contributors be liable for any direct, indirect, incidental, special, exemplary, or consequential damages (including, but not limitted to, procurement of substitute goods or services, loss of use, data or profits, or business interuption) however caused and on any theory of liability, whether in contract, strict liability or tort (including negligence or otherwise) arising in any way out of the use of this software, even if advised of the possibility of such damage.
#
#   \author    Michal Tykac
#   \author    Lada Biedermannová
#   \author    Jiří Černý
#   \version   0.0.1
#   \date      SEP 2020
######################################################

######################################################
# Imports

import time                                           ### For timing

######################################################
# Global variables

_version                                              = "0.0.1"

_logPath                                              = "solvate_log.txt"
_logFile                                              = None

_inputCoordinateFile                                  = ""

_resInputDir                                          = "resData"

_start_time                                           = time.time()

_aaTypes                                              = [ "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU",
                                                          "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
                                                          "PRO", "SER", "THR", "TRP", "TYR", "VAL" ]

_useBackboneAtoms                                     = False
_bestFragmentOnly                                     = False

_RMSDthreshold                                        = 0.5
