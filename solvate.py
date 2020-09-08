# -*- coding: UTF-8 -*-
## @package solvate
#   \file solvate.py
#   \brief Solvate is a Python language module for assigning water molecules to macromolecular structures.
#
#   ...
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

import solvate_globals                                ### Global variables
import solvate_commandLineArgs                        ### Command line arguments parsing
import solvate_log                                    ### Writing log
import solvate_structures                             ### Structures manipulation
import solvate_matchFragments                         ### Fragment matching

######################################################
# Get the args
clArgs                                                = solvate_commandLineArgs.getCLArgs ()

if clArgs.log is not None:
    solvate_globals._logPath                          = str ( clArgs.log[0] )
if clArgs.i is not None:
    solvate_globals._inputCoordinateFile              = str ( clArgs.i[0] )
if clArgs.r is not None:
    solvate_globals._RMSDthreshold                    = float ( clArgs.r[0] )
if clArgs.res is not None:
    solvate_globals._resInputDir                      = str ( clArgs.res[0] )

solvate_globals._useBackboneAtoms                     = clArgs.b

######################################################
# Start the log
solvate_log.startLog                                  ( )

######################################################
# Read in the co-ordinate file
inputCoords                                           = solvate_structures.readInCoordinates ( solvate_globals._inputCoordinateFile )

# Parse co-ordinates for solvate speed-up
resList                                               = solvate_structures.parseOutCoords ( inputCoords )

# Match fragments
solvate_matchFragments.matchFragments                 ( resList )

######################################################
# Terminate
solvate_log.endLog                                    ( )
