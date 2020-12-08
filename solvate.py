# -*- coding: UTF-8 -*-
## @package solvate
#   \file solvate.py
#   \brief ...
#
#   This file
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
#   \version   0.1.0
#   \date      DEC 2020
######################################################

######################################################
# Imports

import solvate

######################################################
# Create and set global settings
settings                                              = solvate.globalSettings ( )
settings.parseCommandLineArguments                    ( )

######################################################
# Start the log
solvate.startLog                                      ( settings )

######################################################
# Read in the structure and parse out required information
resList                                               = solvate.parseInputCoordinates ( settings )

######################################################
# Read in fragments and residues
fragFragments                                         = solvate.getAllFragmentFragments ( settings )
resFragments                                          = {}
if settings.useBackboneAtoms:
    resFragments                                      = solvate.getAllResidueFragments ( settings )

######################################################
# Match fragments to each residue
matchedFrags                                          = solvate.matchFragments ( resList, fragFragments, resFragments,  settings )

######################################################
# Predict all waters for co-ordinates
waters                                                = solvate.predictWaters ( resList, matchedFrags, fragFragments, settings )

######################################################
# Remove clashing waters
( watersNoClash, waterLabels )                        = solvate.removeClashes ( waters, settings )

######################################################
# Cluster waters
waterClusters                                         = solvate.clusterWaters ( watersNoClash, waterLabels, settings )

######################################################
# Combine clusters and add waters
solvate.combineAndAddWaters                           ( waterClusters, waterLabels, settings )

######################################################
# Write out hydrated co-ordinates
solvate.writeOutStructures                            ( waterLabels, settings )

######################################################
# Terminate
solvate.endLog                                        ( settings )
