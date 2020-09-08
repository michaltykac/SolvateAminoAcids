# -*- coding: UTF-8 -*-
#   \file solvate_log.py
#   \brief This file provides the log wriiting functionality for other parts of solvate.
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
import sys                                            ### Access to the argvs and exit

import solvate_globals                                ### Global variables

######################################################
# startLog
def startLog ( ):
    """
    This function opens the logfile and writes the initial information.

    Parameters
    ----------
    NONE

    Returns
    -------
    NONE

    """
    solvate_globals._logFile                          = open ( solvate_globals._logPath, "w" )
    solvate_globals._logFile.write                    ( "Starting solvate (v" + str( solvate_globals._version ) + ") run:\n" )
    solvate_globals._logFile.write                    ( "==============================\n\n" )
    
    solvate_globals._logFile.write                    ( "CALL:\n   " )
    solvate_globals._logFile.write                    ( " ".join ( sys.argv ) + "\n\n" )
    
    solvate_globals._logFile.write                    ( "PARSED ARGS:\n" )
    solvate_globals._logFile.write                    ( "   Log path          : " + str( solvate_globals._logPath ) + "\n" )
    solvate_globals._logFile.write                    ( "   Input file        : " + str( solvate_globals._inputCoordinateFile ) + "\n\n" )
#    solvate_globals._logFile.write                    ( "   Output file       : " + str( outputName ) + "\n" )
#    solvate_globals._logFile.write                    ( "   Fragments loc.    : " + str( fragInputDir ) + "\n" )
#    solvate_globals._logFile.write                    ( "   Residues loc.     : " + str( resInputDir ) + "\n" )
#    solvate_globals._logFile.write                    ( "   RMSD threshhold   : " + str( rmsdThreshold ) + "\n" )
#    solvate_globals._logFile.write                    ( "   Clash distance    : " + str( clashDist ) + "\n" )
#    solvate_globals._logFile.write                    ( "   Clustering dist.  : " + str( minDistForClust ) + "\n" )
#    solvate_globals._logFile.write                    ( "   Only best match   : " + str( not matchAllFrags ) + "\n" )
#    solvate_globals._logFile.write                    ( "   Matched frags out.: " + str( matchedFragsOutput ) + "\n" )
#    solvate_globals._logFile.write                    ( "   Use full structure: " + str( fullStructure ) + "\n" )
#    solvate_globals._logFile.write                    ( "   No hydrogens      : " + str( noHydro ) + "\n" )
#    solvate_globals._logFile.write                    ( "   No waters         : " + str( noWater ) + "\n" )
#    solvate_globals._logFile.write                    ( "   No ligands        : " + str( noLigand ) + "\n" )
#    solvate_globals._logFile.write                    ( "   No H or H2O       : " + str( noWaterOrHydro ) + "\n" )
#    solvate_globals._logFile.write                    ( "   No H or ligand    : " + str( noLigandOrHydro ) + "\n" )
#    solvate_globals._logFile.write                    ( "   No H2O or ligand  : " + str( noLigandOrWater ) + "\n" )
#    solvate_globals._logFile.write                    ( "   Protein only      : " + str( protOnly ) + "\n" )
#    solvate_globals._logFile.write                    ( "   Pre-sort backbone : " + str( useBackbone ) + "\n\n" )

    solvate_globals._logFile.write                    ( "RUN:\n" )
    

######################################################
# endLog
def endLog ( ):
    """
    This function closes the logfile and writes the final information.

    Parameters
    ----------
    NONE

    Returns
    -------
    NONE

    """
    ### Write final message
    elapsed_time                                      = time.time () - solvate_globals._start_time
    solvate_globals._logFile.write                    ( "\nTime taken: %.3f seconds\n" % elapsed_time )
    solvate_globals._logFile.close                    ( )
    
    ### Terminate
    sys.exit                                          ( )

######################################################
# writeLog
def writeLog ( logEntry, indent = 1 ):
    """
    This function writes the supplied string into the log.

    Parameters
    ----------
    str : logEntry
        This is the string that should be written into the log.
        
    int : indent
        How much indented should this log entry be?

    Returns
    -------
    NONE

    """
    ### Initialise
    outputText                                        = " "
    
    ### Deal with indentation
    for iter in range ( 0, indent ):
        outputText                                    = outputText + "... "
        
    ### Write
    solvate_globals._logFile.write                    ( str ( outputText ) + str ( logEntry ) + "\n" )
#    print                                             ( outputText + logEntry )
