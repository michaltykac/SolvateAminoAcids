# -*- coding: UTF-8 -*-
#   \file solvate_globals.py
#   \brief This file provides global variables to all other files.
#
#   This file provides the global_settings class, which contains all the settings and
#   their default values for the rest of the software. This class is used to transfer
#   and make available all the settings to the rest of the files and functions.
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

import sys                                            ### Terminating if fail
import os                                             ### Creating directories
import time                                           ### For timing
import gemmi                                          ### For checking its instances
import solvate.solvate_commandLineArgs as cl          ### For command line parsing

######################################################
# Global variables class
class globalSettings:
    """
    This class holds all the settings variables and makes for simple one-place setting
    of these as well passing these to any function that requires them.
    
    """
    def __init__( self ):
        """
        This is the constructor for the class. It setls all the default values and
        returns an instance of the class.

        Parameters
        ----------
        NONE

        Returns
        -------
        globalSettings: settings
            An instance of the globalSettings class.

        """
        ### Start timing
        self.start_time                               = time.time()
        
        ### Set the version
        self._version                                 = "0.0.2"

        ### Set verbosity
        self.verbose                                  = 1

        ### Set log location and file holder
        self.logPath                                  = "solvate_log.txt"
        self.logFile                                  = None
        
        ### Set input structure
        self.inputCoordinateFile                      = ""
        self.inputCoordinateStructure                 = None
        self.inputCoordinateStructureNoHydro          = None
        self.inputCoordinateStructureNoWaters         = None
        self.inputCoordinateStructureNoLigand         = None
        
        ### Set the hydrated data paths
        self.resInputDir                              = "resData"
        self.fragInputDir                             = "fragData"
        self.matchedFragsPath                         = ""
        
        ### List all supported amino acids
        self.aaTypes                                  = [ "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU",
                                                          "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
                                                          "PRO", "SER", "THR", "TRP", "TYR", "VAL" ]
        
        ### Set the fragment matching options
        self.useBackboneAtoms                         = False
        self.bestFragmentOnly                         = False
        
        ### Set the fragment matching threshold
        self.RMSDthreshold                            = 0.5
        
        ### Set the water clashing options
        self.clashWithinRadius                        = 2.0
        self.noFullStructure                          = False
        self.noHydro                                  = False
        self.noWaters                                 = False
        self.noLigand                                 = False
        
        ### Set the clustering threshold
        self.maxClustDist                             = 1.8

    def setLogPath ( self, logLoc ):
        """
        This mutator function sets the log file path variable.

        Parameters
        ----------
        str : logPath
            The new value (path) for the log file.

        Returns
        -------
        NONE

        """
        self.logPath                                  = logLoc
        
    def setCoordInputFile ( self, inCoordFile ):
        """
        This mutator function sets the file path for the co-ordinate file to be processed by solvate.

        Parameters
        ----------
        str : inCoordFile
            The new value (path) for the co-ordinate file to be processed by solvate.

        Returns
        -------
        NONE

        """
        self.inputCoordinateFile                      = inCoordFile
        
    def setRMSDThreshold ( self, thres ):
        """
        This mutator function sets the RMSD threshold for fragment matching.

        Parameters
        ----------
        float : thres
            The new value (float) for the RMSD threshold for matching fragments.

        Returns
        -------
        NONE

        """
        self.RMSDthreshold                            = thres
        
    def setResiduesPath ( self, resPath ):
        """
        This mutator function sets the path to the folder containing the hydrated residues.

        Parameters
        ----------
        str : resPath
            The new value (path) for the folder containing the hydrated residue files.

        Returns
        -------
        NONE

        """
        self.resInputDir                              = resPath
        
    def setFragmentsPath ( self, fragPath ):
        """
        This mutator function sets the path to the folder containing the hydrated fragments.

        Parameters
        ----------
        str : fragPath
            The new value (path) for the folder containing the hydrated fragment files.

        Returns
        -------
        NONE

        """
        self.fragInputDir                             = fragPath
        
    def setVerbose ( self, verb ):
        """
        This mutator function sets the verbosity settings.

        Parameters
        ----------
        int : verb
            The new value (int) for the verbosity of the standard output.

        Returns
        -------
        NONE

        """
        self.verbose                                  = verb
        
    def setBestFragOnly ( self, bfOnly ):
        """
        This mutator function sets the fragment matching settings regarding whether only the
        best fragment or all threshold passing fragments should be used.

        Parameters
        ----------
        boolean : bfOnly
            The new value (boolean) for whether to use only the best matched fragment or all
            fragments passing the RMSD threshold for fragment matching.

        Returns
        -------
        NONE

        """
        self.bestFragmentOnly                         = bfOnly
        
    def setUseBackbone ( self, useBack ):
        """
        This mutator function sets the fragment matching option regarding whether the
        secondary structure and rotamer should be determined to decide which fragments
        to attempt to match, or whether just matching all fragments and using distances
        is preferred.

        Parameters
        ----------
        boolean : bfOnly
            The new value (boolean) for whether to use the secondary structure and rotamer
            information derived from hydrated residue mathces should be used to determine which
            fragments to match, or whether all fragments should be tried.

        Returns
        -------
        NONE

        """
        self.useBackboneAtoms                         = useBack
        
    def setInputStructure ( self, struct ):
        """
        This mutator function sets the supplied gemmi.Structure object as the structure which
        should be hydrated.

        Parameters
        ----------
        gemmi.Structure : struct
            The gemmi.Structure object for the structure to be hydrated.

        Returns
        -------
        NONE

        """
        ### Check that this is a gemmi object!
        if not isinstance ( struct, gemmi.Structure ):
            sys.exit                                  ( 'Error: The supplied structure ' + str( struct ) + ' is not an instance of gemmi.Structure. Please supply a file parsed by gemmi.' )
            
        ### If so, set
        self.inputCoordinateStructure                 = struct
        
    def setMatchedFragsPath ( self, mfPath ):
        """
        This mutator function sets the path to where hydrated fragment - residue matches should
        be written into.

        Parameters
        ----------
        str : mfPath
            The path to where to save the matches.

        Returns
        -------
        NONE

        """
        ### If so, set
        self.matchedFragsPath                         = mfPath
        
        ### If not empty, create unless exists
        if ( self.matchedFragsPath != "" ) and ( not os.path.isdir ( self.matchedFragsPath ) ):
            os.makedirs                               ( self.matchedFragsPath )
            
    def setClashRadius ( self, clRad ):
        """
        This mutator function sets the radius withing which no atom can be for no clash to be
        detected.

        Parameters
        ----------
        str : clRad
            The radius which must be atom free if there is to be no clash.

        Returns
        -------
        NONE

        """
        ### Set
        self.clashWithinRadius                        = clRad
        
    def setNoFullStructure ( self, nfStr ):
        """
        This mutator function sets the switch deciding whether the full structure contents should
        be used to detect new water molecule clashes.

        Parameters
        ----------
        boolean : nfStr
            The new value for this switch variable.

        Returns
        -------
        NONE

        """
        ### Set
        self.noFullStructure                          = nfStr
        
    def setNoHydro ( self, nHyd ):
        """
        This mutator function sets the switch deciding whether the hydrogen atoms should
        be used to detect new water molecule clashes.

        Parameters
        ----------
        boolean : nHyd
            The new value for this switch variable.

        Returns
        -------
        NONE

        """
        ### Set
        self.noHydro                                  = nHyd
        
    def setNoWaters ( self, nWat ):
        """
        This mutator function sets the switch deciding whether the already present water
        molecules should be used to detect new water molecule clashes.

        Parameters
        ----------
        boolean : nWat
            The new value for this switch variable.

        Returns
        -------
        NONE

        """
        ### Set
        self.noWaters                                 = nWat
        
    def setNoLigand ( self, nLig ):
        """
        This mutator function sets the switch deciding whether the ligand
        molecules should be used to detect new water molecule clashes.

        Parameters
        ----------
        boolean : nLig
            The new value for this switch variable.

        Returns
        -------
        NONE

        """
        ### Set
        self.noLigand                                 = nLig
        
    def setMaxClusterDist ( self, mcd ):
        """
        This mutator function sets the maximum distance that two water molecules can have
        for them to be considered for being in the same water cluster.

        Parameters
        ----------
        float : mcd
            The maximum distance between two waters in the same water cluster.

        Returns
        -------
        NONE

        """
        ### Set
        self.maxClustDist                             = mcd

    def parseCommandLineArguments ( self ):
        """
        This function parses the command line arguments and sets the values in an instance of
        the class accordingly.

        Parameters
        ----------
        NONE

        Returns
        -------
        NONE

        """
        ### Parse the command line arguments
        clArgs                                        = cl.getCLArgs ( self._version )

        ### Save the values, if supplied
        if clArgs.log is not None:
            self.setLogPath                           ( str ( clArgs.log[0] ) )
            
        if clArgs.i is not None:
            self.setCoordInputFile                    ( str ( clArgs.i[0] ) )
            
        if clArgs.r is not None:
            self.setRMSDThreshold                     ( float ( clArgs.r[0] ) )

        if clArgs.res is not None:
            self.setResiduesPath                      ( str ( clArgs.res[0] ) )

        if clArgs.frag is not None:
            self.setFragmentsPath                     ( str ( clArgs.frag[0] ) )

        if clArgs.verbose is not None:
            self.setVerbose                           ( int ( clArgs.verbose[0] ) )
            
        if clArgs.matchedFrags is not None:
            self.setMatchedFragsPath                  ( str ( clArgs.matchedFrags[0] ) )
            
        if clArgs.clRad is not None:
            self.setClashRadius                       ( float ( clArgs.clRad[0] ) )
            
        if clArgs.maxClDist is not None:
            self.setMaxClusterDist                    ( float ( clArgs.maxClDist[0] ) )

        self.setBestFragOnly                          ( clArgs.bestOnly )
        self.setUseBackbone                           ( clArgs.b )
        self.setNoFullStructure                       ( clArgs.noFullStr )
        self.setNoHydro                               ( clArgs.noHydro )
        self.setNoWaters                              ( clArgs.noWaters )
        self.setNoLigand                              ( clArgs.noLigand )

    def getVersion ( self ):
        """
        This acessor function returns the version of the tool.

        Parameters
        ----------
        NONE

        Returns
        -------
        str : version
            The version of the code and thus the tool.

        """
        return                                        ( self._version )
