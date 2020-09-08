# -*- coding: UTF-8 -*-
#   \file solvate_matchFragments.py
#   \brief This file provides ...
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
import sys                                            ### Access to the argvs
import os                                             ### For file path manipulation
import glob                                           ### Listing all files in folder

import solvate_globals                                ### Global variables
import solvate_log                                    ### Writing log
import solvate_structures                             ### Fragment reading
import solvate_maths                                  ### Computing distances and rotations

######################################################
# startLog
def matchFragments ( resList ):
    """
    ...

    Parameters
    ----------
    resList : list
        The list of all residues as returned by parseOutCoords() function.

    Returns
    -------
    NONE

    """
    ### Log progress
    solvate_log.writeLog                              ( "Starting fragment matching", 1 )
    
    ### Initialise variables
    resInd                                            = 0
    
    ### For each amino acid
    for res in resList:
        
        ### How are we matching?
        if solvate_globals._useBackboneAtoms:
        
            ### Initialise variables
            resInd                                    = resInd + 1
            
            ### Report progress
            solvate_log.writeLog                      ( "Processing residue number " + str ( resInd ) + " (" + str ( res[0] ) + ")", 2 )
            
            ### Find all input hydrated residue files
            fragPath                                  = os.path.join ( solvate_globals._resInputDir, res[0] )
            allFragFiles                              = [ fr for fr in glob.glob ( fragPath + "**/*.pdb", recursive = True ) ]
            
            ### If no frags, this is probably an error
            if len ( allFragFiles ) == 0:
                # Print error
                print ( "!!! ERROR !!! Could not find any hydrated residues for residue " + str (  res[0] ) + ". Please use the --res option to supply the path to the hydrated residues folder location." );
                solvate_log.writeLog                  ( "!!! ERROR !!! Could not find any hydrated residues for residue " + str (  res[0] ) + ". Please use the --res option to supply the path to the hydrated residues folder location.", 0 )
                
                # Terminate
                solvate_log.endLog                    ( )
                
            ### For each fragment, find distance to this residue
            distances                                 = []
            for frag in allFragFiles:
            
                # Read it in
                ( protein, waters )                   = solvate_structures.readFragment ( frag )
                
                # Compute distance between hydrated residue fragment and the processed residue
                distances.append                      ( solvate_maths.residueToResFragmentDistance ( res, protein ) )
                
            ### Find best hydrated residue fragment, if it exists
            minimalIndex                              = solvate_maths.findSmallestDistance ( distances, solvate_globals._RMSDthreshold )
            
            ### If no distance is below the threshold
            if minimalIndex == -1:
                solvate_log.writeLog                  ( "Failed to find any matching hydrated residue fragment.", 3 )
                continue
            
            ### Determine best hydrated residue fragment's secondary structure from name
            secStr                                    = dists[minInd][0][4].split( "/" )[ len( dists[minInd][0][4].split( "/" ) ) - 1 ].split( "_" )[1]
            
            ### Report progress
            logFile.write         ( "          Best backbone fitting residue fragment is " + dists[minInd][0][4] + " with RMSD of " + str( minVal ) + " and secondary structure therefore is " + str( secStr ) + "\n" )
            
            print ( distances )
            import sys; sys.exit()
            
            
        else:
        
            ### Using all fragments
            print ( "Will NOT use backbone!" )

