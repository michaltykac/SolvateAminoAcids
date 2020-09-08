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
# matchFragments
def matchFragments ( resList ):
    """
    ...

    Parameters
    ----------
    list : resList
        The list of all residues as returned by parseOutCoords() function.

    Returns
    -------
    NONE

    """
    ### Log progress
    solvate_log.writeLog                              ( "Starting fragment matching", 1 )
    
    ### Initialise variables
    resInd                                            = 0
    resMatches                                        = []
    
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
            
            ### Find the secondary structure
            secStr                                    = findResFragSecondaryStructure ( res, allFragFiles )
            if ( secStr == "NOT_FOUND" ):
                continue;
                
            ### Find the rotamer
            rotamer                                   = findResFragRotamer ( res, allFragFiles, secStr )
            print ( rotamer )
            
        else:
        
            ### Using all fragments
            print ( "Will NOT use backbone!" )


######################################################
# findResFragSecondaryStructure()
def findResFragSecondaryStructure ( res, allFragFiles ):
    """
    This function determines the secondary structure of the best backbone-atoms-fitting hydrated-residue
    fragment and returns it.

    Parameters
    ----------
    list : res
        A list of all atoms parsed for this residue.
        
    list : allFragFiles
        List of all fragments found for this residue

    Returns
    -------
    str : secStr
        The letter coding the secondary structure of the best fitting fragment or "NOT_FOUND" if no match was found.

    """
    ### If no frags, this is probably an error
    if len ( allFragFiles ) == 0:
        # Print error
        solvate_log.writeLog                          ( "!!! ERROR !!! Could not find any hydrated residues for residue " + str (  res[0] ) + ". Please use the --res option to supply the path to the hydrated residues folder location.", 0 )
        
        # Terminate
        solvate_log.endLog                            ( )
        
    ### For each fragment, find distance to this residue
    distances                                         = []
    for frag in allFragFiles:
            
        # Read it in
        ( protein, waters )                           = solvate_structures.readFragment ( frag )
        
        # Compute distance between hydrated residue fragment and the processed residue
        distances.append                              ( solvate_maths.residueToResFragmentDistance_backbone ( res, protein ) )
        
    ### Find best hydrated residue fragment in terms of backbone, if it exists
    minimalIndex                                      = solvate_maths.findSmallestDistance ( distances, solvate_globals._RMSDthreshold )
    
    ### If no distance is below the threshold
    if minimalIndex == -1:
        solvate_log.writeLog                          ( "Failed to find any matching hydrated residue fragment.", 3 )
        return                                        ( "NOT_FOUND" )
    
    ### Determine best hydrated residue fragment's secondary structure from name
    secStr                                            = allFragFiles[minimalIndex].split( "/" )[len ( allFragFiles[minimalIndex].split( "/" ) ) - 1].split ( "_" )[1]
            
    ### Report progress
    solvate_log.writeLog                              ( "Best backbone fitting residue fragment is " + allFragFiles[minimalIndex].split( "/" )[len ( allFragFiles[minimalIndex].split( "/" ) ) - 1] + " with RMSD of " + str( distances[minimalIndex] ) + " and secondary structure therefore is " + str( secStr ), 3 )
    
    ### Done
    return                                            ( secStr )

######################################################
# findResFragRotamer ()
def findResFragRotamer ( res, allFragFiles, secStr ):
    """
    This function determines the rotamer of the best side-chain-atoms fitting hydrated-residue
    fragment and returns it.

    Parameters
    ----------
    list : res
        A list of all atoms parsed for this residue.
        
    list : allFragFiles
        List of all fragments found for this residue
        
    str : secStr
        The secondary structure for which the rotamer is to be determined.

    Returns
    -------
    str : rotamer
        The string coding the rotamer of the best fitting fragment, "NOT_FOUND" if no match was found or "NOT_INCLUDED" for ALA and GLY.

    """
    ### Initialise variables
    rotamer                                           = ""
    
    ### Find rotamer
    if ( res[0] != "ALA" ) and ( res[0]  != "GLY" ):

        ### Re-initialise variables
        distances                                     = []
        forNames                                      = []

        ### Compare to each fragment
        for frag in allFragFiles:

            ### Ignore mismatching secondary structure fragments
            if frag.split( "/" )[ len( frag.split( "/" ) ) - 1 ].split( "_" )[1] != secStr:
                continue
                
            ### Parse the fragment ( again, this could be improved, but would make it harder to read )
            ( protein, waters )                        = solvate_structures.readFragment ( frag )

            ### Find distance
            distances.append                          ( solvate_maths.residueToResFragmentDistance_rotamer ( res, protein ) )
            
        ### Find best hydrated residue fragment in terms of rotamer, if it exists
        minimalIndex                                  = solvate_maths.findSmallestDistance ( distances, solvate_globals._RMSDthreshold )

        ### If no distance is below the threshold
        if minimalIndex == -1:
            solvate_log.writeLog                      ( "Failed to find any matching hydrated residue fragment.", 3 )
            return                                    ( "NOT_FOUND" )
            
        ### Find best rotamer
        rotamer                                       = allFragFiles[minimalIndex].split( "/" )[ len( allFragFiles[minimalIndex].split( "/" ) ) - 1 ].split( "." )[0].split( "_" )[2]

        ### Report progress
        solvate_log.writeLog                          ( "Best rotamer fitting residue fragment is " + allFragFiles[minimalIndex].split( "/" )[len ( allFragFiles[minimalIndex].split( "/" ) ) - 1] + " with RMSD of " + str( distances[minimalIndex] ) + " and the rotamer therefore is " + str( rotamer ), 3 )

    else:

        ### Report progress
        solvate_log.writeLog                          ( "This is either ALA or GLY residue. There are no rotamers for these residues...", 3 )
        
        ### Done
        return                                        ( "NOT_INCLUDED" )
    
    ### Done
    return                                            ( rotamer )
