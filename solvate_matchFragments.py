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
    list : resMatches
        A list with entry for each residue, where the entry is a list of all hydrated fragments that
        were matched to the particular residue using the current settings (backbone, bestOnly,
        threshold, ...)

    """
    ### Log progress
    solvate_log.writeLog                              ( "Starting fragment matching", 1 )
    
    ### Initialise variables
    resInd                                            = 0
    resMatches                                        = []
    
    ### Read in fragments and residues
    fragFragments                                     = solvate_structures.getAllFragmentFragments (  )
    resFragments                                      = {}
    if solvate_globals._useBackboneAtoms:
        resFragments                                  = solvate_structures.getAllResidueFragments (  )
    
    ### For each amino acid
    for res in resList:

        ### How are we matching?
        if solvate_globals._useBackboneAtoms:
        
            ### Initialise variables
            resInd                                    = resInd + 1
            
            ### Report progress
            solvate_log.writeLog                      ( "Processing residue number " + str ( resInd ) + " (" + str ( res[0] ) + ")", 2 )
            
            ### Match the residue to hydrated fragments using backbone step
            resMatches.append                         ( matchFragmentsUsingBackbone ( res, resFragments, fragFragments ) )
            
        else:
        
            ### Initialise variables
            resInd                                    = resInd + 1
            
            ### Report progress
            solvate_log.writeLog                      ( "Processing residue number " + str ( resInd ) + " (" + str ( res[0] ) + ")", 2 )
            
            ### Match the residue to hydrated fragments using backbone step
            resMatches.append                         ( matchFragmentsDirect ( res, fragFragments ) )
            
    ### Done
    return                                            ( resMatches )


######################################################
# findResFragSecondaryStructure()
def findResFragSecondaryStructure ( res, resFragments ):
    """
    This function determines the secondary structure of the best backbone-atoms-fitting hydrated-residue
    fragment and returns it.

    Parameters
    ----------
    list : res
        A list of all atoms parsed for this residue.
        
    dictionary : resFragments
        A dictionary holding the parsed data for each hydrated residue fragment in terms of both
        protein and waters atoms.

    Returns
    -------
    list : secStr
        List of the letters coding the secondary structure of the fitting fragments or empty list if no match
        was found. The number of returned values depends on whether the best fragment only or all passing
        fragments are to be returned.

    """
    ### For each fragment, find distance to this residue
    distances                                         = []
    fragNames                                         = []
    for frag in resFragments.keys():
        
        # Ignore entries with different amino acid
        if frag.split( os.path.sep )[ len ( frag.split( os.path.sep ) ) - 1 ].split( "_" )[0] != res[0]:
            continue
        
        # Compute distance between hydrated residue fragment and the processed residue
        distances.append                              ( solvate_maths.residueToResFragmentDistance_backbone ( res, resFragments[frag]["protein"] ) )
        fragNames.append                              ( frag )
        
    ### Are we searching for best, or all passing?
    if solvate_globals._bestFragmentOnly:
        
        ### Find best hydrated residue fragment in terms of backbone, if it exists
        minimalIndex                                  = solvate_maths.findSmallestDistance ( distances, solvate_globals._RMSDthreshold )
        
        ### If no distance is below the threshold
        if minimalIndex == -1:
            solvate_log.writeLog                      ( "Failed to find any matching hydrated residue fragment.", 3 )
            return                                    ( [] )
        
        ### Determine best hydrated residue fragment's secondary structure from name
        secStr                                        = fragNames[minimalIndex].split( os.path.sep )[len ( fragNames[minimalIndex].split( os.path.sep ) ) - 1].split ( "_" )[1]
        
        ### Report progress
        solvate_log.writeLog                          ( "Best backbone fitting residue fragment is " + fragNames[minimalIndex].split( os.path.sep )[len (    fragNames[minimalIndex].split( os.path.sep ) ) - 1] + " with RMSD of " + str( distances[minimalIndex] ) + " and secondary structure therefore is " + str(    secStr ), 3 )
    
        ### Done
        return                                        ( [ secStr ] )
        
    else:
        
        ### Initialise variables
        secStr                                        = []
        
        ### Find passing indices
        indices                                       = solvate_maths.findPassingDistances ( distances, solvate_globals._RMSDthreshold )
        
        ### If no distance is below the threshold
        if indices == -1:
            solvate_log.writeLog                      ( "Failed to find any matching hydrated residue fragment.", 3 )
            return                                    ( [] )
        
        ### For each index
        for ind in indices:
            
            # Determine the secondary structure for the fragment
            secStr.append                             ( fragNames[ind].split( os.path.sep )[len ( fragNames[ind].split( os.path.sep ) ) - 1].split ( "_" )[1] )
            
            # Report progress
            solvate_log.writeLog                      ( "A backbone fitting residue fragment " + fragNames[ind].split( os.path.sep )[len (    fragNames[ind].split( os.path.sep ) ) - 1] + " with RMSD of " + str( distances[ind] ) + " and secondary structure " + str( fragNames[ind].split( os.path.sep )[len ( fragNames[ind].split( os.path.sep ) ) - 1].split ( "_" )[1] ) + " found.", 3 )
        
        ### Done ( returning unique results )
        return                                        ( list ( set ( secStr ) ) )

######################################################
# findResFragRotamer ()
def findResFragRotamer ( res, resFragments, secStr ):
    """
    This function determines the rotamer of the best side-chain-atoms fitting hydrated-residue
    fragment and returns it.

    Parameters
    ----------
    list : res
        A list of all atoms parsed for this residue.
        
    dictionary : resFragments
        A dictionary holding the parsed data for each hydrated residue fragment in terms of both
        protein and waters atoms.
        
    list : secStr
        The secondary structure for which the rotamer is to be determined.

    Returns
    -------
    str : rotamer
        List of the letters coding the rotamers of the fitting fragments followed by underscore and the
        appropriate secondary structure letter or empty list if no match was found. The number of returned
        values depends on whether the best fragment only or all passing fragments are to be returned. For ALA
        and GLY, "NOT_INCLUDED" will be returned instead.

    """
    ### Initialise variables
    rotamer                                           = ""
    
    ### Find rotamer
    if ( res[0] != "ALA" ) and ( res[0]  != "GLY" ):

        ### Re-initialise variables
        distances                                     = []
        forNames                                      = []

        ### Compare to each fragment
        for frag in resFragments.keys():
        
            ### Ignore mismatching amino acids
            if frag.split(os.path.sep)[ len ( frag.split(os.path.sep) ) - 1 ].split("_")[0] != res[0]:
                continue

            ### Ignore mismatching secondary structure fragments
            if frag.split( os.path.sep )[ len( frag.split( os.path.sep ) ) - 1 ].split( "_" )[1] not in secStr:
                continue

            ### Find distance
            distances.append                          ( solvate_maths.residueToResFragmentDistance_rotamer ( res, resFragments[frag]["protein"] ) )
            forNames.append                           ( frag )
            
        ### Are we searching for best, or all passing?
        if solvate_globals._bestFragmentOnly:
            
            # Find best hydrated residue fragment in terms of rotamer, if it exists
            minimalIndex                              = solvate_maths.findSmallestDistance ( distances, solvate_globals._RMSDthreshold )
    
            # If no distance is below the threshold
            if minimalIndex == -1:
                solvate_log.writeLog                  ( "Failed to find any matching hydrated residue fragment.", 3 )
                return                                ( [] )
            
            # Find best rotamer
            rotamer                                   = forNames[minimalIndex].split( os.path.sep )[ len( forNames[minimalIndex].split( os.path.sep ) ) - 1 ].split( "."    )[0].split( "_" )[2] + "_" + forNames[minimalIndex].split( os.path.sep )[ len( forNames[minimalIndex].split( os.path.sep ) ) - 1 ].split("_")[1];
    
            ### Report progress
            solvate_log.writeLog                      ( "Best rotamer fitting residue fragment is " + forNames[minimalIndex].split( os.path.sep )[len ( forNames[minimalIndex].split( os.path.sep ) ) - 1] + " with RMSD of " + str( distances[minimalIndex] ) + " and the rotamer therefore is " + str( rotamer ), 3 )
            
            ### Done
            return                                    ( [ rotamer ] )
            
        else:
        
            # Initialise variables
            rotamer                                   = []
            
            # Find passing indices
            indices                                   = solvate_maths.findPassingDistances ( distances, solvate_globals._RMSDthreshold )
            
            # If no distance is below the threshold
            if indices == -1:
                solvate_log.writeLog                  ( "Failed to find any matching hydrated residue fragment.", 3 )
                return                                ( [] )
            
            ### For each index
            for ind in indices:
                
                # Determine the secondary structure for the fragment
                rotamer.append                        ( forNames[ind].split( os.path.sep )[ len( forNames[ind].split( os.path.sep ) ) - 1 ].split( "."    )[0].split( "_" )[2] + "_" + forNames[ind].split( os.path.sep )[ len( forNames[ind].split( os.path.sep ) ) - 1 ].split("_")[1] )
                
                # Report progress
                solvate_log.writeLog                  ( "A rotamer fitting residue fragment " + forNames[ind].split( os.path.sep )[len ( forNames[ind].split( os.path.sep ) ) - 1] + " with RMSD of " + str( distances[ind] ) + " and the rotamer " + str( forNames[ind].split( os.path.sep )[ len( forNames[ind].split( os.path.sep ) ) - 1 ].split( "."    )[0].split( "_" )[2] ) + " found.", 3 )

            ### Done ( returning unique results )
            return                                    ( list ( set ( rotamer ) ) )

    else:

        ### Report progress
        solvate_log.writeLog                          ( "This is either ALA or GLY residue. There are no rotamers for these residues...", 3 )
        
        ### Done
        return                                        ( "NOT_INCLUDED" )
    
    ### Done
    return                                            ( rotamer )

######################################################
# matchFragmentsUsingBackbone ()
def matchFragmentsUsingBackbone ( res, resFragments, fragFragments ):
    """
    ...

    Parameters
    ----------
    list : res
        A list of all atoms parsed for this residue.
        
    dictionary : resFragments
        A dictionary holding the parsed data for each hydrated residue fragment in terms of both
        protein and waters atoms.
        
    dictionary : fragFragments
        A dictionary holding the parsed data for each hydrated fragment in terms of both
        protein and waters atoms.

    Returns
    -------
    list : matched
        A list of all fragments that were matched to the input residue using the global criteria ( RMSD threshold,
        onlyBest, ... )

    """
    ### Find the secondary structure
    secStr                                            = findResFragSecondaryStructure ( res, resFragments )

    ### Find the rotamer for all passing secondary structures
    rotamer                                           = []
    fragList                                          = []
    for ss in secStr:
    
        # Search for rotamers using this secondary structure
        rotamerHlp                                    = findResFragRotamer ( res, resFragments, secStr )
        
        # Check for ALA and GLY
        if rotamerHlp == "NOT_INCLUDED":
        
            # Initialise variables
            fragList                                  = []
            
            # Report progress
            for fl in fragFragments.keys():
            
                # Ignore fragments with wrong amino acid type
                if fl.split( os.path.sep )[ len ( fl.split( os.path.sep ) ) - 1 ].split("_")[0] != res[0]:
                    continue
                    
                # Ignore fragments with wrong secondary structure
                if fl.split( os.path.sep )[ len ( fl.split( os.path.sep ) ) - 1 ].split("_")[1] != ss:
                    continue
            
                solvate_log.writeLog                  ( "Matched residue to fragment " + str ( fl ) + ".", 3 )
                fragList.append                       ( fl )
            
        else:
        
            # Save
            rotamer                                   = rotamer + rotamerHlp

    ### If ALA or GLY, this is it!
    if len ( fragList ) > 0:

        # Save to output
        return                                        ( fragList )

    ### Other residue. Well, deal with rotamers then!
    for rot in rotamer:

        # Report progress and save
        for fl in fragFragments.keys():
        
            # Ignore fragments with wrong amino acid type
            if fl.split( os.path.sep )[ len ( fl.split( os.path.sep ) ) - 1 ].split("_")[0] != res[0]:
                continue
                
            # Ignore fragments with wrong secondary structure
            if fl.split( os.path.sep )[ len ( fl.split( os.path.sep ) ) - 1 ].split("_")[1] != rot.split( "_" )[1]:
                continue
                
            # Ignore fragments with wrong rotamer
            if fl.split( os.path.sep )[ len ( fl.split( os.path.sep ) ) - 1 ].split("_")[2] != rot.split( "_" )[0]:
                continue

            # Passed! Now write up and save
            solvate_log.writeLog                      ( "Matched residue to fragment " + str ( fl ) + ".", 3 )
            fragList.append                           ( fl )
            
    ### Save to output
    return                                            ( fragList )

######################################################
# matchFragmentsDirect ()
def matchFragmentsDirect ( res, fragFragments ):
    """
    ...

    Parameters
    ----------
    list : res
        A list of all atoms parsed for this residue.
        
    dictionary : fragFragments
        A dictionary holding the parsed data for each hydrated fragment in terms of both
        protein and waters atoms.

    Returns
    -------
    list : matched
        A list of all fragments that were matched to the input residue using the global criteria ( RMSD threshold,
        onlyBest, ... )

    """
    ### Initialise variables
    distances                                         = []
    forNames                                          = []
    
    ### For each fragment, get the distances
    for fl in fragFragments.keys():
    
        # Ignore fragments with wrong amino acid type
        if fl.split( os.path.sep )[ len ( fl.split( os.path.sep ) ) - 1 ].split("_")[0] != res[0]:
            continue
        
        distances.append                              ( solvate_maths.residueToResFragmentDistance_direct ( res, fragFragments[fl]["protein"] ) )
        forNames.append                               ( fl )
        
    ### Are we searching for best, or all passing?
    if solvate_globals._bestFragmentOnly:
        
        # Find best hydrated residue fragment in terms of rotamer, if it exists
        minimalIndex                                  = solvate_maths.findSmallestDistance ( distances, solvate_globals._RMSDthreshold )
    
        # If no distance is below the threshold
        if minimalIndex == -1:
            solvate_log.writeLog                      ( "Failed to find any matching hydrated residue fragment.", 3 )
            return                                    ( [] )
        
        # Report progress
        solvate_log.writeLog                          ( "Best fragment fitting the residue is " + forNames[minimalIndex].split( os.path.sep )[len ( forNames[minimalIndex].split( os.path.sep ) ) - 1] + " with RMSD of " + str( distances[minimalIndex] ) + ".", 3 )
        
        ### Done
        return                                        ( [ forNames[minimalIndex] ] )
        
    else:
    
        # Initialise variables
        mtchs                                         = []
        
        # Find passing indices
        indices                                       = solvate_maths.findPassingDistances ( distances, solvate_globals._RMSDthreshold )
        
        # If no distance is below the threshold
        if indices == -1:
            solvate_log.writeLog                      ( "Failed to find any matching hydrated residue fragment.", 3 )
            return                                    ( [] )
        
        ### For each index
        for ind in indices:
            
            # Determine the secondary structure for the fragment
            mtchs.append                              ( forNames[ind] )
            
            # Report progress
            solvate_log.writeLog                      ( "A fragment " + forNames[ind].split( os.path.sep )[len ( forNames[ind].split( os.path.sep ) ) - 1] + " fits the residue with RMSD of " + str( distances[ind] ) + " was found.", 3 )

        ### Done ( returning unique results )
        return                                        ( list ( set ( mtchs ) ) )
