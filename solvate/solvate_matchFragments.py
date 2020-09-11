# -*- coding: UTF-8 -*-
#   \file solvate_matchFragments.py
#   \brief This file provides fragment matching algorithms.
#
#   This file contains the functions and algorithms required for fragment matching to a macromolecular
#   model using several different approaches. These include matching all fragments and then determining the
#   matches using the distances, or firstly determining the residue secondary structure and rotamer using
#   matches to hydrated residues and then using this information matching the hydrated fragments. Also, the
#   algorithms here can use either only the best matched fragment, or all fragments passing the RMSD distance
#   threshold.
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
#   \version   0.0.2
#   \date      SEP 2020
######################################################

######################################################
# Imports

import time                                           ### For timing
import sys                                            ### Access to the argvs
import os                                             ### For file path manipulation
import glob                                           ### Listing all files in folder

import solvate.solvate_log as solvate_log             ### Writing log
import solvate.solvate_structures as solvate_structures ### Fragment reading
import solvate.solvate_maths as solvate_maths         ### Computing distances and rotations

######################################################
# matchFragments
def matchFragments ( resList, fragFragments, resFragments, settings ):
    """
    This function takes all supplied residues and finds the optimal/threshold passing
    hydrated fragments for each of the residues. It can use either only the best aligned
    fragment, or all passing the threshold.
    
    Moreover, the function can either align all fragments and decide from there, or first
    determine the secondary structure and rotamer (again, either the best or all passing)
    and use this information to decide on the best hydrated fragment(s) to return.
    
    It then returns a list of aligned fragments for each residue along with the procrustes
    transformations required to obtain this optimal alignment.

    Parameters
    ----------
    list : resList
        The list of all residues as returned by parseOutCoords() function.
        
    list : fragFragments
        The list of all found and parsed hydrated fragment files with their
        protein and waters atom placed into a dictionary.
        
    list : resFragments
        The list of all found and parsed hydrated residue files with their
        protein and waters atom placed into a dictionary.
        
    solvate_globals.globalSettings : settings
        Instance of the settings class contaning all the options and values.

    Returns
    -------
    list : resMatches
        A list with entry for each residue, where the entry is a list of all hydrated fragments that
        were matched to the particular residue using the current settings (backbone, bestOnly,
        threshold, ...)

    """
    ### Log progress
    solvate_log.writeLog                              ( "Starting fragment matching", settings, 1 )
    
    ### Initialise variables
    resInd                                            = 0
    resMatches                                        = []
    
    ### For each amino acid
    for res in resList:

        ### How are we matching?
        if settings.useBackboneAtoms:

            ### Initialise variables
            resInd                                    = resInd + 1

            ### Report progress
            solvate_log.writeLog                      ( "Processing residue number " + str ( resInd ) + " (" + str ( res[0] ) + ")", settings, 2 )

            ### Match the residue to hydrated fragments using backbone step
            resMatches.append                         ( matchFragmentsUsingBackbone ( res, resFragments, fragFragments, settings ) )

        else:

            ### Initialise variables
            resInd                                    = resInd + 1

            ### Report progress
            solvate_log.writeLog                      ( "Processing residue number " + str ( resInd ) + " (" + str ( res[0] ) + ")", settings, 2 )

            ### Match the residue to hydrated fragments using backbone step
            resMatches.append                         ( matchFragmentsDirect ( res, fragFragments, settings ) )

    ### Log progress
    solvate_log.writeLog                              ( "Fragment matching complete", settings, 2 )

    ### Done
    return                                            ( resMatches )


######################################################
# findResFragSecondaryStructure()
def findResFragSecondaryStructure ( res, resFragments, settings ):
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
        
    solvate_globals.globalSettings : settings
        Instance of the settings class contaning all the options and values.

    Returns
    -------
    list : secStr
        List of the letters coding the secondary structure of the fitting fragments or empty list if no match
        was found. The number of returned values depends on whether the best fragment only or all passing
        fragments are to be returned.

    """
    ### Report progress
    solvate_log.writeLog                              ( "Started secondary structure search.", settings, 3 )
    
    ### For each fragment, find distance to this residue
    distances                                         = []
    fragNames                                         = []
    for frag in resFragments.keys():
        
        # Ignore entries with different amino acid
        if frag.split( os.path.sep )[ len ( frag.split( os.path.sep ) ) - 1 ].split( "_" )[0] != res[0]:
            continue
        
        # Compute distance between hydrated residue fragment and the processed residue
        distances.append                              ( solvate_maths.residueToResFragmentDistance_backbone ( res, resFragments[frag]["protein"], settings ) )
        fragNames.append                              ( frag )
        
    ### Are we searching for best, or all passing?
    if settings.bestFragmentOnly:
        
        ### Find best hydrated residue fragment in terms of backbone, if it exists
        minimalIndex                                  = solvate_maths.findSmallestDistance ( distances, settings.RMSDthreshold )
        
        ### If no distance is below the threshold
        if minimalIndex == -1:
            solvate_log.writeLog                      ( "Failed to find any matching hydrated residue fragment.", settings, 3 )
            return                                    ( [] )
        
        ### Determine best hydrated residue fragment's secondary structure from name
        secStr                                        = fragNames[minimalIndex].split( os.path.sep )[len ( fragNames[minimalIndex].split( os.path.sep ) ) - 1].split ( "_" )[1]
        
        ### Report progress
        solvate_log.writeLog                          ( "Best backbone fitting residue fragment is " + fragNames[minimalIndex].split( os.path.sep )[len (    fragNames[minimalIndex].split( os.path.sep ) ) - 1] + " with RMSD of " + str( distances[minimalIndex][0] ) + " and secondary structure therefore is " + str(    secStr ), settings, 4 )
    
        ### Done
        return                                        ( [ secStr ] )
        
    else:
        
        ### Initialise variables
        secStr                                        = []
        
        ### Find passing indices
        indices                                       = solvate_maths.findPassingDistances ( distances, settings.RMSDthreshold )
        
        ### If no distance is below the threshold
        if indices == -1:
            solvate_log.writeLog                      ( "Failed to find any matching hydrated residue fragment.", settings, 4 )
            return                                    ( [] )
        
        ### For each index
        for ind in indices:
            
            # Determine the secondary structure for the fragment
            secStr.append                             ( fragNames[ind].split( os.path.sep )[len ( fragNames[ind].split( os.path.sep ) ) - 1].split ( "_" )[1] )
            
            # Report progress
            solvate_log.writeLog                      ( "A backbone fitting residue fragment " + fragNames[ind].split( os.path.sep )[len (    fragNames[ind].split( os.path.sep ) ) - 1] + " with RMSD of " + str( distances[ind][0] ) + " and secondary structure " + str( fragNames[ind].split( os.path.sep )[len ( fragNames[ind].split( os.path.sep ) ) - 1].split ( "_" )[1] ) + " found.", settings, 4 )
        
        ### Done ( returning unique results )
        return                                        ( list ( set ( secStr ) ) )

######################################################
# findResFragRotamer ()
def findResFragRotamer ( res, resFragments, secStr, settings ):
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
        
    solvate_globals.globalSettings : settings
        Instance of the settings class contaning all the options and values.

    Returns
    -------
    dictionary : { rotamer : tForm }
        Key: List of the letters coding the rotamers of the fitting fragments followed by underscore and the
        appropriate secondary structure letter.
        Value: The dictionary of procrustes transformations.
        Alternatively, an empty dict if no match was found. For ALA
        and GLY, "NOT_INCLUDED" will be returned instead.

    """
    ### Report progress
    solvate_log.writeLog                              ( "Started rotamer search.", settings, 3 )
    
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
            distances.append                          ( solvate_maths.residueToResFragmentDistance_rotamer ( res, resFragments[frag]["protein"], settings ) )
            forNames.append                           ( frag )
            
        ### Are we searching for best, or all passing?
        if settings.bestFragmentOnly:
            
            # Find best hydrated residue fragment in terms of rotamer, if it exists
            minimalIndex                              = solvate_maths.findSmallestDistance ( distances, settings.RMSDthreshold )
    
            # If no distance is below the threshold
            if minimalIndex == -1:
                solvate_log.writeLog                  ( "Failed to find any matching hydrated residue fragment.", settings, 4 )
                return                                ( [] )
            
            # Find best rotamer
            rotamer                                   = forNames[minimalIndex].split( os.path.sep )[ len( forNames[minimalIndex].split( os.path.sep ) ) - 1 ].split( "."    )[0].split( "_" )[2] + "_" + forNames[minimalIndex].split( os.path.sep )[ len( forNames[minimalIndex].split( os.path.sep ) ) - 1 ].split("_")[1];
    
            ### Report progress
            solvate_log.writeLog                      ( "Best rotamer fitting residue fragment is " + forNames[minimalIndex].split( os.path.sep )[len ( forNames[minimalIndex].split( os.path.sep ) ) - 1] + " with RMSD of " + str( distances[minimalIndex][0] ) + " and the rotamer therefore is " + str( rotamer ), settings, 4 )
            
            ### Done
            return                                    ( [ { rotamer : distances[minimalIndex][1] } ] )
            
        else:
        
            # Initialise variables
            rotamer                                   = []
            
            # Find passing indices
            indices                                   = solvate_maths.findPassingDistances ( distances, settings.RMSDthreshold )
            
            # If no distance is below the threshold
            if indices == -1:
                solvate_log.writeLog                  ( "Failed to find any matching hydrated residue fragment.", settings, 4 )
                return                                ( [] )
            
            ### For each index
            for ind in indices:
                
                # Determine the secondary structure for the fragment
                rotamer.append                        ( { forNames[ind].split( os.path.sep )[ len( forNames[ind].split( os.path.sep ) ) - 1 ].split( "."    )[0].split( "_" )[2] + "_" + forNames[ind].split( os.path.sep )[ len( forNames[ind].split( os.path.sep ) ) - 1 ].split("_")[1] : distances[ind][1] } )
                
                # Report progress
                solvate_log.writeLog                  ( "A rotamer fitting residue fragment " + forNames[ind].split( os.path.sep )[len ( forNames[ind].split( os.path.sep ) ) - 1] + " with RMSD of " + str( distances[ind][0] ) + " and the rotamer " + str( forNames[ind].split( os.path.sep )[ len( forNames[ind].split( os.path.sep ) ) - 1 ].split( "."    )[0].split( "_" )[2] ) + " found.", settings, 4 )

            ### Done
            return                                    ( rotamer )

    else:

        ### Report progress
        solvate_log.writeLog                          ( "This is either ALA or GLY residue. There are no rotamers for these residues...", settings, 4 )
        
        ### Done
        return                                        ( "NOT_INCLUDED" )
    
    ### Done
    return                                            ( rotamer )

######################################################
# matchFragmentsUsingBackbone ()
def matchFragmentsUsingBackbone ( res, resFragments, fragFragments, settings ):
    """
    This function first matches all appropriate hydrated residues backbone atoms to
    the supplied residue backbone atoms. From this, it selects either the best hydrated
    residue or all passing hydrated residues and obtains their secondary structure assignment.
    
    It then takes all appropriate hydrated residues with the correct secondary structure as
    already determined and aligns their sidechains against the supplied residue, again either
    taking the best alignment or all passing. From these, it determines the best rotamer(s).
    
    Now, using the secondary structure and rotamer information, it aligns the appropriate
    hydrated fragment to the supplied residue and it passes the threshold, it will return
    a list of all passing fragments and the procrustes transformations required to reach the
    optimal overlay.

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
        
    solvate_globals.globalSettings : settings
        Instance of the settings class contaning all the options and values.

    Returns
    -------
    list of dictionaries : matched
        A list of all fragments that were matched to the input residue using the global criteria ( RMSD threshold,
        onlyBest, ... ). Each list entry is a dictionary with keys being the matched filenames and values being the
        the dictionary of procrustes transformations.

    """
    ### Find the secondary structure
    secStr                                            = findResFragSecondaryStructure ( res, resFragments, settings )
    
    ### Find the rotamer for all passing secondary structures
    rotamer                                           = []
    fragList                                          = []
    for ss in secStr:
    
        # Search for rotamers using this secondary structure
        rotamerHlp                                    = findResFragRotamer ( res, resFragments, secStr, settings )

        # Check for ALA and GLY
        if rotamerHlp == "NOT_INCLUDED":
        
            # Initialise variables
            dists                                     = []
            forNames                                  = []
            
            # Report progress
            for fl in fragFragments.keys():
            
                # Ignore fragments with wrong amino acid type
                if fl.split( os.path.sep )[ len ( fl.split( os.path.sep ) ) - 1 ].split("_")[0] != res[0]:
                    continue
                    
                # Ignore fragments with wrong secondary structure
                if fl.split( os.path.sep )[ len ( fl.split( os.path.sep ) ) - 1 ].split("_")[1] != ss:
                    continue
            
                # Compute the procrustes
                dists.append                          ( solvate_maths.residueToResFragmentDistance_direct ( res, fragFragments[fl]["protein"], settings ) )
                forNames.append                       ( fl )
                
            # Deal with dists
            if settings.bestFragmentOnly:
                
                # Find best hydrated residue fragment in terms of rotamer, if it exists
                minimalIndex                          = solvate_maths.findSmallestDistance ( dists, settings.RMSDthreshold )
            
                # If no distance is below the threshold
                if minimalIndex == -1:
                    solvate_log.writeLog              ( "Failed to find any matching hydrated residue fragment.", settings, 4 )
                    return                            ( [] )
                
                # Report progress
                solvate_log.writeLog                  ( "Best fragment fitting the residue is " + forNames[minimalIndex].split( os.path.sep )[len ( forNames[minimalIndex].split( os.path.sep ) ) - 1] + " with RMSD of " + str( dists[minimalIndex][0] ) + ".", settings, 4 )
                
                ### Done
                fragList.append                       ( { forNames[minimalIndex]:  dists[minimalIndex][1] } )
                
            else:
            
                # Initialise variables
                mtchs                                 = []
                
                # Find passing indices
                indices                               = solvate_maths.findPassingDistances ( dists, settings.RMSDthreshold )
                
                # If no distance is below the threshold
                if indices == -1:
                    solvate_log.writeLog              ( "Failed to find any matching hydrated residue fragment.", settings, 4 )
                    return                            ( [] )
                
                ### For each index
                for ind in indices:
                    
                    # Determine the secondary structure for the fragment
                    mtchs.append                      ( { forNames[ind]: dists[ind][1] } )
                    
                    # Report progress
                    solvate_log.writeLog              ( "A fragment " + forNames[ind].split( os.path.sep )[len ( forNames[ind].split( os.path.sep ) ) - 1] + " fits the residue with RMSD of " + str( dists[ind][0] ) + " was found.", settings, 4 )

                ### Done ( returning unique results )
                fragList.append                       ( { forNames[ind]:  dists[ind][1] } )
            
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
            if fl.split( os.path.sep )[ len ( fl.split( os.path.sep ) ) - 1 ].split("_")[1] != str( list ( rot )[0] ).split( "_" )[1]:
                continue
                
            # Ignore fragments with wrong rotamer
            if fl.split( os.path.sep )[ len ( fl.split( os.path.sep ) ) - 1 ].split("_")[2] != str( list ( rot )[0] ).split( "_" )[0]:
                continue

            # Passed! Now write up and save
            solvate_log.writeLog                      ( "Matched residue to fragment " + str ( fl ) + ".", settings, 3 )
            fragList.append                           ( { fl : rot[list( rot )[0]] } )
            
    ### Save to output
    return                                            ( fragList )

######################################################
# matchFragmentsDirect ()
def matchFragmentsDirect ( res, fragFragments, settings ):
    """
    This function computes the distances between all appropriate hydrated fragments and the
    supplied residue and then it selects either the best matching hydrated fragment, or all
    hydrated fragments passing under the threshold. It then returns the matched hydrated fragments
    as well as the operations required for their procrustes overlay with the residue.

    Parameters
    ----------
    list : res
        A list of all atoms parsed for this residue.
        
    dictionary : fragFragments
        A dictionary holding the parsed data for each hydrated fragment in terms of both
        protein and waters atoms.
        
    solvate_globals.globalSettings : settings
        Instance of the settings class contaning all the options and values.

    Returns
    -------
    list of dictionaries : matched
        A list of all fragments that were matched to the input residue using the global criteria ( RMSD threshold,
        onlyBest, ... ). Each list entry is a dictionary with keys being the matched filenames and values being the
        the dictionary of procrustes transformations.

    """
    ### Initialise variables
    distances                                         = []
    forNames                                          = []
    
    ### For each fragment, get the distances
    flIter                                            = 0
    for fl in fragFragments.keys():
    
        # Ignore fragments with wrong amino acid type
        if fl.split( os.path.sep )[ len ( fl.split( os.path.sep ) ) - 1 ].split("_")[0] != res[0]:
            continue
        
        distances.append                              ( solvate_maths.residueToResFragmentDistance_direct ( res, fragFragments[fl]["protein"], settings ) )
        forNames.append                               ( fl )
        flIter                                        = flIter + 1
        
    ### Report progress
    solvate_log.writeLog                              ( "Found " + str( flIter ) + " hydrated fragments for this residue.", settings, 3 )
        
    ### Are we searching for best, or all passing?
    if settings.bestFragmentOnly:
        
        # Find best hydrated residue fragment in terms of rotamer, if it exists
        minimalIndex                                  = solvate_maths.findSmallestDistance ( distances, settings.RMSDthreshold )
    
        # If no distance is below the threshold
        if minimalIndex == -1:
            solvate_log.writeLog                      ( "Failed to find any matching hydrated residue fragment.", settings, 4 )
            return                                    ( [] )
        
        # Report progress
        solvate_log.writeLog                          ( "Best fragment fitting the residue is " + forNames[minimalIndex].split( os.path.sep )[len ( forNames[minimalIndex].split( os.path.sep ) ) - 1] + " with RMSD of " + str( distances[minimalIndex][0] ) + ".", settings, 4 )
        
        ### Done
        return                                        ( [ { forNames[minimalIndex]:  distances[minimalIndex][1] } ] )
        
    else:
    
        # Initialise variables
        mtchs                                         = []
        
        # Find passing indices
        indices                                       = solvate_maths.findPassingDistances ( distances, settings.RMSDthreshold )
        
        # If no distance is below the threshold
        if indices == -1:
            solvate_log.writeLog                      ( "Failed to find any matching hydrated residue fragment.", settings, 4 )
            return                                    ( [] )
        
        ### For each index
        for ind in indices:
            
            # Determine the secondary structure for the fragment
            mtchs.append                              ( { forNames[ind]: distances[ind][1] } )
            
            # Report progress
            solvate_log.writeLog                      ( "A fragment " + forNames[ind].split( os.path.sep )[len ( forNames[ind].split( os.path.sep ) ) - 1] + " fits the residue with RMSD of " + str( distances[ind][0] ) + " was found.", settings, 4 )

        ### Done ( returning unique results )
        return                                        ( mtchs )
