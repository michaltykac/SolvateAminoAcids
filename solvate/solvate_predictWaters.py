# -*- coding: UTF-8 -*-
#   \file solvate_predictWaters.py
#   \brief This file provides water molecules predicting algorithms.
#
#   This file contains function for achieving several modifications. There are functions for
#   obtaining all waters from the hydrated fragments matched to individual input resudues, including
#   the correct rotation and translation of these to to input residue space.
#
#   Moreover, there are functions allowing detection of clashes between water molecules predicted as
#   described above and various "versions" of the input structure. These "versions" include the original
#   input structure, the structure with no hydrogens, with no already determined water molecules or with
#   no ligands.
#
#   Finally, there are also functions required to take the resulting water molecules list with no clashes
#   and clustering these in order to remove (or rather combine) very similar water molecules into a single
#   most likely water molecule position.
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

import numpy                                          ### Maths
import solvate.solvate_log as solvate_log             ### Writing log
import solvate.solvate_maths as solvate_maths         ### Computing distances and rotations
import solvate.solvate_structures as solvate_structures ### Fragment match file creation
import gemmi                                          ### File reading and writing
import os                                             ### Path separator
from sklearn.cluster import DBSCAN                    ### DBSCAN alorithm

######################################################
# predictWaters ()
def predictWaters ( resList, matchedFrags, fragFragments, settings ):
    """
    This function takes the parsed residues, the list of all matched fragments to each residue and
    the list of all parsed fragments and proceeds to rotate and translate all the water molecules
    so that they would have the same position relative to the input structure residue as the matched
    fragment, effectively adding the waters from hydrated fragments to the residues.
    
    If the matched fragments folder path is set, then this function will also write a PDB file for
    each residue-fragment match to allow re-viewing how good the matches are.
    
    Warning: There is no checking for the water molecules to the co-ordinates. It is expected that
    this will be done later, beforer wirting these water into anywhere.

    Parameters
    ----------
    list : resList
        The list of all residues as returned by parseOutCoords() function.
        
    list : matchedFrags
        A list with entry for each residue, where the entry is a list of all hydrated fragments that
        were matched to the particular residue.
        
    list : fragFragments
        The list of all found and parsed hydrated fragment files with their
        protein and waters atom placed into a dictionary.
        
    solvate_globals.globalSettings : settings
        Instance of the settings class contaning all the options and values.

    Returns
    -------
    list : waters
        List of all water molecules obtained from the assigned hydrated matched fragments for each
        residue. These water molecules are rotated and translated to be in appropriate positions
        to the residue position.

    """
    ### Log progress
    solvate_log.writeLog                              ( "Starting water prediction", settings, 1 )
    
    ### Initialise variables
    predictedWaters                                   = []

    ### For each residue
    resCounter                                        = 1
    for res in range ( 0, len ( resList ) ):
    
        ### Report progress
        solvate_log.writeLog                          ( "Starting water prediction for residue " + str( resCounter ) + " (" + str( resList[res][0] ) + ")", settings, 2 )
        
        ### Initialise variables
        predictedWatersForRes                         = []
        
        ### For each matched fragment
        for frag in range ( 0, len ( matchedFrags[res] ) ):
        
            ### For each key in dict (should be only one)
            for frName in matchedFrags[res][frag].keys():
            
                ### Check for water atoms
                if len( fragFragments[frName]["waters"] ) == 0:
                
                    ### Report progress
                    solvate_log.writeLog              ( "No waters found for fragment " + str( frName ), settings, 3 )
                    
                ### Report progress
                solvate_log.writeLog                  ( "Rotating and translating waters from fragment " + str( frName ), settings, 3 )
                
                ### If required, open the fragment file and add the residue as another chain
                resModel                              = solvate_structures.prepareMatchFile ( resList[res], fragFragments[frName], matchedFrags[res][frag][frName], settings )
                
                ### For each water atom
                matchedWaters                         = []
                for wIt in range ( 0, len ( fragFragments[frName]["waters"] ) ):
                
                    ### Rotate and translate each water atom
                    waterAtomPos                      = numpy.array ( [ fragFragments[frName]["waters"][wIt][1],
                                                                        fragFragments[frName]["waters"][wIt][2],
                                                                        fragFragments[frName]["waters"][wIt][3] ] )
                    waterAtomPos                      = solvate_maths.rotateAndTranslateAtom ( waterAtomPos, matchedFrags[res][frag][frName] )
                    
                    ### Save also the name, occupancy and B-factor
                    waterAtomPos                      = numpy.append ( waterAtomPos, numpy.array ( [ fragFragments[frName]["waters"][wIt][0],
                                                                                                     fragFragments[frName]["waters"][wIt][4],
                                                                                                     fragFragments[frName]["waters"][wIt][5] ] ) )

                    ### Save the predicted water atom
                    predictedWatersForRes.append      ( waterAtomPos )
                    matchedWaters.append              ( waterAtomPos )
                    
                
                ### Add waters and fragment to structure and write it, if required
                matchName                             = str ( resCounter ) + "_" + str ( resList[res][0] ) + "-" + str ( frName.split( os.path.sep )[ len ( frName.split( os.path.sep ) ) - 1 ].split(".")[0] ) + ".pdb"
                solvate_structures.completeAndWriteMatchFile ( resModel, matchedWaters, matchName, settings )
                    
                
        ### Save the rotated and translated waters for this match
        predictedWaters.append                        ( predictedWatersForRes )
        
        ### Prepare for next residue
        resCounter                                    = resCounter + 1

    ### Log progress
    solvate_log.writeLog                              ( "Water prediction complete", settings, 2 )

    ### Done
    return                                            ( predictedWaters )

######################################################
# removeClashes ()
def removeClashes ( waters, settings ):
    """
    This function takes the list of all predicted water molecules and attempts to detect any clashes
    between any such molecule and the original input structure, or any "version" of it. Here, by "version"
    we mean the input structure with no hydorgens, no already detected waters or no ligands. It then returns
    a list of water molecules with no clashes for each of the required "version".

    Parameters
    ----------
    list : waters
        List of all water molecules obtained from the assigned hydrated matched fragments for each
        residue. These water molecules are rotated and translated to be in appropriate positions
        to the residue position.
        
    solvate_globals.globalSettings : settings
        Instance of the settings class contaning all the options and values.

    Returns
    -------
    list : noClashWaters
        List containing a list of water molecules which have no clashes (as defined by the settings)
        to the various versions of the input structure. The versions are: full structure, no hydrogens,
        no waters, no ligands. If some of the versions were not requested, the list for that version will
        be empty.
        
    list : waterLabels
        List of labels describing each of the noClashWaters list entries (i.e. what was used and not used
        for clash detection)

    """
    ### Log progress
    solvate_log.writeLog                              ( "Removing clashes between new waters and structure contents", settings, 1 )
    
    ### Initialise variables
    noClashWaters                                     = []

    ### Should full structure contents be used to detect clashes?
    if not settings.noFullStructure:
        
        ### Get list of waters not clashing wiht the full structure
        noClashWaters.append                          ( removeClashingWaters ( waters, settings.inputCoordinateStructure, "full structure", settings ) )
    
    else:
    
        ### Return empty list if full structure no clashes are not required
        noClashWaters.append                          ( [] )
        
    ### Should hydrogens be removed for clashes detection?
    if settings.noHydro:
        
        ### Produce structure with no hydrogens
        settings.inputCoordinateStructureNoHydro      = gemmi.read_structure ( settings.inputCoordinateFile )
        for model in settings.inputCoordinateStructureNoHydro:
            for chain in model:
                for res in chain:
                    delAtInds                         = []
                    for atIndex, at in enumerate ( res ):
                        if at.is_hydrogen():
                            delAtInds.append          ( atIndex )
                    for atIt in sorted ( delAtInds, reverse = True ):
                        del res[atIt]
        
        ### Get list of waters, ignoring clashes with hydrogens
        noClashWaters.append                          ( removeClashingWaters ( waters, settings.inputCoordinateStructureNoHydro, "ignore hydrogens", settings ) )
    
    else:
    
        ### Return empty list if hydrogen clashes ignoring is not required
        noClashWaters.append                          ( [] )
        
    ### Should already existing waters be removed for clashes detection?
    if settings.noWaters:
        
        ### Produce structure with no waters
        settings.inputCoordinateStructureNoWaters     = gemmi.read_structure ( settings.inputCoordinateFile )
        for model in settings.inputCoordinateStructureNoWaters:
            for chain in model:
                delResInds                            = []
                for resInd, res in enumerate ( chain ):
                    if res.is_water():
                        delResInds.append             ( resInd )
                for rIt in sorted ( delResInds, reverse = True ):
                    del chain[rIt]
        
        ### Get list of waters, ignoring clashes with hydrogens
        noClashWaters.append                          ( removeClashingWaters ( waters, settings.inputCoordinateStructureNoWaters, "ignore existing waters", settings ) )
    
    else:
    
        ### Return empty list if hydrogen clashes ignoring is not required
        noClashWaters.append                          ( [] )
        
    ### Should already existing waters be removed for clashes detection?
    if settings.noLigand:
        
        ### Produce structure with no ligands
        settings.inputCoordinateStructureNoLigand     = gemmi.read_structure ( settings.inputCoordinateFile )
        for model in settings.inputCoordinateStructureNoLigand:
            for chain in model:
                delResInds                            = []
                for resInd, res in enumerate ( chain ):
                    if res in chain.get_ligands():
                        delResInds.append             ( resInd )
                for rIt in sorted ( delResInds, reverse = True ):
                    del chain[rIt]

        ### Get list of waters, ignoring clashes with hydrogens
        noClashWaters.append                          ( removeClashingWaters ( waters, settings.inputCoordinateStructureNoLigand, "ignore ligands", settings ) )
    
    else:
    
        ### Return empty list if hydrogen clashes ignoring is not required
        noClashWaters.append                          ( [] )
        
    ### Log progress
    solvate_log.writeLog                              ( "Clashes removed for all requested structures", settings, 2 )
        
    ### Done
    waterLabels                                       = [ "full structure", "ignore hydrogens", "ignore existing waters", "ignore ligands" ]
    return                                            ( noClashWaters, waterLabels )

######################################################
# removeClashingWaters ()
def removeClashingWaters ( waters, structure, structureState, settings ):
    """
    This function takes a structure, which may have different parts removed in order to make clash
    detection not use such removed parts. It then searches each water from any matched hydrated fragment
    for having a neighbouring atom within the threshold radius and only water molecules with no neighbour in
    such threshold radius are kept.
    
    It then proceeds to return a list containing a list of all passing water molecules for each residue.

    Parameters
    ----------
    list : waters
        List of all water molecules obtained from the assigned hydrated matched fragments for each
        residue. These water molecules are rotated and translated to be in appropriate positions
        to the residue position.
        
    gemmi.Structure : structure
        The structure all of which contents will be used to check for clashes.
        
    str : structureState
        String describing the type of contents that were removed from the structure for
        purposes of clashes detection.
        
    solvate_globals.globalSettings : settings
        Instance of the settings class contaning all the options and values.

    Returns
    -------
    list : noClashWaters
        List of lists, where each residue has its own list of atoms passing the clash detection.

    """
    ### Log progress
    solvate_log.writeLog                              ( "Removing clashes between new waters and structure " + str( structureState ), settings, 2 )
    
    ### Initialise variables
    neighSearches                                     = []
    noClashWaters                                     = []
    
    ### For each model
    for model in structure:
    
        ### Prepare the model for the cell-lists neighbour-search
        neighSearches.append                          ( gemmi.NeighborSearch ( model, structure.cell, settings.clashWithinRadius ).populate() )
        
    ### For each residue
    resCtr                                            = 0
    for res in waters:
    
        ### Update variables
        resCtr                                        = resCtr + 1
        resNCWaters                                   = []
        
        ### Log progress
        solvate_log.writeLog                          ( "Working on residue " + str ( resCtr ), settings, 3 )
    
        ### For each water molecule
        for watM in res:
        
            ### Initialise variables
            noNeighbours                              = True
        
            ### For each model
            for neighSearch in neighSearches:
        
                ### Search for neighbours
                neighbours                            = neighSearch.find_atoms ( gemmi.Position ( float ( watM[0] ),
                                                                                                  float ( watM[1] ),
                                                                                                  float ( watM[2] ) ),
                                                                                 '\0',
                                                                                 settings.clashWithinRadius )
            
                ### Check for number of neighbours
                if len ( neighbours ) > 0:
                    noNeighbours                      = False
                
            ### If no model has any neighbours
            if noNeighbours:
                resNCWaters.append                    ( watM )
            
        ### Add all waters for this residue to return list
        noClashWaters.append                          ( resNCWaters )
        
        ### Log progress
        solvate_log.writeLog                          ( "Found " + str ( len ( resNCWaters ) ) + " water molecules with no neighbours.", settings, 4 )
        
    ### Done
    return                                            ( noClashWaters )

######################################################
# clusterWaters ()
def clusterWaters ( noClashWaters, waterLabels, settings ):
    """
    This function takes the list of predicted water molecules which do not clash with any protein atoms
    and it proceeds to cluster these using the DBSCAN algorithm. It then returns the labels for each
    water molecule as well as the total number of clusters for each member of the input list.

    Parameters
    ----------
    list : noClashWaters
        List containing a list of water molecules which have no clashes (as defined by the settings)
        to the various versions of the input structure. The versions are: full structure, no hydrogens,
        no waters, no ligands. If some of the versions were not requested, the list for that version will
        be empty.
        
    list : waterLabels
        List of labels describing each of the noClashWaters list entries (i.e. what was used and not used
        for clash detection)
        
    solvate_globals.globalSettings : settings
        Instance of the settings class contaning all the options and values.

    Returns
    -------
    list : clusteredWaters
        A list of cluster labels, cluster count, water positions and molecule occupancies and B factors for
        each structure type available in the noClashWaters input variable.

    """
    ### Log progress
    solvate_log.writeLog                              ( "Clustering predicted non-clashing water molecules", settings, 1 )
    
    ### Initialise variables
    clusteredWaters                                   = []
    
    ### For each structure
    stIt                                              = 0
    for strType in noClashWaters:
        
        ### Do nothing if not required
        if len ( strType ) == 0:
            continue
        
        ### Log progress
        solvate_log.writeLog                          ( "Started clustering predicted non-clashing water molecules for waters from " + str( waterLabels[stIt] ), settings, 2 )
        
        ### Find how many waters we have
        noWaters                                      = 0
        for res in strType:
            for at in res:
                noWaters                              = noWaters + 1
        
        ### Create numpy array
        watPositions                                  = numpy.zeros ( ( noWaters, 3 ), dtype=float)
        watOccs                                       = numpy.zeros ( ( noWaters, 1 ), dtype=float )
        watBFacts                                     = numpy.zeros ( ( noWaters, 1 ), dtype=float )
        
        ### For each water position
        noWaters                                      = 0
        for res in strType:
            for at in res:
            
                ### Save to new array
                watPositions[noWaters][0]             = at[0]
                watPositions[noWaters][1]             = at[1]
                watPositions[noWaters][2]             = at[2]
                
                ### Save occupancies and B factors as well
                watOccs[noWaters]                     = at[4]
                watBFacts[noWaters]                   = at[5]
                
                ### Prepare for next iteration
                noWaters                              = noWaters + 1
        
        ### Cluster using DBSCAN
        clusterObj                                    = DBSCAN ( eps = settings.maxClustDist, min_samples = 1 ).fit ( watPositions )
        
        ### Get clusters
        labels                                        = clusterObj.labels_
        nClusters                                     = len ( set ( labels ) ) - ( 1 if -1 in labels else 0 )
        
        ### Log progress
        solvate_log.writeLog                          ( "Found " + str ( noWaters ) + " waters for structure type " + str( waterLabels[stIt] ) + ", which clustered to " + str ( nClusters ) + " clusters.", settings, 3 )
        stIt                                          = stIt + 1
        
        ### Save results
        clusteredWaters.append                        ( [ labels, nClusters, watPositions, watOccs, watBFacts ] )
    
    ### Log progress
    solvate_log.writeLog                              ( "Completed water molecules clustering", settings, 2 )
    
    ### Done
    return                                            ( clusteredWaters )
