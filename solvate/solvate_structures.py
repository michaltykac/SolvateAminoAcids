# -*- coding: UTF-8 -*-
#   \file solvate_structures.py
#   \brief This file contains functions for reading and writing structure data.
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

import sys                                            ### In the case complete stop is needed
import os                                             ### Checking for file existence
import glob                                           ### Listing all files in folder
import gemmi                                          ### Library of choice for structure reading/writing
import numpy                                          ### Get numpy.array's

import solvate.solvate_log as solvate_log             ### For writing the log
import solvate.solvate_maths as solvate_maths         ### Atom rotation and translation

######################################################
# parseInputCoordinates ()
def parseInputCoordinates ( settings ):
    """
    This function combines the readInCoordinates and parseOutCoords functions for user
    convenience, taking the settings object and returning the resList required for
    fragment matching later.

    Parameters
    ----------
    solvate_globals.globalSettings : settings
        Instance of the settings class contaning all the options and values.

    Returns
    -------
    list : retList
        This is a list containing all the information solvate requires from the structure.

    """
    ### Read in the structure
    inputCoords                                       = readInCoordinates ( settings )
    
    ### Save for later
    settings.setInputStructure                        ( inputCoords )
    
    ### Get residue list
    resList                                           = parseOutCoords ( inputCoords, settings )
    
    ### Done
    return                                            ( resList )

######################################################
# readInCoordinates ()
def readInCoordinates ( settings ):
    """
    This function reads in a co-ordinates file and returns the gemmi class representing them.

    Parameters
    ----------
    solvate_globals.globalSettings : settings
        Instance of the settings class contaning all the options and values.

    Returns
    -------
    gemmi.Structure: structure
        This is the gemmi class representing the parsed structure.

    """
    ### Check if path exists
    if not os.path.isfile ( settings.inputCoordinateFile ):
        sys.exit ( 'Error: The supplied file ' + str ( settings.inputCoordinateFile ) + ' does not seem to exist. Please supply a different file.' )
        
    ### Log progress
    solvate_log.writeLog                              ( "Starting to read structure " + str( settings.inputCoordinateFile ), settings, 1 )
    
    ### Read in the structure
    structure                                         = gemmi.read_structure ( settings.inputCoordinateFile )
    
    ### Log progress
    solvate_log.writeLog                              ( "Structure parsed", settings, 2 )
    
    ### Done
    return                                            ( structure )
    
######################################################
# parseOutCoords ()
def parseOutCoords ( structure, settings ):
    """
    This function takes a structure and parses out the co-ordinates in a way that is simpler and faster to
    iterate through than using the full structure. This is also where alternative locations are dealt with.

    Parameters
    ----------
    gemmi.Structure: structure
        This is the gemmi class representing an already parsed co-ordinates.
        
    solvate_globals.globalSettings : settings
        Instance of the settings class contaning all the options and values.

    Returns
    -------
    list : retList
        This is a list containing all the information solvate requires from the structure.

    """
    ### Log progress
    solvate_log.writeLog                              ( "Starting to parse structure for speed-up", settings, 1 )
    
    ### Local variables
    noAtm                                             = 0
    noRes                                             = 0
    noCha                                             = 0
    noMod                                             = 0
    resList                                           = []
    
    ### For each model
    for model in structure:
    
        ### For each subchain
        for subchain in model:

            ### For each residue
            for res in subchain.first_conformer():
            
                ### Get information
                name  = res.name
                atoms = []
                
                ### Ignore anything but amino acids
                if name not in settings.aaTypes:
                    continue
            
                ### For each atom
                for atom in res:
                    atoms.append ( [ atom.name, atom.pos[0], atom.pos[1], atom.pos[2], atom.occ, atom.b_iso ] )
                    noAtm                             = noAtm + 1
                
                ### Save
                resList.append ( [ name, atoms ] )
                noRes                                 = noRes + 1
                
            noCha                                     = noCha + 1
            
        noMod                                         = noMod + 1

    ### Log progress
    solvate_log.writeLog                              ( "Found " + str ( noMod ) + " models, " + str ( noCha ) + " chains with total of " + str ( noRes ) + " residues containing total of " + str ( noAtm ) + " atoms.", settings, 3 )

    ### Log progress
    solvate_log.writeLog                              ( "Structure parsed", settings, 2 )
            
    ### Done
    return                                            ( resList )

######################################################
# readFragment ()
def readFragment ( fragPath, settings ):
    """
    This function reads a hydrated residue file and parses the required information from it.

    Parameters
    ----------
    str : fragPath
        This is the path to the fragment to be read.
        
    solvate_globals.globalSettings : settings
        Instance of the settings class contaning all the options and values.

    Returns
    -------
    list : protein
        This is a list containing all the information solvate requires about the protein atoms of the fragment.
        
    list : waters
        This is a list containing all the information solvate requires about the water atoms of the fragment.

    """
    ### Log progress
    solvate_log.writeLog                              ( "Starting to process fragment " + str ( fragPath ), settings, 4 )
    
    ### Read in the fragment
    fragStr                                           = gemmi.read_structure ( fragPath )
    
    ### Initialise variables
    waters                                            = []
    protein                                           = []
    
    ### For each chain in the fragment (assuming a single model here!)
    for subch in fragStr[0].subchains():
        
        ### For each residue in the chain
        for res in subch.first_conformer():
            
            ### For each atom in the residue
            for at in res:
                
                ### If the residue is water
                if res.is_water():
                    
                    ### Add atom to water's list
                    waters.append                     ( [ at.name, at.pos[0], at.pos[1], at.pos[2], at.occ, at.b_iso ] )
            
                ### For each residue w  hich is not water
                else:
                
                    ### Add atom to protein's list
                    protein.append                    ( [ at.name, at.pos[0], at.pos[1], at.pos[2], at.occ, at.b_iso ] )
              
    ### Report progres
    solvate_log.writeLog                              ( "Read in " + str ( len ( protein ) ) + " protein atoms and " + str ( len ( waters ) ) + " waters.", settings, 5 )
              
    ### Done
    return                                            ( protein, waters )

######################################################
# getAllResidueFragments ()
def getAllResidueFragments ( settings ):
    """
    This function reads in all hydrated residue fragments from the supplied folder and parses these into
    format readily accessible by solvate.

    Parameters
    ----------
    solvate_globals.globalSettings : settings
        Instance of the settings class contaning all the options and values.

    Returns
    -------
    dictionary : ret
        This dictionary contains all files as keys to sub-dictionaries, which have two keys each, waters and protein. These
        hold the atom data for the respective atoms.

    """
    ### Initialise variables
    ret                                               = {}
    
    ### Report progress
    solvate_log.writeLog                              ( "Starting residue fragment reading from path " + str( settings.resInputDir ) + ".", settings, 2 )
    
    ### Find all input hydrated residue files
    allFragFiles                                      = [ fr for fr in glob.glob ( os.path.join ( os.path.join ( settings.resInputDir, "*" ), "*.pdb" ), recursive = True ) ]
    
    ### If no frags, this is probably an error
    if len ( allFragFiles ) == 0:
        # Print error
        solvate_log.writeLog                          ( "!!! ERROR !!! Could not find any hydrated residues. Please use the --res option to supply the path to the hydrated residues folder location.", settings, 0 )
        
        # Terminate
        solvate_log.endLog                            ( settings )
        
    ### Read in all residue fragments
    for frag in allFragFiles:
            
        # Read it in
        ( protein, waters )                           = readFragment ( frag, settings )
        
        # Save it
        ret[frag]                                     = {}
        ret[frag]["protein"]                          = protein
        ret[frag]["waters"]                           = waters
        
    ### Report progress
    solvate_log.writeLog                              ( "All " + str( len ( allFragFiles ) ) + " residue fragments parsed.", settings, 3 )
    
    ### Done
    return                                            ( ret )

######################################################
# getAllFragmentFragments ()
def getAllFragmentFragments ( settings ):
    """
    This function reads in all hydrated fragments from the supplied folder and parses these into
    format readily accessible by solvate.

    Parameters
    ----------
    solvate_globals.globalSettings : settings
        Instance of the settings class contaning all the options and values.

    Returns
    -------
    dictionary : ret
        This dictionary contains all files as keys to sub-dictionaries, which have two keys each, waters and protein. These
        hold the atom data for the respective atoms.

    """
    ### Initialise variables
    ret                                               = {}
    
    ### Report progress
    solvate_log.writeLog                              ( "Starting fragment reading from path " + str( settings.fragInputDir ) + ".", settings, 2 )
    
    ### Find all input hydrated residue files
    allFragFiles                                      = [ fr for fr in glob.glob ( os.path.join ( os.path.join ( settings.fragInputDir, "*" ), "*.pdb" ), recursive = True ) ]
    
    ### If no frags, this is probably an error
    if len ( allFragFiles ) == 0:
        # Print error
        solvate_log.writeLog                          ( "!!! ERROR !!! Could not find any hydrated fragments. Please use the --frag option to supply the path to the hydrated fragments folder location.", settings, 0 )
        
        # Terminate
        solvate_log.endLog                            ( settings )
        
    ### Read in all residue fragments
    for frag in allFragFiles:
            
        # Read it in
        ( protein, waters )                           = readFragment ( frag, settings )
        
        # Save it
        ret[frag]                                     = {}
        ret[frag]["protein"]                          = protein
        ret[frag]["waters"]                           = waters
        
    ### Report progress
    solvate_log.writeLog                              ( "All " + str( len ( allFragFiles ) ) + " fragments parsed.", settings, 3 )
    
    ### Done
    return                                            ( ret )

######################################################
# prepareMatchFile ()
def prepareMatchFile ( residue, fragment, transform, settings ):
    """
    This function creates a gemmi.Model object and fills it with two chains. The first
    chain call "A" will contain all original residue atoms, while the second chain "B"
    will contain the hydrated fragment protein atoms rotated and translated to best match
    the residue atoms.
    
    If such model will not be writte, the function returns None instead to save time and
    memory.

    Parameters
    ----------
    list : residue
        A list of all the atoms found in this particular residue.
        
    dictionary : fragment
        A dictionary with two entries, "protein" containing all protein atoms and
        "waters" containing all water atoms.
        
    dictionary : transform
        The result of the procrustes analysis (t) containing the "rotation" entry
        with the optimal rotation matrix and "translation" entry with the optimal
        translation vector.
    
    solvate_globals.globalSettings : settings
        Instance of the settings class contaning all the options and values.

    Returns
    -------
    gemmi.Model : model
        Gemmi Model object, which now contains two chains, the "A" chain with the original
        residue atoms and the "B" model with the rotated and translated hydrated-fragment
        atoms.

    """
    ### If required
    if settings.matchedFragsPath != "":
        
        ### Initialise residue chain variables
        resModel                                      = gemmi.Model ( "1" )
        resChain                                      = gemmi.Chain ( "A" )
        resRes                                        = gemmi.Residue ( )
        
        ### For each residue atom
        for at in residue[1]:
        
            ### Set the atom info
            resAtom                                   = gemmi.Atom ( )
            resAtom.name                              = at[0]
            resAtom.pos                               = gemmi.Position ( at[1], at[2], at[3] )
            resAtom.occ                               = at[4]
            resAtom.b_iso                             = at[5]
            resRes.add_atom                           ( resAtom )
            
        ### Set the residue info and add to chain and model
        resRes.name                                   = residue[0]
        resRes.seqid                                  = gemmi.SeqId ( 1, ' ' )
        resRes.entity_type                            = gemmi.EntityType.Polymer
        resChain.add_residue                          ( resRes )
        resModel.add_chain                            ( resChain )
        
        ### Initialise fragment chain variables
        fragChain                                     = gemmi.Chain ( "B" )
        fragRes                                       = gemmi.Residue ( )
        
        ### For each fragment protein atom
        for at in fragment["protein"]:
        
            ### Rotate and translate each protein atom
            fragAtomPos                               = numpy.array ( [ at[1],
                                                                        at[2],
                                                                        at[3] ] )
            fragAtomPos                               = solvate_maths.rotateAndTranslateAtom ( fragAtomPos, transform )
                        
            ### Set the atom info
            fragAtom                                  = gemmi.Atom ( )
            fragAtom.name                             = at[0]
            fragAtom.pos                              = gemmi.Position ( fragAtomPos[0], fragAtomPos[1], fragAtomPos[2] )
            fragAtom.occ                              = at[4]
            fragAtom.b_iso                            = at[5]
            fragRes.add_atom                          ( fragAtom )
            
        ### Set the residue information and add to chain and model and structure
        fragRes.name                                  = residue[0]
        fragRes.seqid                                 = gemmi.SeqId ( 1, ' ' )
        fragRes.entity_type                           = gemmi.EntityType.Polymer
        fragChain.add_residue                         ( fragRes )
        resModel.add_chain                            ( fragChain )
        
        ### Done
        return                                        ( resModel )
        
    else:
        ### Writing not required. Just return None
        model                                         = None
        return                                        ( model )
        
######################################################
# completeAndWriteMatchFile ()
def completeAndWriteMatchFile ( resModel, matchedWaters, matchFileName, settings ):
    """
    This function takes the prepared model, adds the now rotated and translated water molecules
    and proceeds to write the file to the disk. If, however, this is not required, then it does
    nothing and just ends.

    Parameters
    ----------
    gemmi.Model : model
        Gemmi Model object, which contains two chains, the "A" chain with the original
        residue atoms and the "B" chain with the rotated and translated hydrated-fragment
        atoms.
    
    list : matchedWaters
        A list of all water atoms matched to this residue, with their positions rotated
        and translated to fir the residue already.
        
    str : matchFileName
        A filename to where the file should be save to.
        
    solvate_globals.globalSettings : settings
        Instance of the settings class contaning all the options and values.

    Returns
    -------
    NONE

    """
    ### If required

    if resModel is not None:
    
        ### Create water chain and add it
        watChain                                      = gemmi.Chain ( "W" )
                    
        ### For each water atom
        resCtr                                        = 1
        for at in matchedWaters:
                    
            ### Set the atom info
            watAtom                                   = gemmi.Atom ( )
            watAtom.name                              = str ( at[3] )
            watAtom.pos                               = gemmi.Position ( float ( at[0] ) , float ( at[1] ), float ( at[2] ) )
            watAtom.occ                               = float ( at[4] )
            watAtom.b_iso                             = float ( at[5] )
            watAtom.charge                            = 0
            watAtom.element                           = gemmi.Element ( "O" )
        
            ### Set the residue information and add to chain
            watRes                                    = gemmi.Residue ( )
            watRes.add_atom                           ( watAtom )
            watRes.name                               = "HOH"
            watRes.seqid                              = gemmi.SeqId ( resCtr, ' ' )
            watRes.entity_type                        = gemmi.EntityType.Water
            watChain.add_residue                      ( watRes )
            resCtr                                    = resCtr + 1
            
        ### Complete the chain and add to model and to structure
        resModel.add_chain                            ( watChain )
        matchStr                                      = gemmi.Structure ( )
        matchStr.add_model                            ( resModel )
        
        ### Write fragment-residue match
        matchStr.write_pdb                            ( os.path.join ( settings.matchedFragsPath, matchFileName ) )
        
        ### Report progress
        solvate_log.writeLog                          ( "Written matched hydrated fragment to residue alignement file with water to " + str( matchFileName ), settings, 4 )
        
        ### Done
        return
