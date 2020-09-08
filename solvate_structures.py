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
#   \version   0.0.1
#   \date      SEP 2020
######################################################

######################################################
# Imports

import sys                                            ### In the case complete stop is needed
import os                                             ### Checking for file existence
import glob                                           ### Listing all files in folder
import gemmi                                          ### Library of choice for structure reading/writing

import solvate_globals                                ### Global variables
import solvate_log                                    ### For writing the log

######################################################
# readInCoordinates ()
def readInCoordinates ( filePath ):
    """
    This function reads in a co-ordinates file and returns the gemmi class representing them.

    Parameters
    ----------
    str : filePath
        This is the path and filename to the input file that should be read in.

    Returns
    -------
    gemmi.Structure: structure
        This is the gemmi class representing the parsed structure.

    """
    ### Check if path exists
    if not os.path.isfile ( filePath ):
        sys.exit ( 'Error: The supplied file ' + str ( filePath ) + ' does not seem to exist. Please supply a different file.' )
        
    ### Log progress
    solvate_log.writeLog                              ( "Starting to read structure " + str( filePath ), 1 )
    
    ### Read in the structure
    structure                                         = gemmi.read_structure ( filePath )
    
    ### Log progress
    solvate_log.writeLog                              ( "Structure parsed", 2 )
    
    ### Done
    return                                            ( structure )
    
######################################################
# parseOutCoords ()
def parseOutCoords ( structure ):
    """
    This function takes a structure and parses out the co-ordinates in a way that is simpler and faster to
    iterate through than using the full structure. This is also where alternative locations are dealt with.

    Parameters
    ----------
    gemmi.Structure: structure
        This is the gemmi class representing an already parsed co-ordinates.

    Returns
    -------
    list : retList
        This is a list containing all the information solvate requires from the structure.

    """
    ### Log progress
    solvate_log.writeLog                              ( "Starting to parse structure for speed-up", 1 )
    
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
                if name not in solvate_globals._aaTypes:
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
    solvate_log.writeLog                              ( "Found " + str ( noMod ) + " models, " + str ( noCha ) + " chains with total of " + str ( noRes ) + " residues containing total of " + str ( noAtm ) + " atoms.", 3 )

    ### Log progress
    solvate_log.writeLog                              ( "Structure parsed", 2 )
            
    ### Done
    return                                            ( resList )

######################################################
# readFragment ()
def readFragment ( fragPath ):
    """
    This function reads a hydrated residue file and parses the required information from it.

    Parameters
    ----------
    str : fragPath
        This is the path to the fragment to be read.

    Returns
    -------
    list : protein
        This is a list containing all the information solvate requires about the protein atoms of the fragment.
        
    list : waters
        This is a list containing all the information solvate requires about the water atoms of the fragment.

    """
    ### Log progress
    solvate_log.writeLog                              ( "Starting to process fragment " + str ( fragPath ), 3 )
    
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
    solvate_log.writeLog                              ( "Read in " + str ( len ( protein ) ) + " protein atoms and " + str ( len ( waters ) ) + " waters.", 4 )
              
    ### Done
    return                                            ( protein, waters )

######################################################
# getAllResidueFragments ()
def getAllResidueFragments ( ):
    """
    This function reads in all hydrated residue fragments from the supplied folder and parses these into
    format readily accessible by solvate.

    Parameters
    ----------
    NONE

    Returns
    -------
    dictionary : ret
        This dictionary contains all files as keys to sub-dictionaries, which have two keys each, waters and protein. These
        hold the atom data for the respective atoms.

    """
    ### Initialise variables
    ret                                               = {}
    
    ### Report progress
    solvate_log.writeLog                              ( "Starting residue fragment reading from path " + str( solvate_globals._resInputDir ) + ".", 2 )
    
    ### Find all input hydrated residue files
    allFragFiles                                      = [ fr for fr in glob.glob ( os.path.join ( os.path.join ( solvate_globals._resInputDir, "*" ), "*.pdb" ), recursive = True ) ]
    
    ### If no frags, this is probably an error
    if len ( allFragFiles ) == 0:
        # Print error
        solvate_log.writeLog                          ( "!!! ERROR !!! Could not find any hydrated residues. Please use the --res option to supply the path to the hydrated residues folder location.", 0 )
        
        # Terminate
        solvate_log.endLog                            ( )
        
    ### Read in all residue fragments
    for frag in allFragFiles:
            
        # Read it in
        ( protein, waters )                           = readFragment ( frag )
        
        # Save it
        ret[frag]                                     = {}
        ret[frag]["protein"]                          = protein
        ret[frag]["waters"]                           = waters
        
    ### Report progress
    solvate_log.writeLog                              ( "Residue fragments parsed.", 3 )
    
    ### Done
    return                                            ( ret )

######################################################
# getAllFragmentFragments ()
def getAllFragmentFragments ( ):
    """
    This function reads in all hydrated fragments from the supplied folder and parses these into
    format readily accessible by solvate.

    Parameters
    ----------
    NONE

    Returns
    -------
    dictionary : ret
        This dictionary contains all files as keys to sub-dictionaries, which have two keys each, waters and protein. These
        hold the atom data for the respective atoms.

    """
    ### Initialise variables
    ret                                               = {}
    
    ### Report progress
    solvate_log.writeLog                              ( "Starting fragment reading from path " + str( solvate_globals._fragInputDir ) + ".", 2 )
    
    ### Find all input hydrated residue files
    allFragFiles                                      = [ fr for fr in glob.glob ( os.path.join ( os.path.join ( solvate_globals._fragInputDir, "*" ), "*.pdb" ), recursive = True ) ]
    
    ### If no frags, this is probably an error
    if len ( allFragFiles ) == 0:
        # Print error
        solvate_log.writeLog                          ( "!!! ERROR !!! Could not find any hydrated fragments. Please use the --frag option to supply the path to the hydrated fragments folder location.", 0 )
        
        # Terminate
        solvate_log.endLog                            ( )
        
    ### Read in all residue fragments
    for frag in allFragFiles:
            
        # Read it in
        ( protein, waters )                           = readFragment ( frag )
        
        # Save it
        ret[frag]                                     = {}
        ret[frag]["protein"]                          = protein
        ret[frag]["waters"]                           = waters
        
    ### Report progress
    solvate_log.writeLog                              ( "Fragments parsed.", 3 )
    
    ### Done
    return                                            ( ret )
