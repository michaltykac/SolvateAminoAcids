# -*- coding: UTF-8 -*-
#   \file solvate_predictWaters.py
#   \brief This file provides water molecules predicting algorithms.
#
#   ...
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

import numpy                                          ### Maths
import solvate.solvate_log as solvate_log             ### Writing log
import solvate.solvate_maths as solvate_maths         ### Computing distances and rotations
import gemmi                                          ### File reading and writing
import os

######################################################
# predictWaters ()
def predictWaters ( resList, matchedFrags, fragFragments, settings ):
    """
    ...

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
    NONE

    """
    ### Log progress
    solvate_log.writeLog                              ( "Starting water prediction", settings, 1 )

    ### For each residue
    resCounter                                        = 1
    for res in range ( 0, len ( resList ) ):
    
        ### Report progress
        solvate_log.writeLog                          ( "Starting water prediction for residue " + str( resList[res][0] ), settings, 2 )
        
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
                matchStr                              = None
                if settings.matchedFragsPath != "":
                    
                    ### Create new structure
                    matchStr                          = gemmi.Structure ( )
                    
                    ### Initialise residue chain variables
                    resModel                          = gemmi.Model ( "1" )
                    resChain                          = gemmi.Chain ( "A" )
                    resRes                            = gemmi.Residue ( )
                    
                    ### For each residue atom
                    for at in resList[res][1]:
                    
                        ### Set the atom info
                        resAtom                       = gemmi.Atom ( )
                        resAtom.name                  = at[0]
                        resAtom.pos                   = gemmi.Position ( at[1], at[2], at[3] )
                        resAtom.occ                   = at[4]
                        resAtom.b_iso                 = at[5]
                        resRes.add_atom               ( resAtom )
                        
                    ### Set the residue info and add to chain and model
                    resRes.name                       = resList[res][0]
                    resRes.seqid                      = gemmi.SeqId ( 1, ' ' )
                    resRes.entity_type                = gemmi.EntityType.Polymer
                    resChain.add_residue              ( resRes )
                    resModel.add_chain                ( resChain )
                    
                    ### Initialise fragment chain variables
                    fragChain                         = gemmi.Chain ( "B" )
                    fragRes                           = gemmi.Residue ( )
                    
                    ### For each fragment protein atom
                    for at in fragFragments[frName]["protein"]:
                    
                        ### Rotate and translate each protein atom
                        fragAtomPos                   = numpy.array ( [ at[1],
                                                                        at[2],
                                                                        at[3] ] )
                        fragAtomPos                   = numpy.matmul ( fragAtomPos, matchedFrags[res][frag][frName]["rotation"] )
                        fragAtomPos                   = fragAtomPos + matchedFrags[res][frag][frName]["translation"]
                        
                        ### Set the atom info
                        fragAtom                      = gemmi.Atom ( )
                        fragAtom.name                 = at[0]
                        fragAtom.pos                  = gemmi.Position ( fragAtomPos[0], fragAtomPos[1], fragAtomPos[2] )
                        fragAtom.occ                  = at[4]
                        fragAtom.b_iso                = at[5]
                        fragRes.add_atom              ( fragAtom )
                        
                    ### Set the residue information and add to chain and model and structure
                    fragRes.name                      = resList[res][0]
                    fragRes.seqid                     = gemmi.SeqId ( 1, ' ' )
                    fragRes.entity_type               = gemmi.EntityType.Polymer
                    fragChain.add_residue             ( fragRes )
                    resModel.add_chain                ( fragChain )
                
                ### For each water atom
                matchedWaters                         = []
                for wIt in range ( 0, len ( fragFragments[frName]["waters"] ) ):
                
                    ### Rotate and translate each water atom
                    waterAtomPos                      = numpy.array ( [ fragFragments[frName]["waters"][wIt][1],
                                                                        fragFragments[frName]["waters"][wIt][2],
                                                                        fragFragments[frName]["waters"][wIt][3] ] )
                    waterAtomPos                      = numpy.matmul ( waterAtomPos, matchedFrags[res][frag][frName]["rotation"] )
                    waterAtomPos                      = waterAtomPos + matchedFrags[res][frag][frName]["translation"]
                    waterAtomPos                      = numpy.append ( waterAtomPos, numpy.array ( [ fragFragments[frName]["waters"][wIt][0],
                                                                                                     fragFragments[frName]["waters"][wIt][4],
                                                                                                     fragFragments[frName]["waters"][wIt][5] ] ) )
                    matchedWaters.append              ( waterAtomPos )
                
                ### Add waters and fragment to structure and write it, if required
                if matchStr is not None:
                
                    ### Create water chain and add it
                    watChain                          = gemmi.Chain ( "W" )
                    
                    ### For each water atom
                    resCtr                            = 1
                    for at in matchedWaters:
                    
                        ### Set the atom info
                        watAtom                       = gemmi.Atom ( )
                        watAtom.name                  = str ( at[3] )
                        watAtom.pos                   = gemmi.Position ( float ( at[0] ) , float ( at[1] ), float ( at[2] ) )
                        watAtom.occ                   = float ( at[4] )
                        watAtom.b_iso                 = float ( at[5] )
                    
                        ### Set the residue information and add to chain
                        watRes                        = gemmi.Residue ( )
                        watRes.add_atom               ( watAtom )
                        watRes.name                   = "HOH"
                        watRes.seqid                  = gemmi.SeqId ( resCtr, ' ' )
                        watRes.entity_type            = gemmi.EntityType.Water
                        watChain.add_residue          ( watRes )
                        resCtr                        = resCtr + 1
                        
                    
                    resModel.add_chain                ( watChain )
                    matchStr.add_model                ( resModel )
                    
                    ### Write fragment-residue match
                    matchName                         = str ( resCounter ) + "_" + str ( resList[res][0] ) + "-" + str ( frName.split( os.path.sep )[ len ( frName.split( os.path.sep ) ) - 1 ].split(".")[0] ) + ".pdb"
                    matchStr.write_pdb                ( os.path.join ( settings.matchedFragsPath, matchName ) )
                    
                    ### Report progress
                    solvate_log.writeLog                  ( "Written matched hydrated fragment to residue alignement file with water to " + str( frName ), settings, 4 )
                    
                ### Save the rotated and translated waters for this match
                
        resCounter                                    = resCounter + 1

