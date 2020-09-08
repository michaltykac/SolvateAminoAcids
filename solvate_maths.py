# -*- coding: UTF-8 -*-
#   \file solvate_maths.py
#   \brief This file provides the specific computations functionality for other parts of solvate.
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

import numpy                                          ### Maths

import solvate_globals                                ### Global variables
import solvate_log                                    ### For writing the log

######################################################
# residueToResFragmentDistance_backbone ()
def residueToResFragmentDistance_backbone ( res, resFrag ):
    """
    This function computes the procrustes optimal overlay of the backbone atoms of the input residue and
    hydrated residue fragment and then proceeds to obtain the RMSD of these backbone atoms.

    Parameters
    ----------
    list : res
        A list containing the positions and names of the compared residue atoms.
        
    list : resFrag
        A list containing the positions and names of the hydrated residue frament atoms.

    Returns
    -------
    float : dist
        The RMSD distance between the supplied residues's backbone atoms.

    """
    ### Reduce residue to N, Ca, C and O (backbone) atoms
    resAts                                            = []
    resAtsNames                                       = []
    for rIt in range( 0, len( res[1] ) ):
        if ( res[1][rIt][0] == "C" ) or ( res[1][rIt][0] == "CA" ) or ( res[1][rIt][0] == "O" ) or ( res[1][rIt][0] == "N" ):
            resAts.append                             ( [ res[1][rIt][1], res[1][rIt][2], res[1][rIt][3] ] )
            resAtsNames.append                        ( res[1][rIt][0] )

    ### Sanity check
    if len( resAts ) != 4:
        solvate_log.writeLog                          ( "!!! ERROR !!! Could not find all backbone atoms (C, CA, N and O) in the residue described by " + str( res ), 0 )
        
        # Terminate
        solvate_log.endLog                            ( )
    
    ### Reduce hydrated residue fragment to N, Ca, C and O (backbone) atoms
    fraAts                                            = []
    for fIt in range( 0, len( resFrag ) ):
        for rIt in range ( 0, len( resAtsNames ) ):
            if resFrag[fIt][0] == resAtsNames[rIt]:
                fraAts.append                         ( [ resFrag[fIt][1], resFrag[fIt][2], resFrag[fIt][3] ] )
        
    ### Sanity check
    if len( fraAts ) != 4:
        solvate_log.writeLog                          ( "!!! ERROR !!! Could not find all backbone atoms (C, CA, N and O) in the supplied hydrated fragment " + str ( resFrag ), 0 )
        
        # Terminate
        solvate_log.endLog                            ( )
    
    ### Convert to numpy arrays
    resNumpy                                          = numpy.array ( [ numpy.array ( xi ) for xi in resAts ] )
    fraNumpy                                          = numpy.array ( [ numpy.array ( xi ) for xi in fraAts ] )
    
    ### Compute procrustes analysis
    ( dist, z, t )                                    = procrustes ( resNumpy, fraNumpy )
    
    ### Get RMSD
    rmsdDist                                          = getRMSD ( resAts, z )
    
    ### Report log
    solvate_log.writeLog                              ( "Found residue to hydrated residue fragment distance of " + str( rmsdDist ), 4 )
    
    ### Done
    return                                            ( rmsdDist )
    
######################################################
# residueToResFragmentDistance_rotamer ()
def residueToResFragmentDistance_rotamer ( res, resFrag ):
    """
    This function computes the procrustes optimal overlay of the side-chain atoms of the input residue and
    hydrated residue fragment and then proceeds to obtain the RMSD of these side-chain atoms.

    Parameters
    ----------
    list : res
        A list containing the positions and names of the compared residue atoms.
        
    list : resFrag
        A list containing the positions and names of the hydrated residue frament atoms.

    Returns
    -------
    float : dist
        The RMSD distance between the supplied residues's side-chain atoms.

    """
    ### Initialise variables
    resAts                                            = []
    resAtsNames                                       = []
    fraAts                                            = []
    
    ### Parse out the side-chain atoms for the residue
    for rIt in range( 0, len( res[1] ) ):
    
        # Deal with various atom names and lengths for different amino acids
        if ( ( ( res[0] == "VAL" ) or ( res[0] == "ILE" ) ) and ( ( res[1][rIt][0] == "C" ) or ( res[1][rIt][0] == "CA" ) or ( res[1][rIt][0] == "CB" ) or ( res[1][rIt][0] == "CG1" ) ) ) or \
           ( ( ( res[0] == "THR" )                        ) and ( ( res[1][rIt][0] == "C" ) or ( res[1][rIt][0] == "CA" ) or ( res[1][rIt][0] == "CB" ) or ( res[1][rIt][0] == "OG1" ) ) ) or \
           ( ( ( res[0] == "SER" )                        ) and ( ( res[1][rIt][0] == "C" ) or ( res[1][rIt][0] == "CA" ) or ( res[1][rIt][0] == "CB" ) or ( res[1][rIt][0] == "OG"  ) ) ) or \
           ( ( ( res[0] == "CYS" )                        ) and ( ( res[1][rIt][0] == "C" ) or ( res[1][rIt][0] == "CA" ) or ( res[1][rIt][0] == "CB" ) or ( res[1][rIt][0] == "SG"  ) ) ) or \
           ( ( ( res[0] == "ARG" ) or ( res[0] == "ASN" ) or ( res[0] == "ASP" ) or ( res[0] == "GLN" ) or ( res[0] == "GLU" ) or ( res[0] == "HIS" ) or ( res[0] == "LEU" ) or \
               ( res[0] == "LYS" ) or ( res[0] == "MET" ) or ( res[0] == "PHE" ) or ( res[0] == "PRO" ) or ( res[0] == "TRP" ) or ( res[0] == "TYR" ) \
                                                          ) and ( ( res[1][rIt][0] == "C" ) or ( res[1][rIt][0] == "CA" ) or ( res[1][rIt][0] == "CB" ) or ( res[1][rIt][0] == "CG"  ) ) ):
        
            # Save
            resAts.append                             ( [ res[1][rIt][1], res[1][rIt][2], res[1][rIt][3] ] )
            resAtsNames.append                        ( res[1][rIt][0] )

    ### Sanity check
    if len ( resAts ) != 4:
        solvate_log.writeLog                          ( "!!! ERROR !!! Could not find all side-chain atoms (C, CA, CB and ?) in the supplied residue " + str ( res ), 0 )
        
        # Terminate
        solvate_log.endLog                            ( )

    
    ### Reduce fragment to N, Ca, Cb and X (as above) atoms
    for rIt in range( 0, len( resAtsNames ) ):
    
        # For each fragment atom
        for fIt in range( 0, len( resFrag ) ):
        
            # Get only the same atoms
            if resAtsNames[rIt] == resFrag[fIt][0]:
            
                # Save
                fraAts.append                         ( [ resFrag[fIt][1], resFrag[fIt][2], resFrag[fIt][3] ] )
                
    ### Sanity check
    if len ( fraAts ) != 4:
        solvate_log.writeLog                          ( "!!! ERROR !!! Could not find all side-chain atoms (C, CA, CB and ?) in the supplied hydrated residue fragment " + str ( resFrag ), 0 )
        
        # Terminate
        solvate_log.endLog                            ( )
        
    ### Convert to numpy arrays
    resNumpy                                          = numpy.array ( [ numpy.array ( xi ) for xi in resAts ] )
    fraNumpy                                          = numpy.array ( [ numpy.array ( xi ) for xi in fraAts ] )
    
    ### Compute procrustes analysis
    ( dist, z, t )                                    = procrustes ( resNumpy, fraNumpy )
    
    ### Get RMSD
    rmsdDist                                          = getRMSD ( resAts, z )
    
    ### Report log
    solvate_log.writeLog                              ( "Found residue to hydrated residue fragment distance of " + str( rmsdDist ), 4 )
    
    ### Done
    return                                            ( rmsdDist )

######################################################
# residueToResFragmentDistance_direct ()
def residueToResFragmentDistance_direct ( res, fragInfo ):
    """
    This function computes the procrustes optimal overlay of XXX ... XXX and then proceeds to obtain the RMSD of these atoms.

    Parameters
    ----------
    list : res
        A list containing the positions and names of the compared residue atoms.
        
    list : fragInfo
        A list containing the positions and names of the hydrated frament atoms.

    Returns
    -------
    float : dist
        The RMSD distance between the supplied residue and fragment atoms.

    """
    ### Initialise variables
    resAtoms                                          = []
    fragAtoms                                         = []

    ### Get residue positions in same order as fragment
    for fIt in range ( 0, len ( fragInfo ) ):
    
        # For each residue atom name
        for atIt in range ( 0, len ( res[1] ) ):
        
            # Compare the atom names
            if res[1][atIt][0] == fragInfo[fIt][0]:
            
                # Save atoms in same order
                resAtoms.append                       ( [ res[1][atIt][1], res[1][atIt][2], res[1][atIt][3] ] )
                fragAtoms.append                      ( [ fragInfo[fIt][1], fragInfo[fIt][2], fragInfo[fIt][3] ] )
    
    ### Check for missing atoms
    if len ( fragInfo ) != len ( resAtoms ):
        solvate_log.writeLog                          ( "!!! ERROR !!! Failed to match the fragment atoms of fragment " + str ( fragInfo ) + " to residue " + str( res ), 0 )
    
        # Terminate
        solvate_log.endLog                            ( )
        

    ### Convert to numpy arrays
    resNumpy                                          = numpy.array ( [ numpy.array ( xi ) for xi in resAtoms ] )
    fraNumpy                                          = numpy.array ( [ numpy.array ( xi ) for xi in fragAtoms ] )
    
    ### Compute procrustes analysis
    ( dist, z, t )                                    = procrustes ( resNumpy, fraNumpy )
    
    ### Get RMSD
    rmsdDist                                          = getRMSD ( resAtoms, z )
    
    ### Report log
    solvate_log.writeLog                              ( "Found residue to hydrated fragment distance of " + str( rmsdDist ), 4 )
    
    ### Done
    return                                            ( rmsdDist )

######################################################
# procrustes ()
def procrustes ( X, Y ):
    """
    A port of MATLAB's "procrustes" function to Numpy.
    
    Procrustes analysis determines a linear transformation (translation,
    reflection, orthogonal rotation and scaling) of the points in Y to best
    conform them to the points in matrix X, using the sum of squared errors
    as the goodness of fit criterion.
    
        d, Z, [tform] = procrustes ( X, Y )
        
    In this particular code, scaling has been disabled (in protein comparisons,
    scale plays a role) as was the reflections ( mirror immage of a molecule is
    a different molecule ).
    
    Parameters
    ----------
    numpy.ndarray : X, Y
        Matrices of target and input coordinates. They must have equal
        numbers of  points (rows), but Y may have fewer dimensions
        (columns) than X.
    
    Returns
    -------
    float : d
        The residual sum of squared errors, normalized according to a
        measure of the scale of X, ( ( X - X.mean(0))^2 ).sum( )
    
    matrix : Z
        The matrix of transformed Y-values
    
    dictionary : tform
        A dict specifying the rotation, translation and scaling that
        maps X --> Y
    
    """
    
    ### Find the dimensions of the input matrices
    nx , mx                                           = X.shape
    ny , my                                           = Y.shape
    
    ### Find the column means of the input matrices
    muX                                               = X.mean(0)
    muY                                               = Y.mean(0)
    
    ### Find the x - mean values of the input matrices
    X0                                                = X - muX
    Y0                                                = Y - muY
    
    ### Sum the squares of the x - mean values of the input matrices
    ssX                                               = ( numpy.power ( X0, 2.0 ) ).sum()
    ssY                                               = ( numpy.power ( Y0, 2.0 ) ).sum()
    
    ### Compute the centred Frobenius norm of the input matrices (i.e. square root of sum of x - mean squared)
    normX                                             = numpy.sqrt ( ssX )
    normY                                             = numpy.sqrt ( ssY )
    
    ### Scale x - mean to equal (unit) norm
    X0                                               /= normX
    Y0                                               /= normY
    
    ### If second matrix has more rows than the first, shorten it ( this should really be error... )
    if my < mx:
        Y0                                            = numpy.concatenate ( ( Y0, numpy.zeros ( nx, mx - my ) ), 0 )
    
    ### Find the optimal rotation matrix of Y using least squares
    A                                                 = numpy.dot ( X0.T, Y0 )
    U,s,Vt                                            = numpy.linalg.svd ( A, full_matrices=False )
    V                                                 = Vt.T
    T                                                 = numpy.dot(V, U.T)
    
    ### Does the current solution use a reflection?
    have_reflection                                   = numpy.linalg.det(T) < 0
    
    ### If that's the case, force the other reflection (i.e. solution with no reflections)
    if have_reflection != False:
        V[:,-1]                                      *= -1
        s[-1]                                        *= -1
        T                                             = numpy.dot(V, U.T)
        
    ### Find the trace of singular values
    traceTA                                           = s.sum()
    
    ### With no scaling
    b                                                 = 1
    d                                                 = 1 + ssY/ssX - 2 * traceTA * normY / normX
    Z                                                 = normY * numpy.dot ( Y0, T ) + muX
    
    ### If second matrix has more rows than the first, shorten it ( this should really be error... )
    if my < mx:
        T                                             = T[:my,:]
        
    ### Compute the translation vector
    c                                                 = muX - b * numpy.dot ( muY, T )
    
    ### Save all to transformation list
    tform                                             = {'rotation':T, 'scale':b, 'translation':c}
    
    ### Done
    return d, Z, tform

######################################################
# getRMSD ()
def getRMSD ( coo1, coo2 ):
    """
    This is a simple function for computing RMSD given two matrices with equal dimmensions (n rows and
    exactly 3 columns) of (presumably atomic) co-ordinates.
    
    Parameters
    ----------
    matrix : coo1, coo2
        matrices of target and input coordinates. They must have equal
        numbers of points (rows) and equal number of dimensions (columns
        = 3)
    
    Returns
    -------
    float : dist
        The RMSD distance value between the two matrices.
    
    """
    ### Initialise sum
    sum                                               = 0.0
    
    ### For each row
    for iter in range ( 0, len ( coo1 ) ):
        sum                                          += ( numpy.power ( coo1[iter][0] - coo2[iter][0], 2.0 ) +
                                                          numpy.power ( coo1[iter][1] - coo2[iter][1], 2.0 ) +
                                                          numpy.power ( coo1[iter][2] - coo2[iter][2], 2.0 ) )
                
    ### Done
    return                                            ( numpy.sqrt ( sum  / float ( len ( coo1 ) ) ) )
    

######################################################
# findSmallestDistance ()
def findSmallestDistance ( distances, threshold ):
    """
    This function takes a list of distances and a threshold and returns the index of the
    smallest distance or -1 if no distances is below the threshold.

    Parameters
    ----------
    list : distances
        The list of the distances to be searched.
        
    float : threshold
        The threshold below which any passing distance needs to be.

    Returns
    -------
    int : index
        The index of the smallest distance or -1 if no distance is below the threshold.

    """
    ### Find lowest RMSD and check if it is larger than the threshold
    minVal                                            = numpy.Inf
    minInd                                            = -1
    indCtr                                            = 0
    
    ### Find the best match (lowest distance)
    for val in distances:
        if ( val < minVal ) and ( val < threshold ):
            minVal                                    = val
            minInd                                    = indCtr
        indCtr                                       += 1
    
    ### Return the index
    return                                            ( minInd )

######################################################
# findPassingDistances ()
def findPassingDistances ( distances, threshold ):
    """
    This function takes a list of distances and a threshold and returns all indices of
    distances which pass the threshold, or -1 if no distances is below the threshold.

    Parameters
    ----------
    list : distances
        The list of the distances to be searched.
        
    float : threshold
        The threshold below which any passing distance needs to be.

    Returns
    -------
    list : indices
        The indices of all passing distances or -1 if no distance is below the threshold.

    """
    ### Initialise variables
    minInd                                            = []
    indCtr                                            = 0
    
    ### Find all passing distances
    for val in distances:
        if val < threshold:
            minInd.append                             ( indCtr )
        indCtr                                       += 1
        
    ### Check for at least one
    if len( minInd ) == 0:
        minInd                                        = -1
    
    ### Return the index
    return                                            ( minInd )
