# -*- coding: UTF-8 -*-
#   \file solvate_commandLineArgs.py
#   \brief This file provides functionality for command line parsing.
#
#   This file contains functions and settings used to parse command line arguments. None of these functions should
#   be directly accessible from the python module or the python tool, they are for internal purposes only.
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

import argparse                                       ### Command line arguments parser
import sys                                            ### In the case complete stop is needed
import os                                             ### Checking for file existence

######################################################
# Description
parserDesc                                            = """Solvate is a Python language module for assigning water molecules to macromolecular structures.

Solvate uses a method for fast addition of water molecules to pdb files which may not have any
water molecules or which do not have all water molecules added yet. The method is as described
in

    Černý J., Schneider B. and Biedermannová L. (2017) WatAA: Atlas of Protein Hydration. Exploring
    synergies between data mining and ab initio calculations. Physical Chemistry Chemical Physics 19(26),
    pp. 17094-17102.

Written by: Michal Tykac for BTU AVCR in 2019 and 2020."""

######################################################
# getCLArgs function
def getCLArgs ( version ):
    """
    This function deals with all the command line parsing details and returns the Namespace with all available
    values.

    Parameters
    ----------
    str : version
        String defining the version of the tool.

    Returns
    -------
    Args: Namespace
        This is the args structure which holds the values supplied by the user, which are to be assigned to the correct variables.

    """
    ### Initialise Parser
    parser                                            = argparse.ArgumentParser ( description = parserDesc, formatter_class = argparse.RawTextHelpFormatter, prog = "solvate" )

    ### Supply the input co-ordinate file
    parser.add_argument ( '-i',
                          type                        = str,
                          nargs                       = 1,
                          required                    = True,
                          help                        = 'The location and file name of the co-ordinate file which will have waters added to it.')

#    ### Output file name
#    parser.add_argument ( '-O',
#                          type                        = str,
#                          nargs                       = 1,
#                          help                        = 'The path and file name for the output PDB file. Please do not add the extension, as it will be generated automatically. The default value is \"./solvate_res\".')

    ### RMSD fragment matching threshold
    parser.add_argument ( '--verbose',
                          type                        = int,
                          nargs                       = 1,
                          help                        = 'How verbose should the run be? 0 = nothing 5 = max')

    ### RMSD fragment matching threshold
    parser.add_argument ( '-r',
                          type                        = float,
                          nargs                       = 1,
                          help                        = 'The RMSD value threshold for fragment to be considered \"similar enough\" to the residue. The default value is 0.5.')

    ### Use backbone to determine secondary structure
    parser.add_argument ( '-b',
                          action                      = 'store_true',
                          help                        = 'If this option is used, hydrated residues will be used to determine secondary structure and rotamer, limiting the range of tested hydrated fragments.' )
                          
    ### Input fragments folder
    parser.add_argument ( '--frag',
                          type                        = str,
                          nargs                       = 1,
                          help                        = 'The path to the folder containing the hydrated fragment input files. Defaults to \"fragData\".')
                          
    ### Input matched residues location
    parser.add_argument ( '--res',
                          type                        = str,
                          nargs                       = 1,
                          help                        = 'The path to the folder containing the hydrated residues input files. Defaults to \"resData\".' )
   
    parser.add_argument ( '--bestOnly',
                          action                      = 'store_true',
                          help                        = 'If used, only the best fragment will be used as a match to the residue.' )
   
    ### Change the log file name
    parser.add_argument ( '--log',
                          type                        = str,
                          nargs                       = 1,
                          help                        = 'The path to where the log will be written to. Defaults to \"solvate_log.txt\".' )
    
    ### Version
    parser.add_argument('--version', action = 'version', version = '%(prog)s ' + str( version ) )

    ### Parse the input command line arguments
    args                                              = parser.parse_args()
    
    ### Check for mandatory arguments
    if args.i is None or args.i[0] == "" or args.i[0] == " ":
        sys.exit ( 'Error: Missing input file. Please supply a co-ordinate file using the -i command line option or use -h for help.' )
        
    if not os.path.isfile ( args.i[0] ):
        sys.exit ( 'Error: The supplied co-ordinate file ' + args.i[0] + ' does not seem to exist. Please supply a different file.' )
    
    return ( args )
