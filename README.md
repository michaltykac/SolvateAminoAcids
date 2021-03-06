# SolvateAminoAcids

SolvateAminoAcids is a Python language module for predicting water molecules in macromolecular structures.

## Description

SolvateAminoAcids uses a method for fast addition of water molecules to pdb files which may not have any water molecules or which do not have all water molecules added yet. The method is as described in:

    Černý J., Schneider B. and Biedermannová L. (2017) WatAA: Atlas of Protein Hydration. Exploring synergies between data mining and ab initio calculations. Physical Chemistry Chemical Physics 19(26), pp. 17094-17102.
    
## Index

- [Description](#description)
- [Index](#index)
- [Obtaining SolvateAminoAcids](#obtaining-solvateaminoacids)
- [Usage](#usage)
    -[Help dialogue](#help-dialogue)
    -[Default run](#default-run)
    -[More options](#more-options)
- [Description](#description)
- [Dependencies](#dependencies)
- [People behind SolvateAminoAcids](#people-behind-solvateaminoacids)
    - [Authors](#authors)
    - [Copyright](#copyright)
    - [License](#license)

## Obtaining SolvateAminoAcids

The latest version of SolvateAminoAcids can be cloned from its GitHub repository 

    https://github.com/michaltykac/SolvateAminoAcids

## Usage

### Help dialogue

After the GitHub repository has been cloned onto your machine, solvate can be started by running the top-level python script file called *solvate.py* from the command line.

    python ./solvate.py -h
    
This command will display the *help* screen for running the SolvateAminoAcids script in the standard GNU format. 

### Default run

Now, assuming that we have a macromolecular co-ordinates structure in the PDB format, for example 1LFA ( available from http://www.rcsb.org/structure/1LFA ), that you have copied into a folder called *solvateaa_test*, you can supply this structure to the solvate script by issuing the following command:

    python /path/to/SolvateAminoAcids/folded/solvate.py -i ./1lfa.pdb --frag /path/to/SolvateAminoAcids/folded/fragData --res /path/to/SolvateAminoAcids/folded/resData
    
where the *--frag* and the *--res* command line options could be avoided if the current working directory is the same as */path/to/SolvateAminoAcids/folded*. 

The result of this command will be a default run of SolvateAminoAcids script, which will give following output

    ... Starting to read structure ./1lfa.pdb
    ... Starting to parse structure for speed-up
    ... Starting fragment matching
    ... Starting water prediction
    ... Removing clashes between new waters and structure contents
    ... Clustering predicted non-clashing water molecules
    ... Started combining and adding waters.
    ... Writing out the modified co-ordinate files
    
And which will create the *solvate_log.txt* file containing the log for the run and also the *solvate_fullStructure_waters.pdb* file, which is a macromolecular co-ordinate file with identical protein structure to the input file (*1lfa.pdb*), but containing the predicted water molecules. The following picture shows the original structure on the left and the new structure with predicted waters on the right.

![](https://github.com/michaltykac/SolvateAminoAcids/blob/solvate_amino_acids/fig1_1lfa.png)

### More options

While the default run shown above is capable of providing predictions as described in the paper, the SolvateAminoAcids module has more features. Firstly, there are several thresholds that SolvateAminoAcids uses to make decisions; while care has been taken to select a reasonable default values, the user is invited to change these to their liking. These thresholds are as follows:

- ```r``` Fragment RMSD threshold - only hydrated fragments with RMSD threshold below this value will be considered similar enough.
- ```clRad``` The clash radius - there must be no other atom within this radius for the water molecule to be added.
- ```maxClDist``` Maximum cluster distance - when water molecules are clustered, only water molecules closer than this distance will be clustered together.

Moreover, SolvateAminoAcids also has two distinct ways of selecting which fragment should be considered to hydrate each residue. The default way is to try all fragments available for the processed residue and select all fragments with RMSD value to the residue under the threshold ```r```. Should the user require only the best fragment to be used, this can be arranged using the command line option ```bestOnly```, while if the user requires SolvateAminoAcids to first determine the backbone secondary-structure and side-chain rotamers and then testing only hydrated fragments conforming to these constraints, this can be achieved by supplying the command line option ```b```.

Finally, SolvateAminoAcids will be default use the structure as supplied, including ligands, already present waters and hydrogens. Should the user require the structure to be modified in any of the following ways, it can be requested using the appropriate command line optino. Please note that multiple structure modifications do not combine, but rather each either allows or disables a particular modified structure run; multiple structure modifications can be processed with a single SolvateAminoAcids run.

- ```noFullStr``` By default, SolvateAminoAcids will use all structure contents to detect clashes of any predicted waters. This option turns that off.
- ```noHydro``` Should SolvateAminoAcids produce another set of results ignoring the hydrogen atoms clashes with the predicted waters?
- ```noWaters``` Should SolvateAminoAcids produce another set of results ignoring the already present water molecule clashes with the predicted waters?
- ```noLigand``` Should SolvateAminoAcids produce another set of results ignoring the ligand molecule clashes with the predicted waters?


## Dependencies

SolvateAminoAcids requires a few python modules to be installed in order to work. They can all be installed from their respective PyPi repositories:

- **numpy** ( *python -m pip install numpy* )
- **sklearn** ( *python -m pip install sklearn* )
- **gemmi** ( *python -m pip install gemmi* )

## People behind SolvateAminoAcids

### Authors

- Michal Tykac
- Lada Biedermannová
- Jiří Černý

### Copyright

Copyright by the [Authors](#authors) and individual contributors. All rights reserved.

### License

SolvateAminoAcids is available under the 3-clause BSD license stated as follows:

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
    1) Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
    2) Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
    3) Neither the name of Michal Tykac, Lada Biedermannová or Jiří Černý, nor the names of this code's contributors may be used to endorse or promote products derived from this software without specific prior written permission.

This software is provided by the copyright holder and contributors "as is" and any express or implied warranties, including, but not limitted to, the implied warranties of merchantibility and fitness for a particular purpose are disclaimed. In no event shall the copyright owner or the contributors be liable for any direct, indirect, incidental, special, exemplary, or consequential damages (including, but not limitted to, procurement of substitute goods or services, loss of use, data or profits, or business interuption) however caused and on any theory of liability, whether in contract, strict liability or tort (including negligence or otherwise) arising in any way out of the use of this software, even if advised of the possibility of such damage.
