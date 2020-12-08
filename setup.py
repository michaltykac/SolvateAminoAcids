# -*- coding: UTF-8 -*-
#   \file setup.py
#   \brief This is the python script driving the pip installation of the SolvateAminoAcids module.
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
#   \date      OCT 2020
######################################################

### Include requirements
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import os
import sys
import setuptools

### Define version
__version__                                           = '0.1.0'

### Set README.md as long description
this_directory                                        = os.path.abspath ( os.path.dirname ( __file__ ) )
with open ( os.path.join ( this_directory, 'README.md' ) ) as f:
    long_description                                  = f.read ( )

### Setup function call
setup ( name                                          = 'SolvateAminoAcids',
        version                                       = __version__,
        author                                        = 'Michal Tykac, Lada Biedermannova, Jiri Cerny',
        author_email                                  = 'Michal.Tykac@gmail.com',
        url                                           = 'https://github.com/michaltykac/SolvateAminoAcids',
        description                                   = 'SolvateAminoAcids is a Python language module for predicting water molecules in macromolecular structures.',
        long_description_content_type                 = "text/markdown",
        long_description                              = long_description,
        packages                                      = ['solvate'],
        setup_requires                                = [ 'numpy', 'gemmi' ],
        zip_safe                                      = False,
        classifiers                                   = [ 'Development Status :: 3 - Alpha',
                                                          'Intended Audience :: Science/Research',
                                                          'License :: OSI Approved :: BSD License',
                                                          'Operating System :: POSIX :: Linux',
                                                          'Operating System :: MacOS',
                                                          'Programming Language :: Python',
                                                          'Topic :: Scientific/Engineering :: Bio-Informatics' ],
        keywords                                      = 'bioinformatics computational-biology structural-biology water-prediction solvatation',
        project_urls                                  = { 'Github': 'https://github.com/michaltykac/SolvateAminoAcids',
                                                          'Bug Reports': 'https://github.com/michaltykac/SolvateAminoAcids/issues' }
      )
