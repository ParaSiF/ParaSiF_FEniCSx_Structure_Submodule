"""
##############################################################################
# Parallel Partitioned Multi-Physics Simulation Framework (ParaSiF)          #
#                                                                            #
# Copyright (C) 2025 The ParaSiF Development Team                            #
# All rights reserved                                                        #
#                                                                            #
# This software is licensed under the GNU General Public License version 3   #
#                                                                            #
# ** GNU General Public License, version 3 **                                #
#                                                                            #
# This program is free software: you can redistribute it and/or modify       #
# it under the terms of the GNU General Public License as published by       #
# the Free Software Foundation, either version 3 of the License, or          #
# (at your option) any later version.                                        #
#                                                                            #
# This program is distributed in the hope that it will be useful,            #
# but WITHOUT ANY WARRANTY; without even the implied warranty of             #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              #
# GNU General Public License for more details.                               #
#                                                                            #
# You should have received a copy of the GNU General Public License          #
# along with this program.  If not, see <http://www.gnu.org/licenses/>.      #
##############################################################################

"""

#_________________________________________________________________________________________
#
#%% Import packages
#_________________________________________________________________________________________

from dolfinx import *
import configparser
import numpy as np

#_________________________________________________________________________________________
#
#%% Import configure file
#_________________________________________________________________________________________

config = configparser.ConfigParser()
config.read('./structureFSISetup/structureInputPara.ini')

OBeamXtemp=float(config['GEOMETRY']['OBeamX'])
OBeamYtemp=float(config['GEOMETRY']['OBeamY'])
OBeamZtemp=float(config['GEOMETRY']['OBeamZ'])
XBeamtemp=float(config['GEOMETRY']['XBeam'])
YBeamtemp=float(config['GEOMETRY']['YBeam'])
ZBeamtemp=float(config['GEOMETRY']['ZBeam'])

#_________________________________________________________________________________________
#
#%% Define SubDomains classes
#%% for defining parts of the boundaries and the interior of the domain
#_________________________________________________________________________________________

class SubDomains:
    def Fixed (self , x ):
        return np.isclose(x[0], (OBeamXtemp))

    def Flex (self , x ):
        return np.isclose(x[0], (OBeamXtemp + XBeamtemp))

    def Symmetry (self , x ):
        return np.logical_or(np.logical_or(np.isclose(x[2], (OBeamZtemp + ZBeamtemp)), np.isclose(x[2], (OBeamZtemp))), np.logical_or(np.isclose(x[1], (OBeamYtemp)), np.isclose(x[1], (OBeamYtemp + YBeamtemp))))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FILE END  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#