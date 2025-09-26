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

    @file DOFCoordMapping.py
    
    @author W. Liu
    
    @brief DOF Coord Mapping file of the structure code.

"""

#_________________________________________________________________________________________
#
#%% Import packages
#_________________________________________________________________________________________
from dolfinx import *
import numpy as np

class DOFCoordMapping:

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #%% Define DOFs extract function
    #%% and
    #%% DOFs-Coordinates mapping function
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def dofs_to_xyz(self, FunctionSpace, dimension):
        # Convert dofs to coordinates
        return FunctionSpace.tabulate_dof_coordinates()

    def dofs_list(self, MeshFunction, FunctionSpace, boundary):
        return list(self.flexdofs)

    def xyz_np(self, dofs_list, FunctionSpace, dimension):
        xyz_np = np.zeros((len(dofs_list), dimension))
        for i, p in enumerate(dofs_list):
            xyz_np[i] = self.dofs_to_xyz(FunctionSpace, dimension)[p]
        return xyz_np
