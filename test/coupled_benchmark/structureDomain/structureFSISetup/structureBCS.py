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
import numpy as np

#_________________________________________________________________________________________
#
#%% Define boundary conditions
#_________________________________________________________________________________________

class boundaryConditions:
    def DirichletMixedBCs(self, MixedVectorFunctionSpace, boundaries, marks):
        #  !! OUTDATED FUNCTION, NEED UPDATED TO FENICS-X !!
        bc1 = DirichletBC(MixedVectorFunctionSpace.sub(0), ((0.0,0.0,0.0)),boundaries, marks)
        bc2 = DirichletBC(MixedVectorFunctionSpace.sub(1), ((0.0,0.0,0.0)),boundaries, marks)
        return bc1, bc2
    def DirichletBCs(self, VectorFunctionSpace, boundary_dofs):
        bc3 = fem.dirichletbc(np.zeros(3),
                              boundary_dofs,
                              V=VectorFunctionSpace)
        return bc3

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FILE END  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#