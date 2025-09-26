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
#%% Import configure file
#_________________________________________________________________________________________

import configparser

config = configparser.ConfigParser()
config.read('./structureFSISetup/structureInputPara.ini')

#_________________________________________________________________________________________
#
#%% If iMUICoupling, initialise MPI by mpi4py/MUI for parallelised computation
#_________________________________________________________________________________________

if config['MUI'].getboolean('iMUICoupling'):
    import sys
    from mpi4py import MPI
    import mui4py
    import petsc4py
    import os

    # App common world claims
    LOCAL_COMM_WORLD = mui4py.mpi_split_by_app()
    # MUI parameters
    dimensionMUI = 3
    data_types = {"dispX": mui4py.FLOAT64,
                  "dispY": mui4py.FLOAT64,
                  "dispZ": mui4py.FLOAT64,
                  "forceX": mui4py.FLOAT64,
                  "forceY": mui4py.FLOAT64,
                  "forceZ": mui4py.FLOAT64}
    # MUI interface creation
    domain = "structureDomain"
    config3d = mui4py.Config(dimensionMUI, mui4py.FLOAT64)

    iface = ["threeDInterface0"]
    ifaces3d = mui4py.create_unifaces(domain, iface, config3d)
    ifaces3d["threeDInterface0"].set_data_types(data_types)

    # Necessary to avoid hangs at PETSc vector communication
    petsc4py.init(comm=LOCAL_COMM_WORLD)

    # Define local communicator rank
    rank = LOCAL_COMM_WORLD.Get_rank()

    # Define local communicator size
    size = LOCAL_COMM_WORLD.Get_size()

#_________________________________________________________________________________________
#
#%% Import packages
#_________________________________________________________________________________________

from dolfinx import *
import structureFSISetup
import structureFSISolver

#_________________________________________________________________________________________
#
#%% Create instances for sub-domains and boundary condition
#_________________________________________________________________________________________

# Create sub-domain instances
subDomains = structureFSISetup.structureSubDomain.SubDomains()
# Create boundary condition instances
BCs = structureFSISetup.structureBCS.boundaryConditions()

#_________________________________________________________________________________________
#
#%% Create solver instances
#_________________________________________________________________________________________

solver = structureFSISolver.structureFSISolver.StructureFSISolver(config, subDomains, BCs)

#_________________________________________________________________________________________
#
#%% Solving
#_________________________________________________________________________________________

if config['MUI'].getboolean('iMUICoupling'):
    solver.solve(LOCAL_COMM_WORLD, ifaces3d)
else:
    solver.solve()

#_________________________________________________________________________________________
#
#%% Ploting
#_________________________________________________________________________________________

import matplotlib.pyplot as plt
import numpy as np

x = []
y = []
xBench = []
yBench = []
xCpp = []
yCpp = []

for i in np.arange(0.1, 100.1, 0.1):
  x.append(i)
for line in open('structureResults/tip-displacementY_0.txt', 'r'):
    lines = [i for i in line.split()]
    y.append(float(lines[0]))

for line in open('dataInput/Slone_et_al.txt', 'r'):
    lines = [i for i in line.split()]
    xBench.append(float(lines[0]))
    yBench.append(float(lines[1]))

plt.title("Y-Disp Compare")
plt.xlabel('Time [s]')
plt.ylabel('Y-Disp [m]')
plt.plot(xBench, yBench, label = 'Slone et al. 2003', marker= 'o', linestyle='None', c = 'b')
plt.plot(x, y, label = 'Present FEniCS Output', linestyle='-', c = 'g')
plt.xticks(np.arange(0, 101, step=20))
plt.yticks(np.arange(-0.15, 0.16, step=0.05))
plt.legend(loc='upper right')
plt.savefig('../result_compare.png')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FILE END  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#