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

    @file structureFSISolver.py

    @author W. Liu

    @brief The main structure solver class of the structure code.

"""

#_________________________________________________________________________________________
#
#%% Import packages
#_________________________________________________________________________________________

from dolfinx import *
from mpi4py import MPI
import configparser
import math
import numpy as np
import os
import datetime
import socket
import sys
import structureFSISolver

import structureFSISolver.functions.cfgPrsFn
import structureFSISolver.functions.lameParm
import structureFSISolver.functions.meshBoundarySubdomian
import structureFSISolver.functions.meshMotion
import structureFSISolver.functions.compilerOpt
import structureFSISolver.functions.utility
import structureFSISolver.functions.checkpointLogCtrl
import structureFSISolver.functions.couplingMUIFn
import structureFSISolver.functions.DOFCoordMapping
import structureFSISolver.functions.facetAreas
import structureFSISolver.functions.facetAreas
import structureFSISolver.functions.timeMarching
import structureFSISolver.solvers.linearElasticSolver
import structureFSISolver.solvers.hyperElasticSolver

#_________________________________________________________________________________________
#
#%% Main Structure Solver Class
#_________________________________________________________________________________________

class StructureFSISolver(structureFSISolver.functions.cfgPrsFn.readData,
                         structureFSISolver.functions.meshBoundarySubdomian.meshBoundarySubdomian,
                         structureFSISolver.functions.timeMarching.timeMarching,
                         structureFSISolver.functions.lameParm.lameParm,
                         structureFSISolver.functions.meshMotion.meshMotion,
                         structureFSISolver.functions.compilerOpt.compilerOpt,
                         structureFSISolver.functions.utility.utility,
                         structureFSISolver.functions.checkpointLogCtrl.checkpointLogCtrl,
                         structureFSISolver.functions.couplingMUIFn.couplingMUIFn,
                         structureFSISolver.functions.DOFCoordMapping.DOFCoordMapping,
                         structureFSISolver.functions.facetAreas.facetAreas,
                         structureFSISolver.solvers.linearElasticSolver.linearElastic,
                         structureFSISolver.solvers.hyperElasticSolver.hyperElastic):
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #%% Solver initialize
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def __init__(self,
                 Configure,
                 SubDomains,
                 DirichletBoundaryConditions):

        #===========================================
        #%% Obtain files and instances
        #===========================================

        # Obtain configure file
        self.cfg = Configure
        # Obtain sub-domain instances
        self.subDomains = SubDomains
        # Obtain dirichlet boundary condition instances
        self.dirichletBCs = DirichletBoundaryConditions

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #%% Main solver function
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def solve(self, arg_COMM_WORLD=None, arg_ifaces3d=None):
        #===========================================
        #%% Setup the wall clock
        #===========================================

        # create an instance of the TicToc wall clock class
        wallClock = structureFSISolver.tictoc.TicToc()
        # Starts the wall clock
        wallClock.tic()

        #===========================================
        #%% Initialise MPI by mpi4py/MUI for
        #%%   parallelised computation
        #===========================================

        if arg_COMM_WORLD is not None:
            self.MUI_Init(arg_COMM_WORLD, arg_ifaces3d)
        else:
            self.MUI_Init()

        #===========================================
        #%% Set target folder
        #===========================================

        # Folder directory
        if self.iAbspath():
            self.outputFolderPath = os.path.abspath(self.outputFolderName())
            self.inputFolderPath = os.path.abspath(self.inputFolderName())
        else:
            self.outputFolderPath = self.outputFolderName()
            self.inputFolderPath = self.inputFolderName()

        #===========================================
        #%% Print log information
        #===========================================

        self.Pre_Solving_Log()

        #===========================================
        #%% Set form compiler options
        #===========================================
        #  !! OUTDATED FUNCTION, NEED UPDATED TO FENICS-X !!
        #self.Set_Compiler_Options()

        #===========================================
        #%% Time marching parameters define
        #===========================================

        self.Time_Marching_Parameters()
        self.Time_Marching_Log()

        #===========================================
        #%% Call solvers
        #===========================================

        if self.solving_method() == 'STVK':
            self.hyperElasticSolve()
        elif self.solving_method() == 'MCK':
            if arg_ifaces3d is not None:
                self.linearElasticSolve(arg_ifaces3d)
            else:
                self.linearElasticSolve()
        else:
            sys.exit("{FENICS} Error, solving method not recognised")
        #===========================================
        #%% Calculate wall time
        #===========================================

        # Finish the wall clock
        simtime = wallClock.toc()
        self.Post_Solving_Log(simtime)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FILE END  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#