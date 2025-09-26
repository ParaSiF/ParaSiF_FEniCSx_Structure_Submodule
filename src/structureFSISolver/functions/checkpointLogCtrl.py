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

    @file checkpointLogCtrl.py
    
    @author W. Liu
    
    @brief checkpoint and Log Control file of the structure code.

"""

#_________________________________________________________________________________________
#
#%% Import packages
#_________________________________________________________________________________________
import datetime
import socket
from dolfinx import *
from dolfinx.io import XDMFFile
import numpy as np
from mpi4py import MPI

class checkpointLogCtrl:

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #%% Print log information
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def Pre_Solving_Log(self):
        if self.rank == 0:
            print ("\n")
            print ("{FENICS} ********** STRUCTURAL-ELASTICITY SIMULATION BEGIN **********")
            print ("\n")
            if self.iDebug():
                print ("{FENICS} ### !!! DEBUG LEVEL ON !!! ###")
                print ("\n")

            if self.iMUICoupling():
                print ("{FENICS} ### !!! MUI COUPLING ON !!! ###")
                print ("\n")

            print ("{FENICS} Current Date and Time: ", datetime.datetime.now())
            print ("{FENICS} System Host Name: ", socket.gethostbyaddr(socket.gethostname())[0])
            print ("\n")

            print ("{FENICS} Solver info: ")
            if self.solving_method() == 'STVK':
                print ("{FENICS} Solver for the problem: ", self.prbsolver())
                print ("{FENICS} Solver for project between domains: ", self.prjsolver())
                print ("{FENICS} Pre-conditioner for the problem: ", self.prbpreconditioner())
                print ("{FENICS} Relative tolerance: ", self.prbRelative_tolerance())
                print ("{FENICS} Absolute tolerance: ", self.prbAbsolute_tolerance())
                print ("{FENICS} Maximum iterations: ", self.prbMaximum_iterations())
                print ("{FENICS} Relaxation parameter: ", self.prbRelaxation_parameter())
                print ("{FENICS} Representation of the compiler: ", self.compRepresentation())
                print ("{FENICS} C++ code optimization: ", self.cppOptimize())
                print ("{FENICS} optimization of the compiler: ", self.optimize())
                print ("{FENICS} Extrapolation: ", self.allow_extrapolation())
                print ("{FENICS} Ghost cell mode: ", self.ghost_mode())
                print ("{FENICS} Error of non convergence: ", self.error_on_nonconvergence())
            elif self.solving_method() == 'MCK':
                print ("{FENICS} Solver for the problem: ", self.prbsolver())
                print ("{FENICS} Solver for project between domains: ", self.prjsolver())
            print ("\n")

            print ("{FENICS} Input parameters: ")
            print ("{FENICS} E: ", self.E_s(), "[Pa]")
            print ("{FENICS} rho: ", self.rho_s(), "[kg/m^3]")
            print ("{FENICS} nu: ", self.nu_s(), "[-]")
            print ("\n")
        else:
            pass

    def Time_Marching_Log(self):
        if self.rank == 0: 
            print ("\n")
            print ("{FENICS} Total time: ", self.T(), " [s]")
            print ("{FENICS} Time step size: ", self.dt(), " [s]")
            print ("{FENICS} Time steps: ", self.Start_Time, " [-]")
            print ("{FENICS} Start time: ", self.Time_Steps, " [s]")
            print ("{FENICS} Numbers of sub-iterations: ", self.num_sub_iteration(), " [-]")
            print ("\n")

    def print_Disp (self, msh, displacement_function):
        # Compute and print the displacement of monitored point
        # Create bounding box for function evaluation
        bb_tree = geometry.bb_tree(msh, 2)

        # Check against standard table value
        p = np.array([self.pointMoniX(),self.pointMoniY(),self.pointMoniZ()], dtype=np.float64)
        cell_candidates = geometry.compute_collisions_points(bb_tree, p)
        cells = geometry.compute_colliding_cells(msh, cell_candidates, p).array
        displacement_function.x.scatter_forward()
        if len(cells) > 0:
            d_DispSum = np.zeros(3)
            d_tempDenominator  = np.array([ self.size,
                                            self.size,
                                            self.size])
            self.LOCAL_COMM_WORLD.Reduce((displacement_function.eval(p, cells[0])),
                                        d_DispSum,op=MPI.SUM,root=0)
            d_Disp = np.divide(d_DispSum,d_tempDenominator)

            if self.rank == 0:
                print ("{FENICS} Monitored point deflection [m]: ", d_Disp)

    def Export_Disp_txt(self, msh, displacement_function):
        if self.iExporttxt():
            # Create bounding box for function evaluation
            bb_tree = geometry.bb_tree(msh, 2)

            # Check against standard table value
            p = np.array([self.pointMoniX(),self.pointMoniY(),self.pointMoniZ()], dtype=np.float64)
            pb = np.array([self.pointMoniXb(),self.pointMoniYb(),self.pointMoniZb()], dtype=np.float64)
            cell_candidates = geometry.compute_collisions_points(bb_tree, p)
            cells = geometry.compute_colliding_cells(msh, cell_candidates, p).array
            displacement_function.x.scatter_forward()
            if len(cells) > 0:
                pointMoniDispSum = np.zeros(3)
                pointMoniDispSumb = np.zeros(3)
                tempDenominator  = np.array([ self.size,
                                                self.size,
                                                self.size])
                self.LOCAL_COMM_WORLD.Reduce((displacement_function.eval(p, cells[0])),
                                            pointMoniDispSum,op=MPI.SUM,root=0)
                self.LOCAL_COMM_WORLD.Reduce((displacement_function.eval(pb, cells[0])),
                                            pointMoniDispSumb,op=MPI.SUM,root=0)
                pointMoniDisp = np.divide(pointMoniDispSum,tempDenominator)
                pointMoniDispb = np.divide(pointMoniDispSumb,tempDenominator)

                for irank in range(self.size):
                    if self.rank == irank:
                        ftxt_dispX = open(self.outputFolderPath + "/tip-displacementX_" + str(irank)+ ".txt", "a")
                        ftxt_dispX.write(str(pointMoniDisp[0]))
                        ftxt_dispX.write("\n")
                        ftxt_dispX.close

                        ftxt_dispY = open(self.outputFolderPath + "/tip-displacementY_" + str(irank)+ ".txt", "a")
                        ftxt_dispY.write(str(pointMoniDisp[1]))
                        ftxt_dispY.write("\n")
                        ftxt_dispY.close

                        ftxt_dispZ = open(self.outputFolderPath + "/tip-displacementZ_" + str(irank)+ ".txt", "a")
                        ftxt_dispZ.write(str(pointMoniDisp[2]))
                        ftxt_dispZ.write("\n")
                        ftxt_dispZ.close

                        ftxt_dispXb = open(self.outputFolderPath + "/tip-displacementXb_" + str(irank)+ ".txt", "a")
                        ftxt_dispXb.write(str(pointMoniDispb[0]))
                        ftxt_dispXb.write("\n")
                        ftxt_dispXb.close

                        ftxt_dispYb = open(self.outputFolderPath + "/tip-displacementYb_" + str(irank)+ ".txt", "a")
                        ftxt_dispYb.write(str(pointMoniDispb[1]))
                        ftxt_dispYb.write("\n")
                        ftxt_dispYb.close

                        ftxt_dispZb = open(self.outputFolderPath + "/tip-displacementZb_" + str(irank)+ ".txt", "a")
                        ftxt_dispZb.write(str(pointMoniDispb[2]))
                        ftxt_dispZb.write("\n")
                        ftxt_dispZb.close

    def Export_Disp_xdmf(self,
                        Current_Time_Step,
                        current_time,
                        mesh,
                        grid_dimension,
                        VectorFunctionSpace,
                        VectorFunctionSpace1,
                        displacement_function):
        # Export post-processing files
        if ((self.rank == 0) and self.iDebug()):
            print ("\n")
            print ("{FENICS} time steps: ", Current_Time_Step,
                    " output_interval: ", self.output_interval(),
                    " %: ", (Current_Time_Step % self.output_interval()))

        if (Current_Time_Step % self.output_interval()) == 0:
            if self.rank == 0:
                print ("\n")
                print ("{FENICS} Export files at ", current_time, " [s] ...   ", end="", flush=True)

            # # Compute stress
            # Vsig = fem.FunctionSpace(mesh, ("Lagrange", self.deg_fun_spc()))
            # sig = fem.Function(Vsig, name="Stress")
            # sigma_dev = self.sigma(displacement_function,grid_dimension) - (1 / 3) * ufl.tr(self.sigma(displacement_function, grid_dimension)) * ufl.Identity(len(displacement_function))
            # sigma_vm = ufl.sqrt((3 / 2) * inner(sigma_dev, sigma_dev))
            # stress_expr = fem.Expression(sigma_vm, Vsig.element.interpolation_points())
            # sig.interpolate(stress_expr)
            #
            # # Save stress solution to file
            # self.stress_file.write_function(sig, current_time)

            # Save displacement solution to file
            displacement_function1 = fem.Function(VectorFunctionSpace1)              # Function for displacement by MCK solving method
            displacement_function1.interpolate(displacement_function)
            self.disp_file.write_function(displacement_function1, current_time)

            # Compute traction
            traction = fem.Function(VectorFunctionSpace1, name="Traction")
            traction.interpolate(self.tF_apply)
            # Save traction solution to file
            self.traction_file.write_function(traction, current_time)
            if self.rank == 0: print ("Done")
        else:
            pass

    def Post_Solving_Log(self, simtime):
        if self.rank == 0:
            print ("\n")
            print ("{FENICS} Current Date and Time: ", datetime.datetime.now())
            print ("\n")
            print ("{FENICS} Total Simulation time: %g [s]" % simtime)
            print ("\n")
            print ("{FENICS} ********** STRUCTURAL-ELASTICITY SIMULATION COMPLETED **********")

    def Create_Post_Process_Files(self, mesh):
        if self.rank == 0: print ("{FENICS} Preparing post-process files ...   ", end="", flush=True)
        self.disp_file=XDMFFile(self.LOCAL_COMM_WORLD, self.outputFolderPath + "/displacement.xdmf", "w")
        self.disp_file.write_mesh(mesh)
        self.stress_file=XDMFFile(self.LOCAL_COMM_WORLD, self.outputFolderPath + "/stress.xdmf", "w")
        self.stress_file.write_mesh(mesh)
        self.traction_file=XDMFFile(self.LOCAL_COMM_WORLD, self.outputFolderPath + "/surface_traction_structure.xdmf", "w")
        self.traction_file.write_mesh(mesh)
        if self.rank == 0: print ("Done")

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #%% Setup checkpoint file
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def Checkpoint_Output_Linear(self,
                                 current_time,
                                 mesh,
                                 d0mck_Functions_previous,
                                 u0mck_Functions_previous,
                                 a_Function_previous,
                                 dmck_Function,
                                 File_Exists=True):
        if File_Exists:
            import os
            if self.rank == 0:
                os.remove(self.outputFolderPath + "/checkpointData_" + str(current_time) +".h5")
            self.LOCAL_COMM_WORLD.Barrier()
        else:
            pass

        hdf5checkpointDataOut = HDF5File(self.LOCAL_COMM_WORLD, self.outputFolderPath + "/checkpointData_" + str(current_time) +".h5", "w")
        hdf5checkpointDataOut.write(mesh, "/mesh")
        hdf5checkpointDataOut.write(u0mck_Functions_previous, "/u0mck", current_time)
        hdf5checkpointDataOut.write(d0mck_Functions_previous, "/d0mck", current_time)
        hdf5checkpointDataOut.write(a_Function_previous, "/a0mck", current_time)
        hdf5checkpointDataOut.write(dmck_Function, "/dmck", current_time)
        hdf5checkpointDataOut.write(self.areaf, "/areaf")
        hdf5checkpointDataOut.close()
        # Delete HDF5File object, closing file
        del hdf5checkpointDataOut

    def Load_Functions_Continue_Run_Linear(self,
                                           d0mck,
                                           u0mck,
                                           a0mck,
                                           dmck):
        if self.iContinueRun():
            hdf5checkpointDataInTemp = HDF5File(self.LOCAL_COMM_WORLD, self.inputFolderPath + "/checkpointData.h5", "r")
            hdf5checkpointDataInTemp.read(d0mck, "/d0mck/vector_0")
            hdf5checkpointDataInTemp.read(u0mck, "/u0mck/vector_0")
            hdf5checkpointDataInTemp.read(a0mck, "/a0mck/vector_0")
            hdf5checkpointDataInTemp.read(dmck, "/dmck/vector_0")
            hdf5checkpointDataInTemp.close()
            # Delete HDF5File object, closing file
            del hdf5checkpointDataInTemp
        else:
            pass

    def Checkpoint_Output_Nonlinear(self,
                                    current_time,
                                    mesh,
                                    ud_Functions_previous,
                                    ud_Functions,
                                    File_Exists=True):
        if File_Exists:
            import os
            if self.rank == 0:
                os.remove(self.outputFolderPath + "/checkpointData_" + str(current_time) +".h5")
            self.LOCAL_COMM_WORLD.Barrier()
        else:
            pass

        hdf5checkpointDataOut = HDF5File(self.LOCAL_COMM_WORLD, self.outputFolderPath + "/checkpointData_" + str(current_time) +".h5", "w")
        hdf5checkpointDataOut.write(mesh, "/mesh")
        hdf5checkpointDataOut.write(ud_Functions_previous, "/u0d0", current_time)
        hdf5checkpointDataOut.write(ud_Functions, "/ud", current_time)
        hdf5checkpointDataOut.write(self.areaf, "/areaf")
        hdf5checkpointDataOut.close()
        # Delete HDF5File object, closing file
        del hdf5checkpointDataOut

    def Load_Functions_Continue_Run_Nonlinear(self, u0d0, ud, sigma_s):
        if self.iContinueRun():
            hdf5checkpointDataInTemp = HDF5File(self.LOCAL_COMM_WORLD, self.inputFolderPath + "/checkpointData.h5", "r")
            hdf5checkpointDataInTemp.read(u0d0, "/u0d0/vector_0")
            hdf5checkpointDataInTemp.read(ud, "/ud/vector_0")
            hdf5checkpointDataInTemp.read(sigma_s, "/sigma_s/vector_0")
            hdf5checkpointDataInTemp.close()
            # Delete HDF5File object, closing file
            del hdf5checkpointDataInTemp
        else:
            pass

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FILE END  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#