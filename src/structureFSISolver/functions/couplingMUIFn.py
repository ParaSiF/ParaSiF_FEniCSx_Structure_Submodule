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

    @file couplingMUIFn.py
    
    @author W. Liu
    
    @brief MUI related coupling functions.

"""

#_________________________________________________________________________________________
#
#%% Import packages
#_________________________________________________________________________________________
import sys
from mpi4py import MPI
import mui4py
import petsc4py
import os

class couplingMUIFn:

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #%% Initialise MPI by mpi4py/MUI for parallelised computation
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def MUI_Init(self, arg_COMM_WORLD=None, arg_ifaces3d=None):
        if self.iMUICoupling():
            # App common world claims
            if arg_COMM_WORLD is not None:
                self.LOCAL_COMM_WORLD = arg_COMM_WORLD
            else:
                self.LOCAL_COMM_WORLD = mui4py.mpi_split_by_app()
            # MUI parameters
            if arg_ifaces3d is not None:
                self.ifaces3d = arg_ifaces3d
            else:
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
                self.ifaces3d = mui4py.create_unifaces(domain, iface, config3d)
                self.ifaces3d["threeDInterface0"].set_data_types(data_types)

                # URI = "mpi://structureDomain/threeDInterface0" 
                # self.iface = mui4py.Uniface(uri=URI, config=config3d) 
                # self.iface.set_data_types(data_types)
                # App common world claims
                # self.LOCAL_COMM_WORLD = mui4py.mpi_split_by_app()

                # Necessary to avoid hangs at PETSc vector communication
                petsc4py.init(comm=self.LOCAL_COMM_WORLD)
        else:
            # App common world claims
            self.LOCAL_COMM_WORLD = MPI.COMM_WORLD

        # Define local communicator rank
        self.rank = self.LOCAL_COMM_WORLD.Get_rank()

        # Define local communicator size
        self.size = self.LOCAL_COMM_WORLD.Get_size()

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #%% Define MUI samplers and commit ZERO step
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def MUI_Sampler_Define(self,
                           function_space,
                           grid_dimension,
                           dofs_fetch_list,
                           dofs_push_list,
                           xyz_fetch,
                           Total_Time_Steps):

        if self.iMUICoupling():
            synchronised=False

            dofs_to_xyz = self.dofs_to_xyz(function_space, grid_dimension)

            send_min_X = sys.float_info.max
            send_min_Y = sys.float_info.max
            send_min_Z = sys.float_info.max

            send_max_X = -sys.float_info.max
            send_max_Y = -sys.float_info.max
            send_max_Z = -sys.float_info.max

            for i, p in enumerate(dofs_push_list):
                if (dofs_to_xyz[p][0] < send_min_X):
                    send_min_X = dofs_to_xyz[p][0]

                if (dofs_to_xyz[p][1] < send_min_Y):
                    send_min_Y = dofs_to_xyz[p][1]

                if (dofs_to_xyz[p][2] < send_min_Z):
                    send_min_Z = dofs_to_xyz[p][2]

                if (dofs_to_xyz[p][0] > send_max_X):
                    send_max_X = dofs_to_xyz[p][0]

                if (dofs_to_xyz[p][1] > send_max_Y):
                    send_max_Y = dofs_to_xyz[p][1]

                if (dofs_to_xyz[p][2] > send_max_Z):
                    send_max_Z = dofs_to_xyz[p][2]

            if (send_max_X < send_min_X):
                print("{** FENICS ERROR **} send_max_X: ", send_max_X, " smaller than send_min_X: ", send_min_X, " at rank: ", self.rank)

            if (send_max_Y < send_min_Y):
                print("{** FENICS ERROR **} send_max_Y: ", send_max_Y, " smaller than send_min_Y: ", send_min_Y, " at rank: ", self.rank)

            if (send_max_Z < send_min_Z):
                print("{** FENICS ERROR **} send_max_Z: ", send_max_Z, " smaller than send_min_Z: ", send_min_Z, " at rank: ", self.rank)

            if (len(dofs_push_list)!=0):
                # Set up sending span
                span_push = mui4py.geometry.Box([send_min_X, send_min_Y, send_min_Z],
                                                [send_max_X, send_max_Y, send_max_Z])

                # Announce the MUI send span
                self.ifaces3d.announce_send_span(0, Total_Time_Steps*self.num_sub_iteration(), span_push, synchronised)
                # self.iface.announce_send_span(0, Total_Time_Steps*self.num_sub_iteration(), span_push, synchronised)

                print("{FENICS} at rank: ", self.rank, " send_max_X: ", send_max_X, " send_min_X: ", send_min_X)
                print("{FENICS} at rank: ", self.rank, " send_max_Y: ", send_max_Y, " send_min_Y: ", send_min_Y)
                print("{FENICS} at rank: ", self.rank, " send_max_Z: ", send_max_Z, " send_min_Z: ", send_min_Z)

            else:
                # Announce the MUI send span
                #self.ifaces3d.announce_send_disable()
                pass

            recv_min_X = sys.float_info.max
            recv_min_Y = sys.float_info.max
            recv_min_Z = sys.float_info.max

            recv_max_X = -sys.float_info.max
            recv_max_Y = -sys.float_info.max
            recv_max_Z = -sys.float_info.max

            # Declare list to store mui::point3d
            point3dList = []
            point3dGlobalID = []

            for i, p in enumerate(dofs_fetch_list):
                if (dofs_to_xyz[p][0] < recv_min_X):
                    recv_min_X = dofs_to_xyz[p][0]

                if (dofs_to_xyz[p][1] < recv_min_Y):
                    recv_min_Y = dofs_to_xyz[p][1]

                if (dofs_to_xyz[p][2] < recv_min_Z):
                    recv_min_Z = dofs_to_xyz[p][2]

                if (dofs_to_xyz[p][0] > recv_max_X):
                    recv_max_X = dofs_to_xyz[p][0]

                if (dofs_to_xyz[p][1] > recv_max_Y):
                    recv_max_Y = dofs_to_xyz[p][1]

                if (dofs_to_xyz[p][2] > recv_max_Z):
                    recv_max_Z = dofs_to_xyz[p][2]

                point_fetch = self.ifaces3d.Point([dofs_to_xyz[p][0],
                                                                        dofs_to_xyz[p][1],
                                                                        dofs_to_xyz[p][2]])
                # point_fetch = self.iface.Point([dofs_to_xyz[p][0],
                #                             dofs_to_xyz[p][1],
                #                             dofs_to_xyz[p][2]])

                point_ID = -999
                for ii, pp in enumerate(xyz_fetch):
                    if (pp[0] == dofs_to_xyz[p][0]):
                        if (pp[1] == dofs_to_xyz[p][1]):
                            if (pp[2] == dofs_to_xyz[p][2]):
                                point_ID = ii
                                break

                if (point_ID<0):
                    print("{** FENICS ERROR **} cannot find point: ", dofs_to_xyz[p][0],
                                                                        dofs_to_xyz[p][1],
                                                                        dofs_to_xyz[p][2],
                                                                        " in Global xyz fetch list")
                point3dList.append(point_fetch)
                point3dGlobalID.append(point_ID)
# Debugging start
            #print("point3dList:")
            #for point in point3dList:
            #    print(point)

            #print("\npoint3dGlobalID:")
            #for point_ID in point3dGlobalID:
            #    print(point_ID)
# Debugging end
            if (recv_max_X < recv_min_X):
                print("{** FENICS ERROR **} recv_max_X: ", recv_max_X, " smaller than recv_min_X: ", recv_min_X, " at rank: ", self.rank)

            if (recv_max_Y < recv_min_Y):
                print("{** FENICS ERROR **} recv_max_Y: ", recv_max_Y, " smaller than recv_min_Y: ", recv_min_Y, " at rank: ", self.rank)

            if (recv_max_Z < recv_min_Z):
                print("{** FENICS ERROR **} recv_max_Z: ", recv_max_Z, " smaller than recv_min_Z: ", recv_min_Z, " at rank: ", self.rank)

            if (len(dofs_fetch_list)!=0):
                # Set up receiving span
                span_fetch = mui4py.geometry.Box([recv_min_X, recv_min_Y, recv_min_Z],
                                                 [recv_max_X, recv_max_Y, recv_max_Z])

                # Announce the MUI receive span
                self.ifaces3d.announce_recv_span(0, Total_Time_Steps*self.num_sub_iteration()*10, span_fetch, synchronised)
                # self.iface.announce_recv_span(0, Total_Time_Steps*self.num_sub_iteration()*10, span_fetch, synchronised)

                print("{FENICS} at rank: ", self.rank, " recv_max_X: ", recv_max_X, " recv_min_X: ", recv_min_X)
                print("{FENICS} at rank: ", self.rank, " recv_max_Y: ", recv_max_Y, " recv_min_Y: ", recv_min_Y)
                print("{FENICS} at rank: ", self.rank, " recv_max_Z: ", recv_max_Z, " recv_min_Z: ", recv_min_Z)

            else:
                # Announce the MUI receive span
                #self.ifaces3d.announce_recv_disable()
                pass

            # Spatial/temporal samplers
            if self.rank == 0: print ("{FENICS} Defining MUI samplers ...   ", end="", flush=True)

            fileAddressLocalMake = self.outputFolderName() + '/RBFMatrix'
            fileAddress=fileAddressLocalMake + '/' + str(self.rank)
            os.makedirs(fileAddressLocalMake)

            if (self.iWriteMatrix()):
                if self.rank == 0:
                    with open(fileAddressLocalMake+'/partitionSize.dat', 'w') as f_ps:
                        f_ps.write("%i\n" % self.size)
            else:    
                print ("{FENICS} Reading RBF matrix from ", self.rank)
                sourcefileAddress=self.inputFolderName() + '/RBFMatrix'

                # search line number of the pointID
                numberOfFolders = 0
                with open(sourcefileAddress +'/partitionSize.dat', 'r') as f_psr:
                    print ("{FENICS} open partitionSize from ", self.rank)
                    for line in f_psr:
                        numberOfFolders = int(line)
                f_psr.close()
                print ("{FENICS} Number of RBF subfolders: ", numberOfFolders, " from ", self.rank)

                numberOfCols=-99
                for i, point_IDs in enumerate(point3dGlobalID):
                    # search line number of the pointID
                    iFolder=0
                    while iFolder < numberOfFolders:
                        line_number = -1
                        result_line_number = -99
                        result_folder_number = -99
                        with open(sourcefileAddress+'/'+str(iFolder)+'/pointID.dat', 'r') as f_pid:
                            for line in f_pid:
                                line_number += 1
                                if str(point_IDs) in line:
                                    result_line_number = line_number
                                    result_folder_number = iFolder
                                    break
                        f_pid.close()
                        iFolder += 1
                        if (result_folder_number >= 0):
                            break

                    if (result_line_number < 0):
                        print ("{** FENICS ERROR **} Cannot find Point ID: ", point_ID, " in file")
                    # Get the line in H matrix and copy to local file
                    with open(sourcefileAddress+'/'+str(result_folder_number)+'/Hmatrix.dat', 'r') as f_h:
                        for i, line in enumerate(f_h):
                            if i == (result_line_number+6):
                                with open(fileAddress+'/Hmatrix.dat', 'a') as f_h_result:
                                    if line[-1] == '\n':
                                        f_h_result.write(line)
                                    else:
                                        f_h_result.write(line+'\n')
                                    if (numberOfCols<0):
                                        numberOfCols=len(line.split(","))
                                f_h_result.close()
                            elif i > (result_line_number+6):
                                break
                    f_h.close()

                with open(fileAddress+'/matrixSize.dat', 'w') as f_size:
                    f_size.write(str(numberOfCols)+","+str(len(point3dGlobalID))+",0,0,"+str(len(point3dGlobalID))+","+str(numberOfCols))
                f_size.close()

            # Best practice suggestion: for a better performance on the RBF method, always switch on the smoothFunc when structure Dofs are more than
            #                           fluid points; Tune the rMUIFetcher to receive a reasonable totForce_Fetch value; Tune the areaListFactor to
            #                           ensure totForce_Fetch and Total_Force_on_structure are the same.
            self.t_sampler = mui4py.TemporalSamplerExact()
# Problem creating RBF matrix folder here!!
            #self.s_sampler = mui4py.SamplerPseudoNearestNeighbor(self.rMUIFetcher())
            self.s_sampler = mui4py.SamplerRbf(self.rMUIFetcher(),
                                               point3dList,
                                               self.basisFunc(),
                                               self.iConservative(),
                                               self.iSmoothFunc(),
                                               self.iWriteMatrix(),
                                               '',
                                               self.cutoffRBF(),
                                               self.cgSolveTolRBF(),
                                               self.cgMaxIterRBF(),
                                               self.pouSizeRBF(),
                                               self.precondRBF(),
                                               self.LOCAL_COMM_WORLD)

            with open(fileAddressLocalMake+'/pointID.dat', 'w') as f_pid:
                for pid in point3dGlobalID:
                    f_pid.write("%i\n" % pid)

            # Commit ZERO step
            self.ifaces3d.commit(0)
            # self.iface.commit(0)
            if self.rank == 0: print ("{FENICS} Commit ZERO step")
        else:
            pass

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #%% Define MUI Fetch and Push
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def MUI_Fetch(self, dofs_to_xyz, dofs_fetch_list, total_Sub_Iteration):
        totForceX = 0.0
        totForceY = 0.0
        totForceZ = 0.0
        temp_vec_function_temp = self.tF_apply_vec

        if self.iparallelFSICoupling():
            fetch_iteration = total_Sub_Iteration-1
        else:
            fetch_iteration = total_Sub_Iteration

        if (fetch_iteration >= 0):
            if self.iMUIFetchMany():
                temp_vec_function_temp[0::3][dofs_fetch_list] = self.ifaces3d.\
                            fetch_many("forceX",
                                       dofs_to_xyz,
                                       fetch_iteration,
                                       self.s_sampler,
                                       self.t_sampler)
                temp_vec_function_temp[1::3][dofs_fetch_list] = self.ifaces3d.\
                            fetch_many("forceY",
                                       dofs_to_xyz,
                                       fetch_iteration,
                                       self.s_sampler,
                                       self.t_sampler)
                temp_vec_function_temp[2::3][dofs_fetch_list] = self.ifaces3d.\
                            fetch_many("forceZ",
                                       dofs_to_xyz,
                                       fetch_iteration,
                                       self.s_sampler,
                                       self.t_sampler)
                # temp_vec_function_temp[0::3][dofs_fetch_list] = self.iface.\
                #             fetch_many("forceX",
                #                        dofs_to_xyz,
                #                        fetch_iteration,
                #                        self.s_sampler,
                #                        self.t_sampler)
                # temp_vec_function_temp[1::3][dofs_fetch_list] = self.iface.\
                #             fetch_many("forceY",
                #                        dofs_to_xyz,
                #                        fetch_iteration,
                #                        self.s_sampler,
                #                        self.t_sampler)
                # temp_vec_function_temp[2::3][dofs_fetch_list] = self.iface.\
                #             fetch_many("forceZ",
                #                        dofs_to_xyz,
                #                        fetch_iteration,
                #                        self.s_sampler,
                #                        self.t_sampler)

                for i, p in enumerate(dofs_fetch_list):
                    if self.iparallelFSICoupling():
                        self.tF_apply_vec[0::3][p] += (temp_vec_function_temp[0::3][p] - \
                                                       self.tF_apply_vec[0::3][p])*self.undRelxCpl()
                        self.tF_apply_vec[1::3][p] += (temp_vec_function_temp[1::3][p] - \
                                                       self.tF_apply_vec[1::3][p])*self.undRelxCpl()
                        self.tF_apply_vec[2::3][p] += (temp_vec_function_temp[2::3][p] - \
                                                       self.tF_apply_vec[2::3][p])*self.undRelxCpl()
                    else:
                        self.tF_apply_vec[0::3][p] = temp_vec_function_temp[0::3][p]
                        self.tF_apply_vec[1::3][p] = temp_vec_function_temp[1::3][p]
                        self.tF_apply_vec[2::3][p] = temp_vec_function_temp[2::3][p]

                    totForceX += self.tF_apply_vec[0::3][p]
                    totForceY += self.tF_apply_vec[1::3][p]
                    totForceZ += self.tF_apply_vec[2::3][p]

                    if (self.areaf_vec[p] == 0):
                        self.tF_apply_vec[0::3][p] = 0.
                        self.tF_apply_vec[1::3][p] = 0.
                        self.tF_apply_vec[2::3][p] = 0.
                    else:
                        self.tF_apply_vec[0::3][p] /= self.areaf_vec[p]
                        self.tF_apply_vec[1::3][p] /= self.areaf_vec[p]
                        self.tF_apply_vec[2::3][p] /= self.areaf_vec[p]

            else:
                if (fetch_iteration >= 0):
                    for i, p in enumerate(dofs_fetch_list):
                        temp_vec_function_temp[0::3][p] = self.ifaces3d.\
                                    fetch("forceX",
                                          dofs_to_xyz[i],
                                          fetch_iteration,
                                          self.s_sampler,
                                          self.t_sampler)

                        temp_vec_function_temp[1::3][p] = self.ifaces3d.\
                                    fetch("forceY",
                                          dofs_to_xyz[i],
                                          fetch_iteration,
                                          self.s_sampler,
                                          self.t_sampler)

                        temp_vec_function_temp[2::3][p] = self.ifaces3d.\
                                    fetch("forceZ",
                                          dofs_to_xyz[i],
                                          fetch_iteration,
                                          self.s_sampler,
                                          self.t_sampler)
                        # temp_vec_function_temp[0::3][p] = self.iface.\
                        #             fetch("forceX",
                        #                   dofs_to_xyz[i],
                        #                   fetch_iteration,
                        #                   self.s_sampler,
                        #                   self.t_sampler)

                        # temp_vec_function_temp[1::3][p] = self.iface.\
                        #             fetch("forceY",
                        #                   dofs_to_xyz[i],
                        #                   fetch_iteration,
                        #                   self.s_sampler,
                        #                   self.t_sampler)

                        # temp_vec_function_temp[2::3][p] = self.iface.\
                        #             fetch("forceZ",
                        #                   dofs_to_xyz[i],
                        #                   fetch_iteration,
                        #                   self.s_sampler,
                        #                   self.t_sampler)
                        # Debugging lines begin
                        #if((total_Sub_Iteration * 0.1) <= 7.0):
                        #    temp_vec_function_temp[1::3][p] = -((total_Sub_Iteration * 0.1)/7.0) * 20
                        #else:
                        #    temp_vec_function_temp[1::3][p] = 0.0
                        # Debugging lines end
                        if self.iparallelFSICoupling():
                            self.tF_apply_vec[0::3][p] += (temp_vec_function_temp[0::3][p] - \
                                                           self.tF_apply_vec[0::3][p])*self.undRelxCpl()
                            self.tF_apply_vec[1::3][p] += (temp_vec_function_temp[1::3][p] - \
                                                           self.tF_apply_vec[1::3][p])*self.undRelxCpl()
                            self.tF_apply_vec[2::3][p] += (temp_vec_function_temp[2::3][p] - \
                                                           self.tF_apply_vec[2::3][p])*self.undRelxCpl()
                        else:
                            self.tF_apply_vec[0::3][p] = temp_vec_function_temp[0::3][p]
                            self.tF_apply_vec[1::3][p] = temp_vec_function_temp[1::3][p]
                            self.tF_apply_vec[2::3][p] = temp_vec_function_temp[2::3][p]

                        totForceX += self.tF_apply_vec[0::3][p]
                        totForceY += self.tF_apply_vec[1::3][p]
                        totForceZ += self.tF_apply_vec[2::3][p]

                        self.tF_apply_vec[0::3][p] /= self.areaf_vec[p]
                        self.tF_apply_vec[1::3][p] /= self.areaf_vec[p]
                        self.tF_apply_vec[2::3][p] /= self.areaf_vec[p]

                    if self.iDebug():
                        print ("{FENICS**} totForce Apply: ", totForceX, "; ",totForceY, "; ",totForceZ,
                                "; at iteration: ", fetch_iteration, " at rank: ", self.rank)

    def MUI_Push(self, dofs_to_xyz, dofs_push, displacement_function, total_Sub_Iteration):
        d_vec_x = displacement_function.x.array[0::3]
        d_vec_y = displacement_function.x.array[1::3]
        d_vec_z = displacement_function.x.array[2::3]

        if self.iMUIPushMany():
            if self.iPushX():
                self.ifaces3d.\
                            push_many("dispX", dofs_to_xyz, (d_vec_x[dofs_push]))
            if self.iPushY():
                self.ifaces3d.\
                            push_many("dispY", dofs_to_xyz, (d_vec_y[dofs_push]))
            if self.iPushZ():
                self.ifaces3d.\
                            push_many("dispZ", dofs_to_xyz, (d_vec_z[dofs_push]))

            a = self.ifaces3d.\
                            commit(total_Sub_Iteration)
            #     self.iface.\
            #                 push_many("dispX", dofs_to_xyz, (d_vec_x[dofs_push]))
            # if self.iPushY():
            #     self.iface.\
            #                 push_many("dispY", dofs_to_xyz, (d_vec_y[dofs_push]))
            # if self.iPushZ():
            #     self.iface.\
            #                 push_many("dispZ", dofs_to_xyz, (d_vec_z[dofs_push]))

            # a = self.iface.\
            #                 commit(total_Sub_Iteration)
        else:
            if self.iPushX():
                for i, p in enumerate(dofs_push):
                    self.ifaces3d.\
                            push("dispX", dofs_to_xyz[i], (d_vec_x[p]))
            if self.iPushY():
                for i, p in enumerate(dofs_push):
                    self.ifaces3d.\
                            push("dispY", dofs_to_xyz[i], (d_vec_y[p]))
            if self.iPushZ():
                for i, p in enumerate(dofs_push):
                    self.ifaces3d.\
                            push("dispZ", dofs_to_xyz[i], (d_vec_z[p]))

            a = self.ifaces3d.\
                            commit(total_Sub_Iteration)
            # if self.iPushX():
            #     for i, p in enumerate(dofs_push):
            #         self.iface.\
            #                 push("dispX", dofs_to_xyz[i], (d_vec_x[p]))
            # if self.iPushY():
            #     for i, p in enumerate(dofs_push):
            #         self.iface.\
            #                 push("dispY", dofs_to_xyz[i], (d_vec_y[p]))
            # if self.iPushZ():
            #     for i, p in enumerate(dofs_push):
            #         self.iface.\
            #                 push("dispZ", dofs_to_xyz[i], (d_vec_z[p]))

            # a = self.iface.\
            #                 commit(total_Sub_Iteration)

        if (self.rank == 0) and self.iDebug():
            print ('{FENICS} MUI commit step: ',total_Sub_Iteration)

        if ((total_Sub_Iteration-self.forgetTStepsMUI()) > 0):
            a = self.ifaces3d.\
                            forget(total_Sub_Iteration-self.forgetTStepsMUI())
            self.ifaces3d.\
                            set_memory(self.forgetTStepsMUI())
            # a = self.iface.\
            #                 forget(total_Sub_Iteration-self.forgetTStepsMUI())
            # self.iface.\
            #                 set_memory(self.forgetTStepsMUI())
            if (self.rank == 0) and self.iDebug():
                print ('{FENICS} MUI forget step: ',(total_Sub_Iteration-self.forgetTStepsMUI()))

    def MUI_Commit_only(self, total_Sub_Iteration):
        a = self.ifaces3d.\
                            commit(total_Sub_Iteration)
        # a = self.iface.\
        #                     commit(total_Sub_Iteration)

        if (self.rank == 0) and self.iDebug():
            print ('{FENICS} MUI commit step: ',total_Sub_Iteration)

        if ((total_Sub_Iteration-self.forgetTStepsMUI()) > 0):
            a = self.ifaces3d.\
                            forget(total_Sub_Iteration-self.forgetTStepsMUI())
            self.ifaces3d.\
                            set_memory(self.forgetTStepsMUI())
            # a = self.iface.\
            #                 forget(total_Sub_Iteration-self.forgetTStepsMUI())
            # self.iface.\
            #                 set_memory(self.forgetTStepsMUI())
            if (self.rank == 0) and self.iDebug():
                print ('{FENICS} MUI forget step: ',(total_Sub_Iteration-self.forgetTStepsMUI()))

