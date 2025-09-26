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

    @file linearElastic.py

    @author W. Liu

    @brief The core linear solver class of the structure code.

"""

#_________________________________________________________________________________________
#
#%% Import packages
#_________________________________________________________________________________________

from dolfinx import *
from dolfinx.fem.petsc import (apply_lifting, assemble_matrix, assemble_vector,
                               set_bc)
import os
import numpy as np
from mpi4py import MPI
from petsc4py import PETSc
from petsc4py.PETSc import ScalarType
import ufl
import structureFSISolver

class linearElastic:

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #%% Main solver function
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def linearElasticSolve(self, iface3d=None):

        #===========================================
        #%% Time marching parameters define
        #===========================================

        t        = self.Start_Time
        t_step   = self.Time_Steps
        i_sub_it = self.Start_Number_Sub_Iteration

        #===========================================
        #%% Solid Mesh input/generation
        #===========================================

        domain = self.Mesh_Generation()
        gdim   = self.Get_Grid_Dimension(domain)

        #===========================================
        #%% Define coefficients
        #===========================================

        # Time lists
        times    = []
        t_sub_it = 0

        # Rayleigh damping coefficients
        alpha_rdc = fem.Constant(domain, self.alpha_rdc())
        beta_rdc  = fem.Constant(domain, self.beta_rdc())

        # Generalized-alpha method parameters
        # alpha_m_gam >= alpha_f_gam >= 0.5 for a better performance
        alpha_m_gam = fem.Constant(domain, self.alpha_m_gam())
        alpha_f_gam = fem.Constant(domain, self.alpha_f_gam())
        gamma_gam   = (1./2.) + alpha_m_gam - alpha_f_gam
        beta_gam    = (1./4.) * (gamma_gam + (1./2.))**2

        #===========================================
        #%% Define function spaces
        #===========================================

        if self.rank == 0: print ("{FENICS} Creating function spaces ...   ")
        Q         =     fem.functionspace(domain, ("Lagrange", self.deg_fun_spc()))            # Function space with updated mesh
        V         =     fem.functionspace(domain, ("Lagrange", self.deg_fun_spc(), (domain.geometry.dim, )))
        V1        =     fem.functionspace(domain, ("Lagrange", 1, (domain.geometry.dim, )))

        if self.rank == 0: print ("{FENICS} Done with creating function spaces")

        #=======================================================
        #%% Define functions, test functions and trail functions
        #=======================================================

        if self.rank == 0: print ("{FENICS} Creating functions, test functions and trail functions ...   ", end="", flush=True)

        # Trial functions
        ddmck  = ufl.TrialFunction(V)        # Trial function for displacement by MCK solving method

        # Test functions
        chi      = ufl.TestFunction(V)           # Test function for displacement by MCK solving method

        # Functions at present time step
        dmck = fem.Function(V)              # Function for displacement by MCK solving method

        # Functions at previous time step
        d0mck  = fem.Function(V)             # Function for displacement by MCK solving method
        u0mck  = fem.Function(V)             # Function for velocity by MCK solving method
        a0mck  = fem.Function(V)             # Function for acceleration by MCK solving method

        self.Load_Functions_Continue_Run_Linear(d0mck,u0mck,a0mck,dmck)

        if self.rank == 0: print ("Done")

        #===========================================
        #%% Define traction forces
        #===========================================

        self.Traction_Define(V)

        #===========================================
        #%% Create facet to cell connectivity
        #===========================================

        tdim = domain.topology.dim
        fdim = tdim - 1
        domain.topology.create_connectivity(fdim, tdim)
        boundary_facets = mesh.exterior_facet_indices(domain.topology)

        #===========================================
        #%% Define SubDomains and boundaries
        #===========================================

        self.Boundaries_Generation_Fixed_Flex_Sym(domain, V)

        ds = self.Get_ds(domain)

        #===========================================
        #%% Define boundary conditions
        #===========================================

        if self.rank == 0: print ("{FENICS} Creating 3D boundary conditions ...   ", end="", flush=True)
        bc1 = self.dirichletBCs.DirichletBCs(V,self.fixeddofs)
        bcs = [bc1]
        if self.rank == 0: print ("Done")

        #===========================================
        #%% Define DOFs and Coordinates mapping
        #===========================================  
        #  !! OUTDATED FUNCTION, NEED UPDATED TO FENICS-X !!
        dofs_fetch_list = self.dofs_list(domain, V, 2)

        dofs_push_list = self.dofs_list(domain, V, 2)

        xyz_fetch = self.xyz_np(self.flexdofs, V, gdim)

        xyz_push = self.xyz_np(self.flexdofs, V, gdim)

        #===========================================
        #%% Define facet areas
        #===========================================
        #  !! OUTDATED FUNCTION, NEED UPDATED TO FENICS-X !!
        self.facets_area_define(domain, Q, self.flexdofs, gdim)

        #===========================================
        #%% Prepare post-process files
        #===========================================
        # !! OUTDATED FUNCTION, NEED UPDATED TO FENICS-X !!
        self.Create_Post_Process_Files(domain)

        #===========================================
        #%% Define the variational FORM
        #%% and
        #%% Jacobin functions of structure
        #===========================================

        if self.rank == 0: print ("{FENICS} Defining variational FORM functions ...   ", end="", flush=True)
        # Define the traction terms of the structure variational form
        tF = ufl.dot(chi, self.tF_apply)

        Form_s_Update_Acce = self.AMCK (ddmck, d0mck, u0mck, a0mck, beta_gam)
        Form_s_Update_velo = self.UMCK (Form_s_Update_Acce, u0mck, a0mck, gamma_gam)

        Form_s_Ga_Acce = self.Generalized_Alpha_Weights(Form_s_Update_Acce,a0mck,alpha_m_gam)
        Form_s_Ga_velo = self.Generalized_Alpha_Weights(Form_s_Update_velo,u0mck,alpha_f_gam)
        Form_s_Ga_disp = self.Generalized_Alpha_Weights(ddmck,d0mck,alpha_f_gam)
        Form_s_M_Matrix = self.rho_s() * ufl.inner(Form_s_Ga_Acce, chi) * ufl.dx
        Form_s_M_for_C_Matrix = self.rho_s() * ufl.inner(Form_s_Ga_velo, chi) * ufl.dx
        Form_s_K_Matrix = ufl.inner(self.elastic_stress(Form_s_Ga_disp,gdim), ufl.sym(ufl.grad(chi))) * ufl.dx
        Form_s_K_for_C_Matrix = ufl.inner(self.elastic_stress(Form_s_Ga_velo,gdim), ufl.sym(ufl.grad(chi))) * ufl.dx
        Form_s_C_Matrix = alpha_rdc * Form_s_M_for_C_Matrix + beta_rdc * Form_s_K_for_C_Matrix
        Form_s_F_Ext = tF * ds(2)

        Form_s = Form_s_M_Matrix + Form_s_C_Matrix + Form_s_K_Matrix - Form_s_F_Ext

        Bilinear_Form = fem.form(ufl.lhs(Form_s))
        Linear_Form   = fem.form(ufl.rhs(Form_s))

        Bilinear_Assemble = fem.petsc.assemble_matrix(Bilinear_Form, bcs=bcs)
        Bilinear_Assemble.assemble()
        Linear_Assemble = fem.petsc.create_vector(Linear_Form)

        if self.rank == 0: print ("Done")

        #===========================================
        #%% Initialize solver
        #===========================================

        solver = PETSc.KSP().create(domain.comm)
        solver.setOperators(Bilinear_Assemble)
        solver.setType(PETSc.KSP.Type.PREONLY)
        solver.getPC().setType(PETSc.PC.Type.LU)

        #===========================================
        #%% Setup checkpoint data
        #===========================================
        # !! OUTDATED FUNCTION, NEED UPDATED TO FENICS-X !!
        #self.Checkpoint_Output_Linear((t-self.dt()), mesh, d0mck, u0mck, a0mck, dmck, False)

        #===========================================
        #%% Define MUI samplers and commit ZERO step
        #===========================================
        # !! OUTDATED FUNCTION, NEED UPDATED TO FENICS-X !!
        self.MUI_Sampler_Define(Q, gdim, self.flexdofs, self.flexdofs, xyz_fetch, t_step)

        #===========================================
        #%% Define time loops
        #===========================================

        # Time-stepping
        while t <= self.T():

            # create an instance of the TicToc wall clock class
            wallClockPerStep = structureFSISolver.tictoc.TicToc()
            # Starts the wall clock
            wallClockPerStep.tic()

            # Update time list    
            times.append(t)
            n_steps = len(times)

            if self.rank == 0: 
                print ("\n")
                print ("{FENICS} Time: ", t, " [s]; Time Step Number: ", n_steps)

            # Change number of sub-iterations if needed
            if self.iChangeSubIter():
                if (t >= self.TChangeSubIter()):
                    present_num_sub_iteration = self.num_sub_iteration_new()
                else:
                    present_num_sub_iteration = self.num_sub_iteration()
            else:
                present_num_sub_iteration = self.num_sub_iteration()

            # Sub-iteration for coupling
            while i_sub_it <= present_num_sub_iteration:

                # Increment of total sub-iterations
                t_sub_it += 1

                if self.rank == 0: 
                    print ("\n")
                    print ("{FENICS} Sub-iteration Number: ", i_sub_it, " Total sub-iterations to now: ", t_sub_it)

                # Fetch and assign traction forces at present time step
                # !! OUTDATED FUNCTION, NEED UPDATED TO FENICS-X !!
                self.Traction_Assign(xyz_fetch, self.flexdofs, t_sub_it, n_steps, t)

                if (not ((self.iContinueRun()) and (n_steps == 1))):
                    # Update the right hand side reusing the initial vector
                    with Linear_Assemble.localForm() as loc_b:
                         loc_b.set(0)

                    # Assemble linear form
                    fem.petsc.assemble_vector(Linear_Assemble, Linear_Form)
                    fem.petsc.apply_lifting(Linear_Assemble, [Bilinear_Form], [bcs])
                    Linear_Assemble.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
                    fem.petsc.set_bc(Linear_Assemble, bcs)
                    # Solving the structure functions inside the time loop
                    solver.setMonitor(lambda _, its, rnorm: print(f"Iteration: {its}, rel. residual: {rnorm}"))
                    solver.solve(Linear_Assemble, dmck.x.petsc_vec)
                    solver.view()
                    dmck.x.scatter_forward()

                    force_X = ufl.dot(self.tF_apply, self.X_direction_vector())*ds(2)
                    force_Y = ufl.dot(self.tF_apply, self.Y_direction_vector())*ds(2)
                    force_Z = ufl.dot(self.tF_apply, self.Z_direction_vector())*ds(2)

                    f_X_a = fem.assemble_scalar(fem.form(force_X))
                    f_Y_a = fem.assemble_scalar(fem.form(force_Y))
                    f_Z_a = fem.assemble_scalar(fem.form(force_Z))

                    print ("{FENICS} Total Force_X on structure: ", f_X_a, " at self.rank ", self.rank)
                    print ("{FENICS} Total Force_Y on structure: ", f_Y_a, " at self.rank ", self.rank)
                    print ("{FENICS} Total Force_Z on structure: ", f_Z_a, " at self.rank ", self.rank)

                else:
                    pass

                # Compute and print the displacement of monitored point
                self.print_Disp(domain,dmck)

                # MUI Push internal points and commit current steps
                # !! OUTDATED FUNCTION, NEED UPDATED TO FENICS-X !!
                if (self.iMUICoupling()):
                    if (len(xyz_push)!=0):
                        self.MUI_Push(xyz_push, dofs_push_list, dmck, t_sub_it)
                    else:
                        self.MUI_Commit_only(t_sub_it)
                else:
                    pass

                # Increment of sub-iterations
                i_sub_it += 1

            # Mesh motion
            # !! OUTDATED FUNCTION, NEED UPDATED TO FENICS-X !!
            #self.Move_Mesh(V, dmck, d0mck, mesh)

            # Data output
            #  !! OUTDATED FUNCTION, NEED UPDATED TO FENICS-X !!
            if (not (self.iQuiet())):
                self.Export_Disp_xdmf(n_steps, t, domain, gdim, V, V1, dmck)
                self.Export_Disp_txt(domain,dmck)
                #self.Checkpoint_Output_Linear(t, mesh, d0mck, u0mck, a0mck, dmck, False)

            # Function spaces time marching
            #  !! OUTDATED FUNCTION, NEED UPDATED TO FENICS-X !!
            amck = self.AMCK (dmck.x.array, d0mck.x.array, u0mck.x.array, a0mck.x.array, beta_gam)
            umck = self.UMCK (amck, u0mck.x.array, a0mck.x.array, gamma_gam)

            a0mck.x.array[:] = amck
            u0mck.x.array[:] = umck
            d0mck.x.array[:] = dmck.x.array

            # Sub-iterator counter reset
            i_sub_it = 1
            # Physical time marching
            t += self.dt()

            # Finish the wall clock
            simtimePerStep = wallClockPerStep.toc()
            if self.rank == 0:
                print ("\n")
                print ("{FENICS} Simulation time per step: %g [s] at timestep: %i" % (simtimePerStep, n_steps))

        #===============================================
        #%% MPI barrier to wait for all solver to finish
        #===============================================

        # Wait for the other solver
        if self.iMUICoupling():
            #self.ifaces3d["threeDInterface0"].barrier(t_sub_it)
            iface3d.barrier(t_sub_it)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FILE END  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
