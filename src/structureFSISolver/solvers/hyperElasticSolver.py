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

    @file hyperElastic.py

    @author W. Liu

    @brief The core hyper elastic solver class of the structure code.

"""

#_________________________________________________________________________________________
#
#%% Import packages
#_________________________________________________________________________________________

from dolfinx import *
from dolfinx.nls.petsc import NewtonSolver
import os
import numpy as np
from mpi4py import MPI
import structureFSISolver
import ufl

class hyperElastic:

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #%% Main solver function
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def hyperElasticSolve(self):

        #===========================================
        #%% Time marching parameters define
        #===========================================

        t        = self.Start_Time
        t_step   = self.Time_Steps
        i_sub_it = self.Start_Number_Sub_Iteration

        #===========================================
        #%% Solid Mesh input/generation
        #===========================================

        mesh = self.Mesh_Generation()
        gdim = self.Get_Grid_Dimension(mesh)
        N    = self.Get_Face_Normal(mesh)

        #===========================================
        #%% Define coefficients
        #===========================================

        # Time step constants
        k = fem.Constant(mesh, self.dt())

        # Time lists
        times    = []
        t_sub_it = 0

        # One-step theta value
        theta = fem.Constant(mesh, self.thetaOS())

        if self.rank == 0:
            print ("\n")
            print ("{FENICS} One-step theta: ", float(theta))
            print ("\n")

        #===========================================
        #%% Define function spaces
        #===========================================

        if self.rank == 0: print ("{FENICS} Creating function spaces ...   ")

#        V_ele     =     ufl.VectorElement("Lagrange", mesh.ufl_cell(), self.deg_fun_spc()) # Displacement & Velocity Vector element

        Q         =     fem.functionspace(mesh, ("Lagrange", self.deg_fun_spc()))            # Function space with updated mesh
#        VV        =     FunctionSpace(mesh, MixedElement([V_ele, V_ele]))            # Mixed (Velocity (w) & displacement (d)) function space
        V         =     fem.functionspace(mesh, ("Lagrange", self.deg_fun_spc(), (mesh.geometry.dim, )))
        V1         =     fem.functionspace(mesh, ("Lagrange", 1, (mesh.geometry.dim, )))
#        T_s_space =     TensorFunctionSpace(mesh, 'Lagrange', self.deg_fun_spc())      # Define nth order structure function spaces

        if self.rank == 0: print ("{FENICS} Done with creating function spaces")

        #=======================================================
        #%% Define functions, test functions and trail functions
        #=======================================================

        if self.rank == 0: print ("{FENICS} Creating functions, test functions and trail functions ...   ", end="", flush=True)

        # Test functions
#        psi, phi = TestFunctions(VV)    # Test functions for velocity and displacement

        # Functions at present time step
#        ud   = Function(VV)               # Functions for velocity and displacement
#        u, d = split(ud)                # Split velocity and displacement functions

        # Functions at previous time step
#        u0d0   = Function(VV)             # Functions for velocity and displacement
#        u0, d0 = split(u0d0)            # Split velocity and displacement functions

        # Define structure traction
#        sigma_s = Function(T_s_space)   # Structure traction normal to structure

#        self.Load_Functions_Continue_Run_Nonlinear(u0d0,ud,sigma_s)
        u = fem.Function(V, name="Displacement")
        u0 = fem.Function(V, name="Displacement0")

        v = ufl.TestFunction(V)
        du = ufl.TrialFunction(V)

        if self.rank == 0: print ("Done")

        #===========================================
        #%% Define traction forces
        #===========================================

        self.Traction_Define(V)

        #===========================================
        #%% Define SubDomains and boundaries
        #===========================================

        self.Boundaries_Generation_Fixed_Flex_Sym(mesh, V)

        ds = self.Get_ds(mesh)

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

        dofs_fetch_list = self.dofs_list(mesh, V, 2)

        xyz_fetch = self.xyz_np(dofs_fetch_list, V, gdim)

        dofs_push_list = self.dofs_list(mesh, V, 2)

        xyz_push = self.xyz_np(dofs_push_list, V, gdim)

        #===========================================
        #%% Define facet areas
        #===========================================

        self.facets_area_define(mesh, Q, self.flexdofs, gdim)

        #===========================================
        #%% Prepare post-process files
        #===========================================

        self.Create_Post_Process_Files(mesh)

        #===========================================
        #%% Define the variational FORM
        #%% and
        #%% Jacobin functions of structure
        #===========================================

        if self.rank == 0: print ("{FENICS} Defining variational FORM and Jacobin functions ...   ", end="", flush=True)

        # Define the traction terms of the structure variational form
        tF = ufl.dot(self.F_(u,gdim).T, self.tF_apply)
        tF_ = ufl.dot(self.F_(u0,gdim).T, self.tF_apply)

        # Define the transient terms of the structure variational form
        Form_s_T = (1/k)*ufl.inner((u-u0), v)*ufl.dx

        # Define the stress terms and convection of the structure variational form
        if self.iNonLinearMethod():
            if self.rank == 0: print ("{FENICS} [Defining non-linear stress-strain relation: Define the First Piola-Kirchhoff stress tensor by the constitutive law of hyper-elastic St. Vernant-Kirchhoff material model (non-linear relation). Valid for large deformations but small strain] ...   ", end="", flush=True)
            Form_s_SC = ufl.inner(theta * self.Piola_Kirchhoff_fst(u,gdim) + (1 - theta) *
                        self.Piola_Kirchhoff_fst(u0,gdim), ufl.grad(v)) * ufl.dx
        else:
            if self.rank == 0: print ("{FENICS} [Defining linear stress-strain relation: Define the First Piola-Kirchhoff stress tensor by Hooke's law (linear relation). Valid for small-scale deformations only] ...   ", end="", flush=True)
            Form_s_SC = ufl.inner(theta * self.Hooke_stress(u,gdim) + (1 - theta) *
                        self.Hooke_stress(u0,gdim), ufl.grad(v)) * ufl.dx

        # Define the body forces and surface tractions terms of the structure variational form
        Form_s_ET = -( theta * self.J_(u,gdim) * ufl.inner( (self.b_for(mesh)), v ) +
                    ( 1 - theta ) * self.J_(u0,gdim) * ufl.inner( (self.b_for(mesh)), v ) ) * ufl.dx
        Form_s_ET -= ( theta * self.J_(u,gdim) * ufl.inner( tF, v ) +
                    ( 1 - theta ) * self.J_(u0,gdim) * ufl.inner( tF_, v ) ) * ds(2)
        Form_s_ET -= ( theta * self.J_(u,gdim) * ufl.inner( ufl.inv(self.F_(u,gdim)) * u, v )+
                    ( 1 - theta ) * (self.J_(u0,gdim)) * ufl.inner(ufl.inv(self.F_(u0,gdim)) * u, v )) * ds(2)

        # Define the final form of the structure variational form
        Form_s = Form_s_T + Form_s_SC + Form_s_ET

        # Define Jacobin functions
        Jaco = ufl.derivative(Form_s, u, du)

        if self.rank == 0: print ("Done")

        #===========================================
        #%% Initialize solver
        #===========================================

        problem = fem.petsc.NonlinearProblem(Form_s, u, bcs)
        solver = NewtonSolver(mesh.comm, problem)
        # Set Newton solver options
        solver.atol = 1e-4
        solver.rtol = 1e-4
        solver.convergence_criterion = "incremental"

        #===========================================
        #%% Setup checkpoint data
        #===========================================

        # self.Checkpoint_Output_Nonlinear((t-self.dt()), mesh, u0, u, sigma_s, False)

        #===========================================
        #%% Define MUI samplers and commit ZERO step
        #===========================================

        self.MUI_Sampler_Define(Q, gdim, dofs_fetch_list, dofs_push_list, xyz_fetch, t_step)

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
                self.Traction_Assign(xyz_fetch, dofs_fetch_list, t_sub_it, n_steps, t)

                if (not ((self.iContinueRun()) and (n_steps == 1))):
                    # Solving the structure functions inside the time loop
                    num_its, converged = solver.solve(u)
                    assert converged

                    u.x.scatter_forward()  # updates ghost values for parallel computations

                    print(
                        f"Time step {n_steps}, Number of iterations {num_its}."
                    )

                    force_X = ufl.dot(tF, self.X_direction_vector())*ds(2)
                    force_Y = ufl.dot(tF, self.Y_direction_vector())*ds(2)
                    force_Z = ufl.dot(tF, self.Z_direction_vector())*ds(2)

                    f_X_a = fem.assemble_scalar(fem.form(force_X))
                    f_Y_a = fem.assemble_scalar(fem.form(force_Y))
                    f_Z_a = fem.assemble_scalar(fem.form(force_Z))

                    print ("{FENICS} Total Force_X on structure: ", f_X_a, " at self.rank ", self.rank)
                    print ("{FENICS} Total Force_Y on structure: ", f_Y_a, " at self.rank ", self.rank)
                    print ("{FENICS} Total Force_Z on structure: ", f_Z_a, " at self.rank ", self.rank)

                else:
                    pass

                # Compute and print the displacement of monitored point
                self.print_Disp(mesh, u)

                # MUI Push internal points and commit current steps
                if (self.iMUICoupling()):
                    if (len(xyz_push)!=0):
                        self.MUI_Push(xyz_push, dofs_push_list, u, t_sub_it)
                    else:
                        self.MUI_Commit_only(t_sub_it)
                else:
                    pass

                # Increment of sub-iterations
                i_sub_it += 1

            # Mesh motion
            #self.Move_Mesh(V, d, d0, mesh)

            # Data output
            if (not (self.iQuiet())):
                self.Export_Disp_xdmf(n_steps, t, mesh, gdim, V, V1, u)
                self.Export_Disp_txt(mesh,u)
                # self.Checkpoint_Output_Nonlinear(t, mesh, u0, u, False)

            # Assign the old function spaces
            #u0.assign(u)
            u0.x.array[:] = u.x.array

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
            self.ifaces3d["threeDInterface0"].barrier(t_sub_it)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FILE END  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#