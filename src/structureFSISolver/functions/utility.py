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

    @file utility.py
    
    @author W. Liu
    
    @brief utility file of the structure code.

"""

#_________________________________________________________________________________________
#
#%% Import packages
#_________________________________________________________________________________________
from dolfinx import fem
from dolfinx.fem import Constant
import numpy as np
import ufl

class utility:

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #%% Define directional vectors
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def X_direction_vector(self):
        # Directional vector in x-axis direction
       return ufl.as_vector([1.0, 0.0, 0.0])
    def Y_direction_vector(self):
        # Directional vector in y-axis direction
       return ufl.as_vector([0.0, 1.0, 0.0])
    def Z_direction_vector(self):
        # Directional vector in z-axis direction
       return ufl.as_vector([0.0, 0.0, 1.0])

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #%% Solid gravitational/body forces define
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def b_for (self, mesh):
        # Body external forces define [N/m^3]
        b_for_ext = fem.Constant(mesh, (self.bForExtX(), self.bForExtY(), self.bForExtZ()))
        # Gravitational force define [N/m^3]
        if self.iGravForce():
            g_force = fem.Constant(mesh, (0.0, (self.rho_s() * (-9.81)), 0.0))
        else:
            g_force = fem.Constant(mesh, (0.0, (0.0 * (-9.81)), 0.0))
        return (b_for_ext + g_force)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #%% Define Stress, force gradient and its
    #%% determination functions
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def I(self, grid_dimension):
        # Define the Identity matrix
        return (ufl.Identity(grid_dimension))

    def F_(self, displacement_function, grid_dimension):
        # Define the deformation gradient
        return (self.I(grid_dimension) + ufl.nabla_grad(displacement_function))

    def J_(self, displacement_function, grid_dimension):
        # Define the determinant of the deformation gradient
        return ufl.det(self.F_(displacement_function,grid_dimension))
        # return np.linalg.det(self.F_(displacement_function,grid_dimension))

    def C(self, displacement_function, grid_dimension):
        # Define the right Cauchy-Green strain tensor
        return (ufl.transpose(self.F_(displacement_function, grid_dimension)) *
                self.F_(displacement_function, grid_dimension))

    def E(self, displacement_function, grid_dimension):
        # Define the non-linear Lagrangian Green strain tensor
        return (0.5 * (self.C(displacement_function, grid_dimension) - self.I(grid_dimension)))

    def epsilon(self, displacement_function, grid_dimension):
        # Define the linear Lagrangian Green strain tensor
        return (0.5 * (ufl.nabla_grad(displacement_function) + ufl.transpose(ufl.nabla_grad(displacement_function))))

    def sigma(self, displacement_function, grid_dimension):
        # Define the linear Lagrangian Green strain tensor
        return (self.lamda_s * ufl.tr(ufl.sym(ufl.grad(displacement_function))) * ufl.Identity(len(displacement_function)) + 2*self.mu_s*self.epsilon(displacement_function, grid_dimension))

    def Piola_Kirchhoff_sec(self, displacement_function, strain_tensor, grid_dimension):
        # Define the Second Piola-Kirchhoff stress tensor by the constitutive law
        #   of hyper-elastic St. Vernant-Kirchhoff material model (non-linear relation).
        #   Valid for large deformations but small strain.
        return (self.lamda_s() * ufl.tr(strain_tensor(displacement_function, grid_dimension)) *
                self.I(grid_dimension) + 2.0 * self.mu_s() *
                strain_tensor(displacement_function, grid_dimension))

    def cauchy_stress(self, displacement_function, strain_tensor, grid_dimension):
        # Define the Cauchy stress tensor
        return ((1 / self.J_(displacement_function, grid_dimension)) *
                (self.F_(displacement_function, grid_dimension)) *
                (self.Piola_Kirchhoff_sec(displacement_function, strain_tensor,grid_dimension)) *
                ufl.transpose(self.F_(displacement_function, grid_dimension)))

    def Piola_Kirchhoff_fst(self, displacement_function, grid_dimension):
        # Define the First Piola-Kirchhoff stress tensor by the constitutive law
        #   of hyper-elastic St. Vernant-Kirchhoff material model (non-linear relation).
        #   Valid for large deformations but small strain.
        return (self.J_(displacement_function, grid_dimension) *
                self.cauchy_stress(displacement_function, self.E,grid_dimension) *
                ufl.transpose(ufl.inv(self.F_(displacement_function, grid_dimension))))

    def Hooke_stress(self, displacement_function, grid_dimension):
        # Define the First Piola-Kirchhoff stress tensor by Hooke's law (linear relation).
        #   Valid for small-scale deformations only.
        return (self.J_(displacement_function, grid_dimension) *
                self.cauchy_stress(displacement_function, self.epsilon, grid_dimension) *
                ufl.transpose(ufl.inv(self.F_(displacement_function, grid_dimension))))

    def elastic_stress(self, displacement_function, grid_dimension):
        # Define the elastic stress tensor
        return (2.0 * self.mu_s() * ufl.sym(ufl.grad(displacement_function)) +
                self.lamda_s() * ufl.tr(ufl.sym(ufl.grad(displacement_function))) * self.I(grid_dimension))

    def Traction_Define(self, VectorFunctionSpace):
        # !! OUTDATED FUNCTION, NEED UPDATED TO FENICS-X !!
        if self.iNonUniTraction():
            if self.rank == 0: print ("{FENICS} Non-uniform traction applied")
            self.tF_apply = fem.Function(VectorFunctionSpace)
            self.tF_apply_vec = self.tF_apply.x.array
        else:
            if self.rank == 0: print ("{FENICS} Uniform traction applied")
            self.tF_magnitude = Constant(0.0 *self.X_direction_vector() +
                                    0.0 *self.Y_direction_vector() +
                                    0.0 *self.Z_direction_vector())
            self.tF_apply = self.tF_magnitude

    def Traction_Assign(self, xyz_fetch, dofs_fetch_list, t_sub_it, n_steps, t):
        # Assign traction forces at present time step
        if self.iNonUniTraction():
            if len(xyz_fetch)!=0:
                # Execute only when there are DoFs need to exchange data in this rank.
                if self.iMUICoupling():
                    self.MUI_Fetch(xyz_fetch, dofs_fetch_list, t_sub_it)
                else:
                    self.Traction_DoF_Assign(xyz_fetch, dofs_fetch_list, t_sub_it, t)
            if (self.iMUIFetchValue()) and (not ((self.iContinueRun()) and (n_steps == 1))):
                # Apply traction components. These calls do parallel communication
                # self.tF_apply.vector().set_local(self.tF_apply_vec)
                # self.tF_apply.vector().apply("insert")
                self.tF_apply.x.array[:] = self.tF_apply_vec
                self.tF_apply.x.scatter_forward()

            else:
                # Do not apply the fetched value, i.e. one-way coupling
                pass
        else:
            if self.rank == 0: print ("{FENICS} Assigning uniform traction forces at present time step ...   ",
                                    end="", flush=True)
            if (t <= self.sForExtEndTime()):
                self.tF_magnitude.assign((Constant((self.sForExtX()) /
                                                   (self.YBeam() * self.ZBeam())) *
                                                   self.X_direction_vector()) +
                                                   (Constant((self.sForExtY()) /
                                                   (self.XBeam()*self.ZBeam())) *
                                                   self.Y_direction_vector()) +
                                                   (Constant((self.sForExtZ()) /
                                                   (self.XBeam()*self.YBeam())) *
                                                   self.Z_direction_vector()))
            else:
                self.tF_magnitude.assign(Constant((0.0)))
            if self.rank == 0:
                print ("Done")

    def Traction_DoF_Assign(self, dofs_to_xyz, dofs_fetch_list, total_Sub_Iteration, t):
        totForceX = 0.0
        totForceY = 0.0
        totForceZ = 0.0
        temp_vec_function_temp = self.tF_apply_vec

        for i, p in enumerate(dofs_fetch_list):
            if (t <= self.sForExtEndTime()):
                if self.iConstantSForExt():
                    self.tF_apply_vec[0::3][p] = (self.sForExtX()) / len(dofs_fetch_list)
                    self.tF_apply_vec[1::3][p] = (self.sForExtY()) / len(dofs_fetch_list)
                    self.tF_apply_vec[2::3][p] = (self.sForExtZ()) / len(dofs_fetch_list)
                else:
                    self.tF_apply_vec[0::3][p] = (t / self.sForExtEndTime()) * (self.sForExtX()) / len(dofs_fetch_list)
                    self.tF_apply_vec[1::3][p] = (t / self.sForExtEndTime()) * (self.sForExtY()) / len(dofs_fetch_list)
                    self.tF_apply_vec[2::3][p] = (t / self.sForExtEndTime()) * (self.sForExtZ()) / len(dofs_fetch_list)
            else:
                if self.iConstantSForExt():
                    self.tF_apply_vec[0::3][p] = (0.0) / len(dofs_fetch_list)
                    self.tF_apply_vec[1::3][p] = (0.0) / len(dofs_fetch_list)
                    self.tF_apply_vec[2::3][p] = (0.0) / len(dofs_fetch_list)
                else:
                    self.tF_apply_vec[0::3][p] = (0.0) * (self.sForExtX()) / len(dofs_fetch_list)
                    self.tF_apply_vec[1::3][p] = (0.0) * (self.sForExtY()) / len(dofs_fetch_list)
                    self.tF_apply_vec[2::3][p] = (0.0) * (self.sForExtZ()) / len(dofs_fetch_list)

            totForceX += self.tF_apply_vec[0::3][p]
            totForceY += self.tF_apply_vec[1::3][p]
            totForceZ += self.tF_apply_vec[2::3][p]

            self.tF_apply_vec[0::3][p] /= self.areaf_vec[p]
            self.tF_apply_vec[1::3][p] /= self.areaf_vec[p]
            self.tF_apply_vec[2::3][p] /= self.areaf_vec[p]

        if self.iDebug():
            print ("{FENICS**} totForce Apply: ", totForceX, "; ",totForceY, "; ",totForceZ,
                    "; at iteration: ", total_Sub_Iteration, " at rank: ", self.rank)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FILE END  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
