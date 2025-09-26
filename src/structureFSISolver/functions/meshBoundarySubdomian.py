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

    @file meshBoundarySubdomian.py
    
    @author W. Liu
    
    @brief mesh boundary and sub-domain related file of the structure code.

"""

#_________________________________________________________________________________________
#
#%% Import packages
#_________________________________________________________________________________________
from dolfinx import *
from dolfinx.mesh import (CellType, GhostMode, create_box,
                          locate_entities_boundary)
import numpy as np
import ufl

class meshBoundarySubdomian:

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #%% Solid Mesh input/generation
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def Mesh_Generation(self):
        # Restart simulation
        if self.iMeshLoad():
            # Load mesh from XDMF file
            if self.rank == 0: print ("{FENICS} Loading XDMF mesh ...   ")
            xdmfContinueRunMeshLoad = io.XDMFFile(self.LOCAL_COMM_WORLD, self.inputFolderPath + "/Structure_FEniCS.xdmf", "r")
            domain = xdmfContinueRunMeshLoad.read_mesh()
            if self.rank == 0: print ("{FENICS} Done with loading XDMF mesh")
        else:
            # Generate mesh
            if self.rank == 0: print ("{FENICS} Generating mesh ...   ")
            domain = mesh.create_box(self.LOCAL_COMM_WORLD,
                                   [[self.OBeamX(), self.OBeamY(), self.OBeamZ()],
                                   [(self.OBeamX()+self.XBeam()),
                                   (self.OBeamY()+self.YBeam()),
                                   (self.OBeamZ()+self.ZBeam())]],
                                   [self.XMesh(), self.YMesh(), self.ZMesh()],
                                   mesh.CellType.tetrahedron, ghost_mode=GhostMode.shared_facet)
            if self.rank == 0: print ("{FENICS} Done with generating mesh")

        if self.iXDMFFileExport() and self.iXDMFMeshExport() and (not self.iMeshLoad()):
            if self.rank == 0: print ("{FENICS} Exporting XDMF mesh ...   ", end="", flush=True)
            xdmfMeshExport = io.XDMFFile(self.LOCAL_COMM_WORLD, self.outputFolderPath + "/Structure_FEniCS.xdmf", "w")
            xdmfMeshExport.write_mesh(domain)
            if self.rank == 0: print ("Done")

        return domain

    def Get_Grid_Dimension(self, domain):
        # Geometry dimensions
        grid_dimension = domain.topology.dim
        return grid_dimension

    def Get_Face_Normal(self, domain):
        face_narmal = ufl.FacetNormal(domain)
        return face_narmal

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #%% Define SubDomains and boundaries
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def Boundaries_Generation_Fixed_Flex_Sym (self, domain, VectorFunctionSpace):

        #Define facets
        if self.iMeshLoad() and self.iSubdomainsImport():
            if self.rank == 0: print ("{FENICS} Loading HDF5 facets ...   ", end="", flush=True)
            if self.rank == 0: print ("Error, not implemented yet")
            if self.rank == 0: print ("Done")

        else:
            if self.rank == 0: print ("{FENICS} Creating facets ...   ", end="", flush=True)
            fdim = domain.topology.dim - 1

            fixed_facets    = mesh.locate_entities_boundary(domain, fdim, self.subDomains.Fixed)
            flex_facets     = mesh.locate_entities_boundary(domain, fdim, self.subDomains.Flex)
            symmetry_facets = mesh.locate_entities_boundary(domain, fdim, self.subDomains.Symmetry)

            if self.rank == 0: print ("Done")

        if self.iXDMFFileExport() and self.iHDF5SubdomainsExport():
            if self.rank == 0: print ("{FENICS} Exporting facets ...   ", end="", flush=True) 
            if self.rank == 0: print ("Error, not implemented yet")
            if self.rank == 0: print ("Done")

        #Define and mark mesh boundaries
        if self.iMeshLoad() and self.iBoundariesImport():
            if self.rank == 0: print ("{FENICS} Loading HDF5 boundaries ...   ", end="", flush=True)
            if self.rank == 0: print ("Error, not implemented yet")
            if self.rank == 0: print ("Done")

        else:
            if self.rank == 0: print ("{FENICS} Creating boundaries ...   ", end="", flush=True)

            self.fixedmt = mesh.meshtags(domain, fdim, fixed_facets, 1)
            self.flexmt = mesh.meshtags(domain, fdim, flex_facets, 2)
            self.symmetrymt = mesh.meshtags(domain, fdim, symmetry_facets, 3)

            self.fixeddofs = fem.locate_dofs_topological(VectorFunctionSpace, fdim, fixed_facets)
            self.flexdofs = fem.locate_dofs_topological(VectorFunctionSpace, fdim, flex_facets)
            self.symmetrydofs = fem.locate_dofs_topological(VectorFunctionSpace, fdim, symmetry_facets)

            marked_facets = np.hstack([fixed_facets, flex_facets,symmetry_facets])
            marked_values = np.hstack([np.full_like(fixed_facets, 1), np.full_like(flex_facets, 2), np.full_like(symmetry_facets, 3)])
            sorted_facets = np.argsort(marked_facets)
            self.facet_tag = mesh.meshtags(domain, fdim, marked_facets[sorted_facets], marked_values[sorted_facets])

            if self.rank == 0: print ("Done")

        if self.iXDMFFileExport() and self.iHDF5BoundariesExport():
            if self.rank == 0: print ("{FENICS} Exporting HDF5 boundaries ...   ", end="", flush=True)
            if self.rank == 0: print ("Error, not implemented yet")
            if self.rank == 0: print ("Done")

        if self.rank == 0: 
            print ("\n")
            print ("{FENICS} Structure Mesh Info: ")
            #print ("{FENICS} (geometry dimension, Dofs): ",VectorFunctionSpace.shape)
            print ("{FENICS} Cells:", domain.topology.index_map(domain.topology.dim).size_local)
            print ("\n")

    def Get_ds(self, domain):
        metadata = {"quadrature_degree": 4}
        return ufl.Measure('ds', domain=domain, subdomain_data=self.facet_tag, metadata=metadata)

    def Get_dx(self, domain):
        metadata = {"quadrature_degree": 4}
        return ufl.Measure("dx", domain=domain, metadata=metadata)
