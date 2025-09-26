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

    @file cfgPrsFn.py
    
    @author W. Liu
    
    @brief configure parameters function file of the structure code.

"""

#_________________________________________________________________________________________
#
#%% Import packages
#_________________________________________________________________________________________

class readData:

        #===========================================
        #%% Debug mode on/off switch
        #===========================================
    def iDebug (self):
        # F-Switch off the debug level codes (if any); T-Switch on the debug level codes (if any).
        return self.cfg['LOGGING'].getboolean('iDebug')

    def iQuiet (self):
        # F-Switch off the Quiet; T-Switch on the Quiet.
        return self.cfg['LOGGING'].getboolean('iQuiet')

        #===========================================
        #%% MUI switches & parameters
        #===========================================
    def iMUICoupling (self):
        # F-Switch off the MUI coupling function; T-Switch on the MUI coupling function.
        return self.cfg['MUI'].getboolean('iMUICoupling')
    def iMUIFetchValue (self):
        # F-Switch off the MUI Fetch; T-Switch on the MUI Fetch.
        return self.cfg['MUI'].getboolean('iMUIFetchValue')
    def iMUIFetchMany (self):
        # F-Use normal fetch function; T-Use fetch_many function.
        return self.cfg['MUI'].getboolean('iMUIFetchMany')
    def iLoadAreaList (self):
        # F-The mesh is not evenly spaced; T-The mesh is evenly spaced.
        return self.cfg['MUI'].getboolean('iLoadAreaList')
    def areaListFactor (self):
        # Factor for area list calculation (float)
        return float(self.cfg['MUI']['areaListFactor'])
    def iMUIPushMany (self):
        # F-Use normal push function; T-Use push_many function.
        return self.cfg['MUI'].getboolean('iMUIPushMany')
    def iPushX (self):
        # F-Not push X; T-Push X.
        return self.cfg['MUI'].getboolean('iPushX')
    def iPushY (self):
        # F-Not push Y; T-Push Y.
        return self.cfg['MUI'].getboolean('iPushY')
    def iPushZ (self):
        # F-Not push Z; T-Push Z.
        return self.cfg['MUI'].getboolean('iPushZ')
    def rMUIFetcher (self):
        # Spatial sampler search radius (float)
        return float(self.cfg['MUI']['rMUIFetcher'])
    def iConservative (self):
        # F-Use RBF spatial sampler consistent mode; T-Use RBF spatial sampler conservative mode.
        return self.cfg['MUI'].getboolean('iConservative')
    def cutoffRBF (self):
        # RBF spatial sampler cutoff value (float)
        return float(self.cfg['MUI']['cutoffRBF'])
    def iWriteMatrix (self):
        # F-The RBF matrix will read from file; T-The RBF matrix will write to file.
        return self.cfg['MUI'].getboolean('iWriteMatrix')
    def iPolynomial (self):
        # F-Switch off the RBF spatial sampler polynomial terms; T-Switch on the RBF spatial sampler polynomial terms.
        return self.cfg['MUI'].getboolean('iPolynomial')
    def basisFunc (self):
        # Select of basis functions of the RBF spatial sampler (integer)
        return int(self.cfg['MUI']['basisFunc'])
    def iSmoothFunc (self):
        # F-Switch off the RBF spatial sampler smooth function; T-Switch on the RBF spatial sampler smooth function.
        return self.cfg['MUI'].getboolean('iSmoothFunc')
    def cgSolveTolRBF (self):
        # RBF spatial sampler cg solver tol value (float)
        return float(self.cfg['MUI']['cgSolveTolRBF'])
    def cgMaxIterRBF (self):
        # RBF spatial sampler cg solver maximum iterator (integer)
        return int(self.cfg['MUI']['cgMaxIterRBF'])
    def pouSizeRBF (self):
        # RBF spatial sampler POU size (integer)
        return int(self.cfg['MUI']['pouSizeRBF'])
    def precondRBF (self):
        # Select of pre-conditionore of the RBF spatial sampler (integer)
        return int(self.cfg['MUI']['precondRBF'])
    def forgetTStepsMUI (self):
        # Numbers of time steps to forget for MUI push (integer)
        return int(self.cfg['MUI']['forgetTStepsMUI'])
    def iparallelFSICoupling (self):
        # F-Serial FSI coupling mode; T-Parallel FSI coupling mode.
        return self.cfg['MUI'].getboolean('iparallelFSICoupling')
    def undRelxCpl (self):
        # Initial under relaxation factor for IQNILS (float)
        return float(self.cfg['MUI']['undRelxCpl'])

        #===========================================
        #%% Global solver define
        #===========================================

    def solving_method (self):
        # define the solving Method; "STVK" "MCK"
        return self.cfg['SOLVER']['solving_method']
    def linear_solver (self):
        # define the linear solver; "LU" "LinearVariational"
        return self.cfg['SOLVER']['linear_solver']
    def nonlinear_solver (self):
        # define the non-linear solver; "snes" "newton"
        return self.cfg['SOLVER']['nonlinear_solver']
    def prbsolver (self):
        # define the linear solver for the problem
        return self.cfg['SOLVER']['prbsolver']
    def prjsolver (self):
        # define the solver for project between domains
        return self.cfg['SOLVER']['prjsolver']
    def prbpreconditioner (self):
        # define the pre-conditioner for the problem
        return self.cfg['SOLVER']['prbpreconditioner']
    def lineSearch (self):
        # define the line search for the snes solver
        return self.cfg['SOLVER']['lineSearch']
    def prbRelative_tolerance (self):
        # define the relative tolerance for the problem
        return float(self.cfg['SOLVER']['prbRelative_tolerance'])
    def prbAbsolute_tolerance (self):
        # define the absolute tolerance for the problem
        return float(self.cfg['SOLVER']['prbAbsolute_tolerance'])
    def prbMaximum_iterations (self):
        # define the maximum iterations for the problem
        return int(self.cfg['SOLVER']['prbMaximum_iterations'])
    def prbRelaxation_parameter (self):
        # define the relaxation parameter for the problem
        return float(self.cfg['SOLVER']['prbRelaxation_parameter'])
    def compRepresentation (self):
        # define the representation of the compiler
        return self.cfg['SOLVER']['compRepresentation']
    def cppOptimize (self):
        # switch on the C++ code optimization
        return self.cfg['SOLVER'].getboolean('cppOptimize')
    def optimize (self):
        # switch on optimization of the compiler
        return self.cfg['SOLVER'].getboolean('optimize')
    def allow_extrapolation (self):
        # switch on extrapolation WARRING: Please set it 'FALSE' for Parallel Interpolation
        return self.cfg['SOLVER'].getboolean('allow_extrapolation')
    def ghost_mode (self):
        # Ghost cell mode: "shared_facet"; "shared_vertex"; "none"
        return self.cfg['SOLVER']['ghost_mode']
    def error_on_nonconvergence (self):
        # switch on error of non convergence
        return self.cfg['SOLVER'].getboolean('error_on_nonconvergence')
    def krylov_maximum_iterations (self):
        # define the maximum iterations for the krylov solver
        return int(self.cfg['SOLVER']['krylov_maximum_iterations'])
    def krylov_prbRelative_tolerance (self):
        # define the relative tolerance for the krylov solver
        return float(self.cfg['SOLVER']['krylov_prbRelative_tolerance'])
    def krylov_prbAbsolute_tolerance (self):
        # define the absolute tolerance for the krylov solver
        return float(self.cfg['SOLVER']['krylov_prbAbsolute_tolerance'])
    def monitor_convergence (self):
        # switch on monitor convergence for the krylov solver
        return self.cfg['SOLVER'].getboolean('monitor_convergence')
    def nonzero_initial_guess (self):
        # switch on nonzero initial guess for the krylov solver
        return self.cfg['SOLVER'].getboolean('nonzero_initial_guess')
    def show_report (self):
        # switch on report for the krylov solver
        return self.cfg['SOLVER'].getboolean('show_report')

        #===========================================
        #%% Global degree orders
        #===========================================

    def deg_fun_spc (self):
        # Function space degree order
        return int(self.cfg['ORDER']['deg_fun_spc'])

        #===========================================
        #%% Target folder input
        #===========================================

        # F-Input/output folder directories are relative paths; T-Input/output folder directories are absolute paths.
    def iAbspath (self):
        return self.cfg['FOLDER'].getboolean('iAbspath')
    def outputFolderName (self):
        return self.cfg['FOLDER']['outputFolderName']
    def inputFolderName (self):
        return self.cfg['FOLDER']['inputFolderName']

        #===========================================
        #%% Solid mechanical parameters input
        #===========================================

    def E_s (self):
        # Young's Modulus [Pa] (5.0e5) (1.4e6) (1.0e4)
        return float(self.cfg['MECHANICAL']['E_s'])
    def rho_s (self):
        # Density of solid [kg/m^3]
        return float(self.cfg['MECHANICAL']['rho_s'])
    def nu_s (self):
        # Poisson ratio [-]
        return float(self.cfg['MECHANICAL']['nu_s'])

        #===========================================
        #%% Solid body external forces input
        #===========================================

    def bForExtX (self):
        # Body external forces in x-axis direction [N/m^3]
        return float(self.cfg['EXTFORCE']['bForExtX'])
    def bForExtY (self):
        # Body external forces in y-axis direction [N/m^3]
        return float(self.cfg['EXTFORCE']['bForExtY'])
    def bForExtZ (self):
        # Body external forces in z-axis direction [N/m^3]
        return float(self.cfg['EXTFORCE']['bForExtZ'])
    def sForExtX (self):
        # Surface external forces in x-axis direction [N/m^2]
        return float(self.cfg['EXTFORCE']['sForExtX'])
    def sForExtY (self):
        # Surface external forces in y-axis direction [N/m^2]
        return float(self.cfg['EXTFORCE']['sForExtY'])
    def sForExtZ (self):
        # Surface external forces in z-axis direction [N/m^2]
        return float(self.cfg['EXTFORCE']['sForExtZ'])
    def sForExtEndTime (self):
        # Surface external forces end time [s]
        return float(self.cfg['EXTFORCE']['sForExtEndTime'])
    def iConstantSForExt (self):
        # F-Linear increment of external surface forces; T-Constant external surface forces.
        return self.cfg['EXTFORCE'].getboolean('iConstantSForExt')

        #===========================================
        #%% Time marching parameter input
        #===========================================

    def T (self):
        # End time [s]
        return float(self.cfg['TIME']['T'])
    def dt (self):
        # Time step size [s]
        return float(self.cfg['TIME']['dt'])
    def num_sub_iteration (self):
        # Numbers of sub-iterations (integer) [-]
        return int(self.cfg['TIME']['num_sub_iteration'])
    def iContinueRun (self):
        # F-Run from initial time step; T-Continue run based on previous results.
        return self.cfg['TIME'].getboolean('iContinueRun')
    def iResetStartTime (self):
        # F-Run from initial time; T-Run from a different time.
        return self.cfg['TIME'].getboolean('iResetStartTime')
    def newStartTime (self):
        # New start time (when iResetStartTime = True) [s]
        return float(self.cfg['TIME']['newStartTime'])
    def iChangeSubIter (self):
        # F-sub-iteration remains the same; T-change the sub-iteration number.
        return self.cfg['TIME'].getboolean('iChangeSubIter')
    def TChangeSubIter (self):
        # Time to change the sub-iteration [s]
        return float(self.cfg['TIME']['TChangeSubIter'])
    def num_sub_iteration_new (self):
        # New numbers of sub-iterations (integer) [-]
        return int(self.cfg['TIME']['num_sub_iteration_new'])

        #===========================================
        #%% Time marching accurate control
        #===========================================

        # One-step theta value, valid only on STVK solver
    def thetaOS (self):
        return float(self.cfg['TIMEMARCHCOEF']['thetaOS'])

        # Rayleigh damping coefficients, valid only on MCK solver
    def alpha_rdc (self):
        return float(self.cfg['TIMEMARCHCOEF']['alpha_rdc'])
    def beta_rdc (self):
        return float(self.cfg['TIMEMARCHCOEF']['beta_rdc'])

        # Generalized-alpha method parameters, valid only on MCK solver
        # alpha_m_gam <= alpha_f_gam <= 0.5 for a better performance
        # Suggested values for alpha_m_gam: 0.0 or 0.4
        # Suggested values for alpha_f_gam: 0.0 or 0.2
    def alpha_m_gam (self):
        return float(self.cfg['TIMEMARCHCOEF']['alpha_m_gam'])
    def alpha_f_gam (self):
        return float(self.cfg['TIMEMARCHCOEF']['alpha_f_gam'])

        #===========================================
        #%% Post-processing parameter input
        #===========================================

    def output_interval (self):
        # Output file intervals (integer) [-]
        return int(self.cfg['POSTPROCESS']['output_interval'])
    def pointMoniX (self):
        # X-axis coordinate of the monitoring point [m]
        return float(self.cfg['POSTPROCESS']['pointMoniX'])
    def pointMoniY (self):
        # Y-axis coordinate of the monitoring point [m]
        return float(self.cfg['POSTPROCESS']['pointMoniY'])
    def pointMoniZ (self):
        # Z-axis coordinate of the monitoring point [m]
        return float(self.cfg['POSTPROCESS']['pointMoniZ'])
    def pointMoniXb (self):
        # X-axis coordinate of the monitoring point [m]
        return float(self.cfg['POSTPROCESS']['pointMoniXb'])
    def pointMoniYb (self):
        # Y-axis coordinate of the monitoring point [m]
        return float(self.cfg['POSTPROCESS']['pointMoniYb'])
    def pointMoniZb (self):
        # Z-axis coordinate of the monitoring point [m]
        return float(self.cfg['POSTPROCESS']['pointMoniZb'])

        #===========================================
        #%% Solid Model dimension input
        #===========================================

    def OBeamX (self):
        # x coordinate of the original point of the beam [m]
        return float(self.cfg['GEOMETRY']['OBeamX'])
    def OBeamY (self):
        # y coordinate of the original point of the beam [m]
        return float(self.cfg['GEOMETRY']['OBeamY'])
    def OBeamZ (self):
        # z coordinate of the original point of the beam [m]
        return float(self.cfg['GEOMETRY']['OBeamZ'])
    def XBeam (self):
        # length of the beam [m]
        return float(self.cfg['GEOMETRY']['XBeam'])
    def YBeam (self):
        # width of the beam [m]
        return float(self.cfg['GEOMETRY']['YBeam'])
    def ZBeam (self):
        # thick of the beam [m]
        return float(self.cfg['GEOMETRY']['ZBeam'])

        #===========================================
        #%% Solid calculation selection
        #===========================================

    def iMeshLoad (self):
        # F-Generate mesh; T-Load mesh from file.
        return self.cfg['CALMODE'].getboolean('iMeshLoad')
    def iNonLinearMethod (self):
        # F-Linear Hooke's law; T-Non-linear St. Vernant-Kirchhoff material model.
        return self.cfg['CALMODE'].getboolean('iNonLinearMethod')
    def iXDMFFileExport (self):
        # F-The HDF5 File Export function closed; T-The HDF5 File Export function opened.
        return self.cfg['CALMODE'].getboolean('iXDMFFileExport')
    def iLoadXML (self):
        # F-Load mesh from HDF5 file; T-Load mesh from XML file (when iMeshLoad = T).
        return self.cfg['CALMODE'].getboolean('iLoadXML')
    def iInteractiveMeshShow (self):
        # F-Do not show the generated mesh; T-Show the generated mesh interactively.
        return self.cfg['CALMODE'].getboolean('iInteractiveMeshShow')
    def iXDMFMeshExport (self):
        # F-The HDF5 Mesh Export function closed; T-The HDF5 Mesh Export function opened (when iHDF5FileExport = T).
        return self.cfg['CALMODE'].getboolean('iXDMFMeshExport')
    def iHDF5SubdomainsExport (self):
        # F-The HDF5 Subdomains Export function closed; T-The HDF5 Subdomains Export function opened (when iHDF5FileExport = T).
        return self.cfg['CALMODE'].getboolean('iHDF5SubdomainsExport')
    def iHDF5BoundariesExport (self):
        # F-The HDF5 Boundaries Export function closed; T-The HDF5 Boundaries Export function opened (when iHDF5FileExport = T).
        return self.cfg['CALMODE'].getboolean('iHDF5BoundariesExport')
    def iSubdomainsImport (self):
        # F-The Subdomains Import function closed; T-The Subdomains Import function opened.
        return self.cfg['CALMODE'].getboolean('iSubdomainsImport')
    def iBoundariesImport (self):
        # F-The Boundaries Import function closed; T-The Boundaries Import function opened.
        return self.cfg['CALMODE'].getboolean('iBoundariesImport')
    def iExporttxt (self):
        # F-The txt export of time list and max displacement closed; T-The txt export of time list and max displacement opened.
        return self.cfg['CALMODE'].getboolean('iExporttxt')
    def iNonUniTraction (self):
        # F-Apply uniform traction force; T-Apply non-uniform traction force.
        return self.cfg['CALMODE'].getboolean('iNonUniTraction')
    def iGravForce (self):
        # F-The gravitational force not included; T-The gravitational force included.
        return self.cfg['CALMODE'].getboolean('iGravForce')

        #===========================================
        #%% Solid Mesh numbers input
        #===========================================

    def XMesh (self):
        # cell numbers along the length of the beam, valid when iMeshLoad=False (integer) [-]
        return int(self.cfg['MESH']['XMesh'])
    def YMesh (self):
        # cell numbers along the width of the beam, valid when iMeshLoad=False (integer) [-]
        return int(self.cfg['MESH']['YMesh'])
    def ZMesh (self):
        # cell numbers along the thick of the beam, valid when iMeshLoad=False (integer) [-]
        return int(self.cfg['MESH']['ZMesh'])
