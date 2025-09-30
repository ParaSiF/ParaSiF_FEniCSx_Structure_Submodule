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

    @file timeMarching.py
    
    @author W. Liu
    
    @brief time marching file of the structure code.

"""

#_________________________________________________________________________________________
#
#%% Import packages
#_________________________________________________________________________________________
import math

class timeMarching:

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #%% Time marching parameters define
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def Time_Marching_Parameters(self):
        if self.iContinueRun():
            hdf5checkpointDataInTemp = HDF5File(self.LOCAL_COMM_WORLD, self.inputFolderPath + "/checkpointData.h5", "r")
            # Read start time [s]
            self.Start_Time = self.dt() + hdf5checkpointDataInTemp.attributes("/ud/vector_0")["timestamp"]
            # Calculate time steps [-]
            self.Time_Steps = math.ceil((self.T() - self.Start_Time)/self.dt())
            # Close file and delete HDF5File object
            hdf5checkpointDataInTemp.close()
            del hdf5checkpointDataInTemp
        else:
            if self.iResetStartTime():
                # Reset start time [s]
                self.Start_Time = self.dt() + self.newStartTime()
                # Calculate time steps [-]
                self.Time_Steps = math.ceil((self.T() - self.Start_Time)/self.dt())
            else:
                # Set start time [s]
                self.Start_Time = self.dt()
                # Calculate time steps [-]
                self.Time_Steps = math.ceil(self.T()/self.dt())
        # Initialise sub-iterations counter
        self.Start_Number_Sub_Iteration = 1

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #%% Define Generalised-alpha method functions
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # Define the acceleration at the present time step
    def Acceleration_March_Term_One(self,
                                    displacement_function,
                                    displacement_previous_function,
                                    velocity_previous_function):
        return (2 * (displacement_function - displacement_previous_function -
                     (self.dt() * velocity_previous_function))/
                     (self.dt()**2))

    def Acceleration_March_Term_Two(self,
                                    acceleration_previous_function,
                                    beta_gam):
        return ((1 - (2 * beta_gam)) * acceleration_previous_function)

    def Acceleration_March_Term_Three(self, beta_gam):
        return (1 / (2 * beta_gam))

    def AMCK(self,
             displacement_function,
             displacement_previous_function,
             velocity_previous_function,
             acceleration_previous_function,
             beta_gam):
        return (self.Acceleration_March_Term_Three(beta_gam) *
                    (self.Acceleration_March_Term_One(displacement_function,
                    displacement_previous_function,
                    velocity_previous_function) -
                    self.Acceleration_March_Term_Two(acceleration_previous_function,beta_gam)))

    # Define the velocity at the present time step
    def Velocity_March_Term_One(self,
                                acceleration_previous_function,
                                gamma_gam):
        return ((1 - gamma_gam) * acceleration_previous_function * self.dt())

    def Velocity_March_Term_Two(self,
                                acceleration_function,
                                gamma_gam):
        return (acceleration_function * gamma_gam * self.dt())

    def UMCK(self,
             acceleration_function,
             velocity_previous_function,
             acceleration_previous_function,
             gamma_gam):
        return (self.Velocity_March_Term_One(acceleration_previous_function, gamma_gam) +
                self.Velocity_March_Term_Two(acceleration_function, gamma_gam) +
                velocity_previous_function)

    # define the calculation of intermediate averages based on generalized alpha method
    def Generalized_Alpha_Weights(self,
                                  present_function,
                                  previous_function,
                                  weights):
        return (weights * previous_function + (1 - weights) * present_function)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FILE END  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#