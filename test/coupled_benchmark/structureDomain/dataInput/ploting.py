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
for line in open('structureDomain/structureResults/tip-displacementY_0.txt', 'r'):
    lines = [i for i in line.split()]
    y.append(float(lines[0]))

for line in open('structureDomain/dataInput/Slone_et_al.txt', 'r'):
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
plt.savefig('result_compare.png')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FILE END  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#