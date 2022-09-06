#/usr/bin/env python 
"""
# NAME:
#   openfoam_plotProbes.py
#
# PURPOSE:
#   To read in the openFOAM7 probes temperature, pressure, and velocity data.  
#
# SYNTAX:
#   python openfoam_plotProbes.py openfoam_file_P openfoam_file_U openfoam_file_T
#
# INPUT:
#   1) openfoam_file  - openfoam probe pressure data *.xy. 
#
# Execution Example:
#   Linux example: python openfoam_plotProbes.py p U T
#
# Modification History:
#   2021/07/16 - Joe O'Brien <joseph.r.obrien@und.edu>:  Created by modifiying convert_FreeCAD.py 
#
# Notes: 
#
# Copyright 2016, 2017, 2018, 2019, 2020 David Delene
#
# This program is distributed under the terms of the GNU General Public License
#
# This file is part of Airborne Data Processing and Analysis (ADPAA).
#
# ADPAA is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ADPAA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ADPAA.  If not, see <http://www.gnu.org/licenses/>.
"""
from datetime import date

import matplotlib.pyplot as plt 
import sys 
import numpy as np
import time 
import re
import os

# Define syntax. 
def help_message():
  print('\n')
  print(' SYNTAX: openfoam_plotProbes.py openfoam_press openfoam_vel openfoam_temp\n')
  print('  INPUT:                                                    ')
  print('       openfoam_press - openfoam pressure file (*.xy)')
  print('       openfoam_vel   - openfoam U velocity file (*.xy)')
  print('       openfoam_temp  - openfoam T velocity file (*.xy)')
  print('  EXAMPLE: python openfoam_plotProbes.py p U T \n')

#--------------------------------------------------------------------------------------------------------------
# A) Input Arguements 
#--------------------------------------------------------------------------------------------------------------

# Check input parameters. 
for param in sys.argv:
  if param.startswith('-h') | param.startswith('-help') | param.startswith('--help') | param.startswith('--h'):
    help_message()
    exit()
# Check to make sure correct number of input parameters are sent. 
if (len(sys.argv) > 1):
  pfile=sys.argv[-3]
  ufile=sys.argv[-2]
  tfile=sys.argv[-1]
  print('openFOAM Pressure    File: ', pfile) 
  print('openFOAM Velocity    File: ', ufile)
  print('openFOAM Temperature File: ', tfile)
else:
  help_message()
  exit()

#-------------------------------------------------------------------------------------------------------------
# B) Read in the File
#-------------------------------------------------------------------------------------------------------------
# Create the read object. 
pobj = open(pfile,"r")
uobj = open(ufile,"r")
tobj = open(tfile,"r")

# Read in the file. 
plines = pobj.readlines()
ulines = uobj.readlines()
tlines = tobj.readlines()

# Define blank lists.
distanceX = []  # Note both files have the same x-axis. 
distanceY = []  # Note both files have the same x-axis. 
distanceZ = []  # Note both files have the same x-axis. 

vel_Ux    = [] # x-direction velocity - [m/s]
vel_Uy    = [] # y-direction velocity - [m/s]
vel_Uz    = [] # z-direction velocity - [m/s]

# Loop over the data and determine the probe locations. 
for i in range(len(plines)):
    if plines[i][0:7] == '# Probe':
        # read the probe distances. 
        distanceX.append(float(plines[i].strip().split('(')[-1].split()[0]))
        distanceY.append(float(plines[i].strip().split('(')[-1].split()[1]))
        distanceZ.append(float(plines[i].strip().split('(')[-1].split()[2].split(')')[0]))
        
# read the pressure data. 
press = plines[-1].strip().split()[1:]
tempK = tlines[-1].strip().split()[1:]
vel   = ulines[-1].strip().split()[1:]

# define a counter
count = 0
# loop through the velocity and strip out components
for i in range(len(vel)):
    if vel[i][0] == '(':
        vel_Ux.append(vel[i][1:])
    elif vel[i][-1] == ')':
        print('yes',vel[i])
        vel_Uz.append(vel[i][0:-1])
    else:
        vel_Uy.append(vel[i])


#-------------------------------------------------------------------------------------------------------------
# C) Define the Case
#-------------------------------------------------------------------------------------------------------------

# Define the present working directory
pwd = os.getcwd()

# Define the case. 
ncase = pwd.split('/')[-4]

# Define the environmental parameters. 
model    = ncase.split('_')[0]
version  = ncase.split('_')[1]
tas      = ncase.split('_')[2][3:]
aoa      = ncase.split('_')[3][3:]
envPress = ncase.split('_')[4].split('T')[0]
envT     = ncase.split('_')[4].split('T')[1]

# Check to see what the initial wind components will be based off of the AOA
if aoa == '0':
   ux_corr = 120
   uy_corr = 0
   uz_corr = 0
elif aoa == 'Neg4':
    ux_corr = 120
    uy_corr = -8.3912
    uz_corr = 0
elif aoa == 'Pos4':
    ux_corr = 120
    uy_corr = 8.3912
    uz_corr = 0
else:
    ux_corr = 120
    uy_corr = 0
    uz_corr = 0

#-------------------------------------------------------------------------------------------------------------
# D) Plot the Data
#-------------------------------------------------------------------------------------------------------------

# Define the figure. 
fig = plt.figure()

#-------------------------------------
# Velocity
#-------------------------------------

# Create the velocity plot. 
ax1 = fig.add_subplot(221)   
# Plot the velocity data.
ax1.plot(np.array(distanceX,dtype=float),np.array(vel_Ux,dtype=float), label='Ux')
ax1.set_xlim([4,0])
ax1.set_xlabel('Distance from Sampling Location [m]')
ax1.set_ylabel(r'Velocity [m/s]')
# Define a second x-axis. 
ax1B = ax1.twinx()
# Plot the velocity data.
ax1B.plot(np.array(distanceX,dtype=float),(np.array(vel_Uy,dtype=float)),'C1', label='Uy')
ax1B.plot(np.array(distanceX,dtype=float),(np.array(vel_Uz,dtype=float)),'g', label='Uz')
ax1B.set_ylabel(r'Velocity [m/s]')
legend = ax1B.legend(loc='upper left')
legend = ax1.legend(loc='lower right')

# Create the differential velocity plot.  
ax2 = fig.add_subplot(222)   
# Plot the velocity data
ax2.plot(np.array(distanceX,dtype=float),(np.array(vel_Ux,dtype=float)-ux_corr)*100, label='Ux')
ax2.plot(np.array(distanceX,dtype=float),(np.array(vel_Uy,dtype=float)-uy_corr)*100, label='Uy')
ax2.plot(np.array(distanceX,dtype=float),(np.array(vel_Uz,dtype=float)-uz_corr)*100, label='Uz')
ax2.set_xlabel('Distance from Sampling Location [m]')
ax2.set_ylabel(r'Differential Velocity [%]')
legend = ax2.legend(loc='upper left')
ax2.set_xlim([4,0])

# adjust subplot parameters
plt.subplots_adjust(wspace=0.35)

#-------------------------------------
# Pressure
#-------------------------------------

# Create the pressure plot. 
ax3 = fig.add_subplot(223)   
# Plot the velocity data.
ax3.plot(np.array(distanceX,dtype=float),np.array(press,dtype=float))
ax3.set_xlim([4,0])
ax3.set_xlabel('Distance from Sampling Location [m]')
ax3.set_ylabel(r'Kinematic Pressure [Pa]')
# Define a second x-axis. 
ax3B = ax3.twinx()
# Plot the temperature data.
ax3B.plot(np.array(distanceX,dtype=float),(np.array(tempK,dtype=float)),'C1', label='TempK')
ax3B.set_ylabel(r'Total Temperature [K]')
##legend = ax1.legend(loc='lower right')

# Create the differential pressure plot.  
ax4 = fig.add_subplot(224)   
# Plot the velocity data
ax4.plot(np.array(distanceX,dtype=float),(np.array(press,dtype=float)-(float(envPress)*100))*100)
ax4.set_xlabel('Distance from Sampling Location [m]')
ax4.set_ylabel(r'Differential Kinematic Pressure [%]')
##legend = ax4.legend(loc='upper left')
ax4.set_xlim([4,0])
# Define a second x-axis. 
ax4B = ax4.twinx()
# Plot the temperature data.
ax4B.plot(np.array(distanceX,dtype=float),(np.array(tempK,dtype=float)-(float(envT)+273))*100,'C1', label='TempK')
ax4B.set_ylabel(r'Differential Temperature [%]')

# Define a title
plt.suptitle(model+' (Version='+version+', TAS='+tas+', AOA='+aoa+', Pressure='+envPress+', Temperature='+envT+')')

# END OF PROGRAM
