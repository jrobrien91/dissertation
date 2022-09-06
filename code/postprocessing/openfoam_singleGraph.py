#/usr/bin/env python 
"""
# NAME:
#   openfoam_singleGraph.py
#
# PURPOSE:
#   To read in the openFOAM7 singleGraph data.  
#
# SYNTAX:
#   python openfoam_singleGraph.py openfoam_file_PT openfoam_file_U 
#
# INPUT:
#   1) openfoam_file_PT  - openfoam singleGraph pressure and temperature data *.xy. 
#   2) openfoam_file_PT  - openfoam singleGraph velocity data *.xy. 
#
# Execution Example:
#   Linux example: python openfoam_plotPress.py line_p_T.xy line_u.xy
#
# Modification History:
#   2020/07/21 - Joe O'Brien <joseph.r.obrien@und.edu>:  Created by modifiying convert_FreeCAD.py 
#
# Notes: 
#
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
  print(' SYNTAX: openfoam_singleGraph.py openfoam_press openfoam_vel\n')
  print('  INPUT:                                                    ')
  print('       openfoam_press - openfoam pressure and temperature file (*.xy)')
  print('       openfoam_vel   - openfoam velocity file (*.xy)')
  print('  EXAMPLE: python openfoam_singleGraph.py line_p_T.xy line_U.xy \n')

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
  pfile=sys.argv[-2]
  ufile=sys.argv[-1]
  print('openFOAM Pressure/Temp File: ', pfile) 
  print('openFOAM Velocity File: ', ufile)
else:
  help_message()
  exit()

#-------------------------------------------------------------------------------------------------------------
# B) Read in the File
#-------------------------------------------------------------------------------------------------------------

# Create the read object. 
pobj = open(pfile,"r")
uobj = open(ufile,"r")

# Read in the file. 
plines = pobj.readlines()
ulines = uobj.readlines()

# Define blank lists.
distance = []  # Note both files have the same x-axis. 
press    = []  # in Pascal
tempK    = []  # Total Temperature in Kelvin
vel_xx   = []  # Instanteous Velocity - x component [m/s]
vel_xy   = []  # Instanteous Velocity - y component [m/s]
vel_xz   = []  # Instanteous Velocity - z component [m/s]

# Loop over the data and fill the arrays. 
for i in range(len(plines)):
    # read the pressure distances. 
    distance.append(float(plines[i].strip().split('\t')[0]))
    press.append(float(plines[i].strip().split('\t')[1]))
    tempK.append(float(plines[i].strip().split('\t')[2]))
    # read in the velocities. 
    vel_xx.append(float(ulines[i].strip().split('\t')[1]))
    vel_xy.append(float(ulines[i].strip().split('\t')[2]))
    vel_xz.append(float(ulines[i].strip().split('\t')[3]))

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
   uy_corr = 1
   uz_corr = 1
elif aoa == 'Neg4':
    ux_corr = 120
    uy_corr = -8.3912
    uz_corr = 1
elif aoa == 'Pos4':
    ux_corr = 120
    uy_corr = 8.3912
    uz_corr = 1
else:
    ux_corr = 120
    uy_corr = 1
    uz_corr = 1

#----------------------------------------------------------------------------------------------------------
# D) Calculate some Parameters
#----------------------------------------------------------------------------------------------------------

vel_mag = []
for i in range(len(vel_xx)):
    vel_mag.append(np.sqrt((vel_xx[i]**2)+(vel_xy[i]**2)+(vel_xz[i]**2)))

#----------------------------------------------------------------------------------------------------------
# E) Plot the Data 
#----------------------------------------------------------------------------------------------------------

# Define the figure. 
fig = plt.figure()

#-------------------------------------
# Velocity
#-------------------------------------

# Create the velocity plot. 
ax1 = fig.add_subplot(321)
# Plot the velocity data.
ax1.plot(np.array(distance),np.array(vel_xx), label='Ux')
ax1.set_xlim([-1,1])
ax1.set_xlabel('Distance from Sampling Location [m]')
ax1.set_ylabel(r'Velocity [m/s]', color='b')
# Define a second x-axis. 
ax1B = ax1.twinx()
# Plot the velocity data.
ax1B.plot(np.array(distance),np.array(vel_xy),'C1', label='Uy')
ax1B.plot(np.array(distance),np.array(vel_xz),'g', label='Uz')
ax1B.set_ylabel(r'Velocity [m/s]', color='C1')
legend = ax1B.legend(loc='upper left')
legend = ax1.legend(loc='lower left')

# Create the velocity plot. 
axA = fig.add_subplot(322)
# Plot the velocity data.
axA.plot(np.array(distance),np.array(vel_mag)/ux_corr, label='U/Uo')
axA.set_xlim([-1,1])
axA.set_xlabel('Distance from Sampling Location [m]')
axA.set_ylabel(r'Velocity [m/s]', color='b')
legend = axA.legend(loc='lower left')

#--------------------------------------
# Pressure
#--------------------------------------

# Create the differential velocity plot.  
ax2 = fig.add_subplot(323)
# Plot the velocity data
ax2.plot(np.array(distance),np.array(press), label='Total Pressure [Pa]')
ax2.set_xlabel('Distance from Sampling Location [m]')
ax2.set_ylabel(r'Total Pressure [Pa]', color='b')
#legend = ax2.legend(loc='upper left')
ax2.set_xlim([-1,1])

# Create the differential velocity plot.  
ax3 = fig.add_subplot(324)
# Plot the velocity data
ax3.plot(np.array(distance),np.array(press)/(float(envPress)*100), label='Total Pressure [Pa]')
ax3.set_xlabel('Distance from Sampling Location [m]')
ax3.set_ylabel(r'Total Pressure [Pa]', color='b')
#legend = ax2.legend(loc='upper left')
ax3.set_xlim([-1,1])

#--------------------------------------
# Temperature
#--------------------------------------

# Create the differential velocity plot.  
ax4 = fig.add_subplot(325)
# Plot the velocity data
ax4.plot(np.array(distance),np.array(tempK), label='Total Temperature [K]')
ax4.set_xlabel('Distance from Sampling Location [m]')
ax4.set_ylabel(r'Total Temperature [K]', color='b')
#legend = ax2.legend(loc='upper left')
ax4.set_xlim([-1,1])

# Create the differential velocity plot.  
ax5 = fig.add_subplot(326)
# Plot the velocity data
ax5.plot(np.array(distance),np.array(tempK)/(float(envT)+273), label='Total Temperature [K]')
ax5.set_xlabel('Distance from Sampling Location [m]')
ax5.set_ylabel(r'Total Temperature [K]', color='b')
#legend = ax2.legend(loc='upper left')
ax5.set_xlim([-1,1])

# Define a title
plt.suptitle(model+' (Version='+version+', TAS='+tas+', AOA='+aoa+', Pressure='+envPress+', Temperature='+envT+')')

# END OF PROGRAM
