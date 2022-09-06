#/usr/bin/env python 
"""
# NAME:
#   create_particlePos_openfoam.py
#
# PURPOSE:
#   To create particle positions for input into openFOAM.  
#
# SYNTAX:
#   python create_particlePos_openfoam.py 
#
# INPUT:
#
# Execution Example:
#   Linux example: python create_particlePos_openfoam.py
#
# Modification History:
#   2021/12/14 - Joe O'Brien <joseph.r.obrien@und.edu>:  Created by modifiying convert_FreeCAD.py 
#
# Notes: 
#
# Copyright 2016, 2017, 2018, 2019, 2020, 2021 David Delene
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
  print(' SYNTAX: create_particlePos_openfoam.py [nParticle=#]\n')
  print('   OPTIONS:                                                    ')
  print('       nParticle = number of particles within the 0.5m initial block [#]')
  print('  EXAMPLE: python create_particlePos_openfoam.py \n')

#---------------------------------------------------------------------------------------------------------------------------
# A) Input Arguements 
#---------------------------------------------------------------------------------------------------------------------------

nParticle = 10              # initial number of desired particles. 

# Check input parameters. 
for param in sys.argv:
  if param.startswith('-h') | param.startswith('-help') | param.startswith('--help') | param.startswith('--h'):
    help_message()
    exit()
  if param.startswith('nParticle=') | param.startswith('nParticle=') | param.startswith('--help') | param.startswith('--h'):
    nParticle = int(sys.argv[-1].split('=')[-1])

#----------------------------------------------------------------------------------------------------------------------------
# B) Create a homogeneous cloud 
#----------------------------------------------------------------------------------------------------------------------------

# create a blank list to hold the 
posX = []
posY = []
posZ = []

# Define a blank counter 
j = 0
k = 0

# Loop over the number of particles 
for i in range(nParticle):
        # Define the coordinate position. Initial is (-10, 0.0000, 0.0)
        # Define the x-coordinate inlet boundary. Will always be the same.
        posX.append(-10.00)
        
        # Define the y-coordinate. Add a mm in initial configuration change per particle. 
        ##tempY = 2.55 + (j*0.001)   # mm particle spacing
        tempY = -0.50 + (j*0.01) # 100 micron particle spacing
        # Define the z-coordinate. Add 500 microns to the initial configuration per particle. 
        ##tempZ = 13.65 + (k*0.001)   # mm particle spacing
        ##tempZ =  -1.0+ (k*0.005) # 100 micron particle spacing 
        tempZ = -0.5 + (k*0.005)
        if tempY < 0.5:
            posY.append(np.around(tempY, 5))
            posZ.append(np.around(tempZ, 5))
            j += 1
        else:
            posY.append(0.5)
            posZ.append(np.around(tempZ, 4))
            k += 1
            j = 0

#-----------------------------------------------------------------------------------------------------------------------------
# C) Write out the particles
#-----------------------------------------------------------------------------------------------------------------------------

nfile = open('particlePositions.txt','w')

# Write out the header data
nfile.write("FoamFile\n")
nfile.write("{\n")
nfile.write("\tversion\t\t2.0;\n")
nfile.write("\tformat\t\tascii;\n")
nfile.write("\tclass\t\tvectorField;\n")
nfile.write("\tlocation\tconstant;\n")
nfile.write("\tobject\t\tpositions;\n")
nfile.write("}\n")
nfile.write("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n")
nfile.write("(\n")

# Write out the particle postiions
for i in range(len(posX)):
    temp = '( '+str(posX[i])+' '+str(posY[i])+' '+str(posZ[i])+')\n'
    nfile.write(temp)

# Close the partheses and file. 
nfile.write(")")
nfile.close()

# END OF PROGRAM
