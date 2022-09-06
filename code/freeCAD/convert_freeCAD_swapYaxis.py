#!/usr/bin/env python 
"""
# NAME:
#   convert_freeCAD_mm2meters.py
#
# PURPOSE:
#   To read in the FreeCAD waveobject formatted file (*.obj) and convert all coordinates from mm to meters.  
#
# SYNTAX:
#   python convert_freeCAD_mm2meters.py freeCAD_file 
#
# INPUT:
#   1) freeCAD_file  - Waveobject geometery file *obj. 
#
# Execution Example:
#   Linux example: python convert_freeCAD_mm2meters.py pmsCanister_v3.obj
#
# Modification History:
#   2020/09/25 - Joe O'Brien <joseph.r.obrien@und.edu>:  Created by modifiying convert_freeCAD_mm2meters.py.
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

import sys 
import numpy as np
import time 
import re

# Define syntax. 
def help_message():
  print('\n')
  print(' SYNTAX: convert_freeCAD_mm2meters.py freeCAD_file\n')
  print('  INPUT:                                                    ')
  print('       freeCAD_file - FreeCAD Waveobject geometry file (*.obj) - mm units')
  print('  EXAMPLE: python convert_freeCAD_mm2meters.py pmsCanister_v3.obj \n')

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
  nfile=sys.argv[-1]
  print('FreeCAD File:', nfile) 
else:
  help_message()
  exit()

#-------------------------------------------------------------------------------------------------------------
# B) Read in the File
#-------------------------------------------------------------------------------------------------------------

# Create the read/write object. 
nobj = open(nfile,"r+")

# Read in the file. 
nlines = nobj.readlines()

# Loop over the lines within the file and convert values from mm to m. 
for i in range(len(nlines)):
    # Check the leading character for comments. 
    if (nlines[i][0] == 'v'):
        print(nlines[i])
        # Strip out the newline feed. 
        nlines[i].strip()
        # Convert the coordinates from mm to m. 
        x = float(nlines[i].split()[1])
        y = float(nlines[i].split()[2])
        if y > 0:
            y = y * (-1)
        else:
            y = abs(y)
        z = float(nlines[i].split()[3])
        # Assemble the new string. 
        nlines[i] = nlines[i].split(' ')[0]+' '+str(x)+' '+str(y)+' '+str(z)+'\n'
        print(nlines[i])

#----------------------------------------------------------------------------------------------------------
# Write out the data
#----------------------------------------------------------------------------------------------------------

# Move the pointer back to the beginning of the file. 
nobj.seek(0)

# Write out the data. 
for i in range(len(nlines)):
    nobj.write(nlines[i])

# Trunacate and close the file. 
nobj.truncate()
nobj.close


# END OF PROGRAM
