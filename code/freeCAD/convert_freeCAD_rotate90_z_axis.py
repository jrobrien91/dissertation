#!/usr/bin/env python 
"""
# NAME:
#   convert_freeCAD_rotate90_z_axis.py
#
# PURPOSE:
#   To read in the FreeCAD waveobject formatted file (*.obj) and rotate the geometory 90deg about the z-axis.  
#
# SYNTAX:
#   python convert_freeCAD_rotate90_ccw.py freeCAD_file 
#
# INPUT:
#   1) freeCAD_file  - Waveobject geometery file *obj. 
#
# Execution Example:
#   Linux example: python convert_freeCAD_rotate90_ccw.py pmsCanister_v4.obj
#
# Modification History:
#   2021/06/28 - Joe O'Brien <joseph.r.obrien@und.edu>:  Created by modifiying convert_freeCAD_mm2meters.py.
#
# Notes: 
#
"""
from datetime import date

import sys 
import numpy as np
import time 
import re

# Define syntax. 
def help_message():
  print('\n')
  print(' SYNTAX: convert_freeCAD_rotate90_ccw.py freeCAD_file\n')
  print('  INPUT:                                                    ')
  print('       freeCAD_file - FreeCAD Waveobject geometry file (*.obj) - mm units')
  print('  EXAMPLE: python convert_freeCAD_rotate90_ccw.py pmsCanister_v4.obj \n')

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
    if (nlines[i][0] == 'v') and (nlines[i][0:2] != 'vt'):
        if (nlines[i][0:2] != 'vn'):
            print("")
            print(nlines[i])
            # Strip out the newline feed. 
            nlines[i].strip()
            # Split the line. 
            ntmp = nlines[i].split()
            # Define x,y coordinates and translate y. 
            x = float(ntmp[1])
            y = float(ntmp[2])
            z = float(ntmp[3])
            # Assemble the new string.
            ##nlines[i] = str(ntmp[0])+' '+str(x)+' '+str((-1)*z)+' '+str(y)+' '+str(ntmp[4])+' '+str(ntmp[5])+\
            ##        ' '+str(ntmp[6])+'\n'
            nlines[i] = str(ntmp[0])+' '+str(z)+' '+str(y)+' '+str((-1)*x)+' '+str(ntmp[4])+' '+str(ntmp[5])+\
                    ' '+str(ntmp[6])+'\n'
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
