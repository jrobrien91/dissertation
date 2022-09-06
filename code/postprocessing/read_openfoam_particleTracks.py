#!/usr/bin/env python
"""
# NAME:
#   read_openfoam_particleTracks.py
#
# PURPOSE:
#   To read in the openFOAM particle tracks *csv file. 
#
# SYNTAX:
#   python read_openfoam_particleTracks.pyi csv_file 
#
#   Keywords:
#      csv_file  - openFOAM lagrangian particle track file witin the comma separated value format. 
#
# EXECUTION EXAMPLE:
#   Linux example: python read_openfoam_particleTracks.py pmsCanister_v8_tas120_aoa0_900T33_particleTracks.csv
#
# MODIFICATION HISTORY:
#   2022/01/25 - Joe O'Brien <joseph.r.obrien@und.edu> : Created 
#
# PROGRAMMING GUIDES: 
#
# NOTES:
#
"""
from timesync  import time_sync

import time 
import sys
import numpy as np
import matplotlib.pyplot as plt 
import datetime 

#----------------------------------------------------------------------------------------------
# I) Define Functions
#----------------------------------------------------------------------------------------------

def help_message():
  print('\n')
  print('Syntax: read_openfoam_particleTracks.py <-h> csv_file\n')
  print('      INPUT:  ')
  print('        csv_file:  - openFOAM lagrangian particle track file with particle characteristics\n') 
  print('    OPTIONS:  ')
  print('               -h:  - Help statement. Print Syntax')
  print('    EXAMPLE: python read_HSRL_APR3_merge.py pmsCanister_v8_tas120_aoa0_900T33_particleTracks.csv \n')


#----------------------------------------------------------------------------------------------
# II) Input Agruments 
#----------------------------------------------------------------------------------------------

# Define the starting time of the code. 
t0 = time.time()

# Check for options within the system arguments
for param in sys.argv:
    if param.startswith('-h'):
        help_message()
        exit()

# Check to make sure there was an input file.
if (len(sys.argv) < 2):
    help_message()
    exit()
else:
  nfile = sys.argv[-1]
  print("\n")
  print("LPT CSV FILE: ", nfile)

#-----------------------------------------------------------------------------------------------
# III) Read In The File
#-----------------------------------------------------------------------------------------------

# Open the file with read option
obj    = open(nfile,'r')
rows = obj.readlines()

# Define tags or keys to create a dictionary
tags = rows[0].strip().split(',')

# Read in the entire file and transpose to form columns
columns = np.loadtxt(rows[1:],delimiter=',').T

# Create a dictionary 
data = dict(zip(tags,columns))

#-----------------------------------------------------------------------------------------------
# IV) Sample Volume Determination
#-----------------------------------------------------------------------------------------------

# Define a new dictionary to hold particles within 
tempx = ['X_pos','Y_pos','Z_pos','size','age','u_vel','v_vel','w_vel']
tempy = [[],[],[],[],[],[],[],[]]
inSample = dict(zip(tempx,tempy))

# Vertices for the cubic sample volume
xmin = -0.35
xmax =  0.35
ymin = -0.06
ymax =  0.06
zmin = -0.60
zmax =  0.60

# Loop over all the particles
for i in range(len(data['d'][:])):
    # Compare against all vertices of the cubic domain
    if (data['Points:0'][i] > xmin) and  (data['Points:0'][i] < xmax):
        if (data['Points:1'][i] > ymin) and  (data['Points:1'][i] < ymax):
            if (data['Points:2'][i] > zmin) and  (data['Points:2'][i] < zmin):
                  print('yes')
                  inSample['X_pos'].append(data['Points:0'][i])
                  inSample['Y_pos'].append(data['Points:1'][i])
                  inSample['Z_pos'].append(data['Points:2'][i])
                  inSample['size'].append(data['size'][i])
                  inSample['age'].append(data['age'][i])
                  inSample['u_vel'].append(data['u_vel'][i])
                  inSample['v_vel'].append(data['v_vel'][i])
                  inSample['z_vel'].append(data['w_vel'][i])

#-----------------------------------------------------------------------------------------------
# END OF PROGRAM 
#-----------------------------------------------------------------------------------------------

# Define the ending time of the program
t1 = time.time()

# Display the timing of the program.
tyme = t1-t0
print('')
if tyme > 60:
  print((tyme/60.), 'min ', tyme, 'sec')
else:
  print(tyme, 'sec')

