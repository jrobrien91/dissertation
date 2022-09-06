#/usr/bin/env python 
"""
# NAME:
#   openfoam_singleGraph_Combined.py
#
# PURPOSE:
#   To read in the openFOAM7 singleGraph data.  
#
# SYNTAX:
#   python openfoam_singleGraph.py openfoam_Combined_file
#
# INPUT:
#   1) openfoam_Combiend_file  - openfoam singleGraph combined file. 
#
# Execution Example:
#   Linux example: python openfoam_plotPress.py line_p_T.xy line_u.xy
#
# Modification History:
#   2020/08/02 - Joe O'Brien <joseph.r.obrien@und.edu>:  Created by modifiying openfoam_singleGraph.py 
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
  print(' SYNTAX: openfoam_singleGraph_Combined.py openfoam_singleGraph_Combined\n')
  print('  EXAMPLE: python openfoam_singleGraph_Combined.py pmsCanister_v7_combined_singleGraph.xy \n')

#----------------------------------------------------------------------------------------------------------------------
# A) Input Arguements 
#----------------------------------------------------------------------------------------------------------------------

# Check input parameters. 
for param in sys.argv:
  if param.startswith('-h') | param.startswith('-help') | param.startswith('--help') | param.startswith('--h'):
    help_message()
    exit()
# Check to make sure correct number of input parameters are sent. 
if (len(sys.argv) > 1):
  nfile=sys.argv[-1]
  print('openFOAM File: ', nfile) 
else:
  help_message()
  exit()

#---------------------------------------------------------------------------------------------------------------------
# B) Open the File
#---------------------------------------------------------------------------------------------------------------------

# Create the read object. 
nobj = open(nfile,"r")

# Read in the file. 
nlines = nobj.readlines()

#---------------------------------------------------------------------------------------------------------------------
# C) Initalize Variables
#---------------------------------------------------------------------------------------------------------------------

# Define blank lists to hold all the data. 
case         = []                             # Geometry case 
version      = []                             # Version # for specific case
tas_init     = []                             # Freestream TAS [m/s]
aoa_init     = []                             # Angle of Attack of the Freestream Velocity 
press_init   = []                             # Freestream Total Pressure [Pa]
tempC_init   = []                             # Freestream Total Temperature [C] 
#-----------------------------------------------------------------------------------------
vel_xx_init  = []                             # Freestream Velocity (x-component) [m/s]
vel_xy_init  = []                             # Freestream Velocity (y-component) [m/s]
vel_xz_init  = []                             # Freestream Velocity (z-component) [m/s]
#-----------------------------------------------------------------------------------------
distance     = []                             # Case dimension (y-z plane; along x-axis) [m]
vel_xx       = []                             # Instanteous Velocity (x-component) [m/s]
vel_xy       = []                             # Instanteous Velocity (y-component) [m/s]
vel_xz       = []                             # Instanteous Velocity (z-component) [m/s]
press        = []                             # Freestream Static Pressure [Pa] 
tempC        = []                             # Freestream Temperature [C] 
Mach         = []                             # Mach Number [#]
rho          = []                             # Air Density [g/m3]
pressTot     = []                             # Total Pressure [Pa]
#-----------------------------------------------------------------------------------------
probeDis     = []                             # Distance from Geometry (case specific) [m]
magU         = []                             # Velocity Magnitude [m/s]
pressDynam   = []                             # Dynamic Pressure based on total pressure calculation [Pa]
tas_calcMach = []                             # True Air Speed based on Mach Number [m/s]
tas_calcPres = []                             # True Air Speed based on Dynamic/Static Pressure [m/s]

#--------------------------------------------------------------------------------------------------------------------
# D) Read in the Data
#--------------------------------------------------------------------------------------------------------------------

# Loop over the data and fill the arrays. 
for i in range(9,len(nlines)):
    # Read the info. 
    case.append(nlines[i].strip().split('\t')[0])
    version.append(nlines[i].strip().split('\t')[1])
    # Figure out the AOA. 
    tempAOA = nlines[i].strip().split('\t')[3]
    if tempAOA == '0':
        aoa_init.append(0)
    elif tempAOA == 'Pos4':
        aoa_init.append(4)
    elif tempAOA == 'Neg4':
        aoa_init.append(-4)
    else:
        print('ERROR: AOA not found')
        sys.exit()
    # Read the initalized values
    tas_init.append(float(nlines[i].strip().split('\t')[4]))
    tempC_init.append(float(nlines[i].strip().split('\t')[5]))
    press_init.append(float(nlines[i].strip().split('\t')[6]))
    vel_xx_init.append(float(nlines[i].strip().split('\t')[7]))
    vel_xy_init.append(float(nlines[i].strip().split('\t')[8]))
    vel_xz_init.append(float(nlines[i].strip().split('\t')[9]))
    # read the pressure distances. 
    distance.append(float(nlines[i].strip().split('\t')[10]))
    vel_xx.append(float(nlines[i].strip().split('\t')[11]))
    vel_xy.append(float(nlines[i].strip().split('\t')[12]))
    vel_xz.append(float(nlines[i].strip().split('\t')[13]))
    press.append(float(nlines[i].strip().split('\t')[14]))
    tempC.append(float(nlines[i].strip().split('\t')[15]))
    Mach.append(float(nlines[i].strip().split('\t')[16]))
    rho.append(float(nlines[i].strip().split('\t')[17]))
    pressTot.append(float(nlines[i].strip().split('\t')[18]))
    # calculated parameters
    magU.append(float(nlines[i].strip().split('\t')[19]))
    #if len(nlines[i].strip().split('\t')[20]) > 0:
    probeDis.append(float(nlines[i].strip().split('\t')[20]))
    #else:
        #probeDis.append(float(nlines[i].strip().split('\t')[21]))
    pressDynam.append(float(nlines[i].strip().split('\t')[21]))
    tas_calcMach.append(float(nlines[i].strip().split('\t')[22]))
    tas_calcPres.append(float(nlines[i].strip().split('\t')[23]))
# Close the file. 
nobj.close()

# Flip to numpy arrays.
case         = np.array(case)
tas_init     = np.ma.masked_less(tas_init, -998)
tempC_init   = np.ma.masked_less(tempC_init, -998)
press_init   = np.ma.masked_less(press_init, -998)
aoa_init     = np.ma.masked_less(aoa_init, -998)
distance     = np.ma.masked_less(distance, -998)
#------------------------------
vel_xx       = np.ma.masked_less(vel_xx, -998)
vel_xy       = np.ma.masked_less(vel_xy, -998)
vel_xz       = np.ma.masked_less(vel_xz, -998)
#-----------------------------
press        = np.ma.masked_less(press, -998)
tempC        = np.ma.masked_less(tempC, -998)
magU         = np.ma.masked_less(magU, -998)
Mach         = np.ma.masked_less(Mach, -998)
rho          = np.ma.masked_less(rho, -998)
probeDis     = np.ma.masked_less(probeDis, -998)
pressTot     = np.ma.masked_less(pressTot, -998)
pressDynam   = np.ma.masked_less(pressDynam, -998)
tas_calcMach = np.ma.masked_less(tas_calcMach, -998)
tas_calcPres = np.ma.masked_less(tas_calcPres, -998)
#---------------------------------------------------------------------------------------------------------------------
# E) Determine Cases
#---------------------------------------------------------------------------------------------------------------------

# Seperate Cases based on Pressure/Temperature combination. 
spress900   = np.where(press_init == 900)
spress800   = np.where(press_init == 800)
# Seperate Pressures based on TAS 
tas100     = np.where(tas_init == 100)
tas120     = np.where(tas_init == 120)
tas140     = np.where(tas_init == 140)
# Seperate Based on AOA 
aoa0       = np.where(aoa_init == 0)
aoa4       = np.where(aoa_init == 4)
aoaN4      = np.where(aoa_init == -4)

# Find the cases. Determine pressure - airspeed combinations. 
tas100_spress900   = np.intersect1d(tas100, spress900)
tas120_spress900   = np.intersect1d(tas120, spress900)
tas140_spress900   = np.intersect1d(tas140, spress900)
#--------------------------------------------------------------
tas100_spress800   = np.intersect1d(tas100, spress800)
tas120_spress800   = np.intersect1d(tas120, spress800)
tas140_spress800   = np.intersect1d(tas140, spress800)
#--------------------------------------------------------------
tas120_aoa0        = np.intersect1d(tas120, aoa0)
tas120_aoaPos4     = np.intersect1d(tas120, aoa4)
tas120_aoaNeg4     = np.intersect1d(tas120, aoaN4)
#--------------------------------------------------------------
tas140_aoa0        = np.intersect1d(tas140, aoa0)
tas140_aoaPos4     = np.intersect1d(tas140, aoa4)
tas140_aoaNeg4     = np.intersect1d(tas140, aoaN4)

#----------------------------------------------------------------------------------------------------------
# E) Plot the Data 
#----------------------------------------------------------------------------------------------------------
"""
parms     = ['AOA_0', 'AOA_4', 'AOA_N4'] 
perPress9 = [tas100_spress900, tas120_spress900, tas140_spress900]
perPress8 = [tas100_spress800, tas120_spress800, tas140_spress800]

# Iterate over the three AOA cases. 
for i in range(len(parms)):
    # Define the figure
    fig = plt.figure()
    # Define the title of the plot.
    fig.suptitle(case[0]+' '+version[0]+ ' ('+parms[i]+')')

    #------------------------------------------------------------------------------------------------------
    # Velocity
    #-------------------------------------------------------------------------------------------------------

    #---------------------------
    # 800 mb, 20C Environmental
    #---------------------------

    # Create the velocity plot. 
    axA = fig.add_subplot(321)
    for j in range(len(perPress8)):
        # Find the AOA subset.
        if parms[i] == 'AOA_0':
            aoa_x      = np.where(aoa_init == 0)
            aoa_subset = np.intersect1d(perPress8[j],aoa_x)
        elif parms[i] == 'AOA_4':
            aoa_x      = np.where(aoa_init == 4)
            aoa_subset = np.intersect1d(perPress8[j],aoa_x)
        elif parms[i] == 'AOA_N4':
            aoa_x      = np.where(aoa_init == -4)
            aoa_subset = np.intersect1d(perPress8[j],aoa_x)
        # Plot the velocity data.
        if j == 0:
            axA.plot(probeDis[aoa_subset],magU[aoa_subset]/100, label='TAS=100')
        elif j == 1:
            axA.plot(probeDis[aoa_subset],magU[aoa_subset]/120, label='TAS=120')
        else:
            axA.plot(probeDis[aoa_subset],magU[aoa_subset]/140, label='TAS=140')
    axA.set_xlim([-0.5,0.5])
    axA.set_xlabel('Distance from Sampling Location [m]')
    axA.set_ylabel(r'Velocity [m/s] - (U/Uo)', color='b')
    axA.set_title('800 mb, 20C Freestream Environemental Values')
    legend = axA.legend(loc='lower left')
    
    #-----------------------------
    # 900 mb, 30C Environmental
    #-----------------------------

    # Create the velocity plot. 
    ax1 = fig.add_subplot(322)
    for j in range(len(perPress9)):
        # Find the AOA subset.
        if parms[i] == 'AOA_0':
            aoa_x      = np.where(aoa_init == 0)
            aoa_subset = np.intersect1d(perPress9[j],aoa_x)
        elif parms[i] == 'AOA_4':
            aoa_x      = np.where(aoa_init == 4)
            aoa_subset = np.intersect1d(perPress9[j],aoa_x)
        elif parms[i] == 'AOA_N4':
            aoa_x      = np.where(aoa_init == -4)
            aoa_subset = np.intersect1d(perPress9[j],aoa_x)
        # Plot the velocity data.
        if j == 0:
            ax1.plot(probeDis[aoa_subset],magU[aoa_subset]/100, label='TAS=100')
        elif j == 1:
            ax1.plot(probeDis[aoa_subset],magU[aoa_subset]/120, label='TAS=120')
        else:
            ax1.plot(probeDis[aoa_subset],magU[aoa_subset]/140, label='TAS=140')
    ax1.set_xlim([-0.5,0.5])
    ax1.set_xlabel('Distance from Sampling Location [m]')
    ax1.set_ylabel(r'Velocity [m/s] - (U/Uo)', color='b')
    ax1.set_title("900 mb, 30C Freestream Environmental Values")
    legend = ax1.legend(loc='lower left')

    #-----------------------------------------------------------------------------------------------------
    # Pressure
    #-----------------------------------------------------------------------------------------------------

    #----------------------------
    # 800 mb, 20C Environmental
    #----------------------------

    # Create the differential velocity plot.  
    axB = fig.add_subplot(323)
    for j in range(len(perPress8)):
        # Find the AOA subset.
        if parms[i] == 'AOA_0':
            aoa_x      = np.where(aoa_init == 0)
            aoa_subset = np.intersect1d(perPress8[j],aoa_x)
        elif parms[i] == 'AOA_4':
            aoa_x      = np.where(aoa_init == 4)
            aoa_subset = np.intersect1d(perPress8[j],aoa_x)
        elif parms[i] == 'AOA_N4':
            aoa_x      = np.where(aoa_init == -4)
            aoa_subset = np.intersect1d(perPress8[j],aoa_x)
        # Plot the velocity data.
        if j == 0:
            axB.plot(probeDis[aoa_subset],press[aoa_subset]/(800*100), label='TAS=100')
        elif j == 1:
            axB.plot(probeDis[aoa_subset],press[aoa_subset]/(800*100), label='TAS=120')
        else:
            axB.plot(probeDis[aoa_subset],press[aoa_subset]/(800*100), label='TAS=140')
    axB.set_xlabel('Distance from Sampling Location [m]')
    axB.set_ylabel(r'Total Pressure [Pa] - (Pt/Po)', color='b')
    legend = axB.legend(loc='upper left')
    axB.set_xlim([-0.5,0.5])
    
    #------------------------------
    # 900 mb, 30C Environmental
    #------------------------------

    # Create the differential velocity plot.  
    ax2 = fig.add_subplot(324)
    for j in range(len(perPress9)):
        # Find the AOA subset.
        if parms[i] == 'AOA_0':
            aoa_x      = np.where(aoa_init == 0)
            aoa_subset = np.intersect1d(perPress9[j],aoa_x)
        elif parms[i] == 'AOA_4':
            aoa_x      = np.where(aoa_init == 4)
            aoa_subset = np.intersect1d(perPress9[j],aoa_x)
        elif parms[i] == 'AOA_N4':
            aoa_x      = np.where(aoa_init == -4)
            aoa_subset = np.intersect1d(perPress9[j],aoa_x)
        # Plot the velocity data.
        if j == 0:
            ax2.plot(probeDis[aoa_subset],press[aoa_subset]/(900*100), label='TAS=100')
        elif j == 1:
            ax2.plot(probeDis[aoa_subset],press[aoa_subset]/(900*100), label='TAS=120')
        else:
            ax2.plot(probeDis[aoa_subset],press[aoa_subset]/(900*100), label='TAS=140')
    ax2.set_xlabel('Distance from Sampling Location [m]')
    ax2.set_ylabel(r'Total Pressure [Pa] - (Pt/Po)', color='b')
    legend = ax2.legend(loc='upper left')
    ax2.set_xlim([-0.5,0.5])


    #---------------------------------------------------------------------------------------------------
    # Temperature
    #----------------------------------------------------------------------------------------------------

    #-------------------------------
    # 800 mb, 20C Environmental 
    #-------------------------------
    # Create the differential velocity plot.  
    axC = fig.add_subplot(325)
    for j in range(len(perPress8)):
        # Find the AOA subset.
        if parms[i] == 'AOA_0':
            aoa_x      = np.where(aoa_init == 0)
            aoa_subset = np.intersect1d(perPress8[j],aoa_x)
        elif parms[i] == 'AOA_4':
            aoa_x      = np.where(aoa_init == 4)
            aoa_subset = np.intersect1d(perPress8[j],aoa_x)
        elif parms[i] == 'AOA_N4':
            aoa_x      = np.where(aoa_init == -4)
            aoa_subset = np.intersect1d(perPress8[j],aoa_x)
        # Plot the velocity data.
        if j == 0:
            axC.plot(probeDis[aoa_subset],tempC[aoa_subset]/(293), label='TAS=100')
        elif j == 1:
            axC.plot(probeDis[aoa_subset],tempC[aoa_subset]/(293), label='TAS=120')
        else:
            axC.plot(probeDis[aoa_subset],tempC[aoa_subset]/(293), label='TAS=140')
    axC.set_xlabel('Distance from Sampling Location [m]')
    axC.set_ylabel(r'Total Temperature [K] - (T/To)', color='b')
    legend = axC.legend(loc='upper left')
    axC.set_xlim([-0.5,0.5])
    
    #-------------------------------
    # 900 mb, 30C Environmental 
    #-------------------------------
    # Create the differential velocity plot.  
    ax3 = fig.add_subplot(326)
    for j in range(len(perPress9)):
        # Find the AOA subset.
        if parms[i] == 'AOA_0':
            aoa_x      = np.where(aoa_init == 0)
            aoa_subset = np.intersect1d(perPress9[j],aoa_x)
        elif parms[i] == 'AOA_4':
            aoa_x      = np.where(aoa_init == 4)
            aoa_subset = np.intersect1d(perPress9[j],aoa_x)
        elif parms[i] == 'AOA_N4':
            aoa_x      = np.where(aoa_init == -4)
            aoa_subset = np.intersect1d(perPress9[j],aoa_x)
        # Plot the velocity data.
        if j == 0:
            ax3.plot(probeDis[aoa_subset],tempC[aoa_subset]/(303), label='TAS=100')
        elif j == 1:
            ax3.plot(probeDis[aoa_subset],tempC[aoa_subset]/(303), label='TAS=120')
        else:
            ax3.plot(probeDis[aoa_subset],tempC[aoa_subset]/(303), label='TAS=140')
    ax3.set_xlabel('Distance from Sampling Location [m]')
    ax3.set_ylabel(r'Total Temperature [K] - (T/To)', color='b')
    legend = ax3.legend(loc='upper left')
    ax3.set_xlim([-0.5,0.5])

    # save the plot. 
    pltTitle = 'OBrien_'+case[0]+'_'+version[0]+'_'+parms[i]
    plt.savefig(pltTitle)
"""
"""
parms     = ['AOA_N4'] 
##perPress9 = [tas120_aoaNeg4]
#perPress9 = [tas120_aoa0]
perPress9 = [tas140_aoaNeg4]

# Iterate over the three AOA cases. 
for i in range(len(parms)):
    # Define the figure
    fig = plt.figure()
    # Define the title of the plot.
    #fig.suptitle('Ratio of Freestream Values through the Instrument Location')

    #------------------------------------------------------------------------------------------------------
    # Velocity
    #-------------------------------------------------------------------------------------------------------
    
    #-----------------------------
    # 900 mb, 30C Environmental
    #-----------------------------

    # Create the velocity plot. 
    ax1 = fig.add_subplot(311)
    for j in range(len(perPress9)):
        # Find the pylon subset.
        extend_x      = np.where(case == 'extendedPylon')
        extend_subset = np.intersect1d(perPress9[j],extend_x)
        navy_x      = np.where(case == 'navyPylon')
        navy_subset = np.intersect1d(perPress9[j],navy_x)
        # Find the velocity subset. 
        vel_extend = np.where(magU[extend_subset].mask != True)
        vel_navy   = np.where(magU[navy_subset].mask != True)
        # Plot the velocity data.
        if j == 0:
            ##ax1.plot(probeDis[extend_subset][vel_extend],magU[extend_subset][vel_extend]/120, label='Extended Pylon, AOA-4$^\circ$')
            ##ax1.plot(probeDis[navy_subset][vel_navy],magU[navy_subset][vel_navy]/120, label='Navy Pylon, AOA-4$^\circ$')
            ##ax1.plot(probeDis[extend_subset][vel_extend],magU[extend_subset][vel_extend]/120, label='Extended Pylon, AOA 0$^\circ$')
            ##ax1.plot(probeDis[navy_subset][vel_navy],magU[navy_subset][vel_navy]/120, label='Navy Pylon, AOA 0$^\circ$')
            ##ax1.plot(probeDis[extend_subset][vel_extend],magU[extend_subset][vel_extend]/120, label='Extended Pylon, AOA+4$^\circ$')
            ##ax1.plot(probeDis[navy_subset][vel_navy],magU[navy_subset][vel_navy]/120, label='Navy Pylon, AOA+4$^\circ$')
            #------------------------------------------------------------------------------------------------------------------------#
            ax1.plot(probeDis[extend_subset][vel_extend],magU[extend_subset][vel_extend]/140, label='Extended Pylon, AOA-4$^\circ$')
            ax1.plot(probeDis[navy_subset][vel_navy],magU[navy_subset][vel_navy]/140, label='Navy Pylon, AOA-4$^\circ$')
            ##ax1.plot(probeDis[extend_subset][vel_extend],magU[extend_subset][vel_extend]/140, label='Extended Pylon, AOA 0$^\circ$')
            ##ax1.plot(probeDis[navy_subset][vel_navy],magU[navy_subset][vel_navy]/140, label='Navy Pylon, AOA 0$^\circ$')
            ##ax1.plot(probeDis[extend_subset][vel_extend],magU[extend_subset][vel_extend]/140, label='Extended Pylon, AOA+4$^\circ$')
            ##ax1.plot(probeDis[navy_subset][vel_navy],magU[navy_subset][vel_navy]/140, label='Navy Pylon, AOA+4$^\circ$')
        #if j == 1:
        #    ax1.plot(probeDis[extend_subset][vel_extend],magU[extend_subset][vel_extend]/120, label='Extended, AOA +4')
        #    ax1.plot(probeDis[navy_subset][vel_navy],magU[navy_subset][vel_navy]/120, label='Navy, AOA +4')
        #else:
        #if j == 2:
        #    ax1.plot(probeDis[extend_subset][vel_extend],magU[extend_subset][vel_extend]/120, label='Extended, AOA -4')
        #    ax1.plot(probeDis[navy_subset][vel_navy],magU[navy_subset][vel_navy]/120, label='Navy, AOA -4')
    ax1.set_xlim([-0.5,0.5])
    ax1.set_xlabel('Distance from Canister [m]')
    ax1.set_ylabel(r'Velocity (U/Uo)', color='b')
    ##ax1.set_title("900 mb, 30C, 120 m/s Freestream Environmental Values")
    ax1.set_title("800 mb, 20C, 140 m/s Freestream Environmental Values")
    #ax1.set_title("800 mb, 20C, 140 m/s Freestream Environmental Values")
    legend = ax1.legend(loc='lower left')
    ax1.grid()

    #-----------------------------------------------------------------------------------------------------
    # Pressure
    #-----------------------------------------------------------------------------------------------------

    #------------------------------
    # 900 mb, 30C Environmental
    #------------------------------

    # Create the differential velocity plot.  
    ax2 = fig.add_subplot(312)
    for j in range(len(perPress9)):
        extend_x      = np.where(case == 'extendedPylon')
        extend_subset = np.intersect1d(perPress9[j],extend_x)
        navy_x      = np.where(case == 'navyPylon')
        navy_subset = np.intersect1d(perPress9[j],navy_x)
        # Find the velocity subset. 
        p_extend = np.where(press[extend_subset].mask != True)
        p_navy   = np.where(press[navy_subset].mask != True)
        # Plot the velocity data.
        if j == 0:
            ##ax2.plot(probeDis[extend_subset][p_extend],press[extend_subset][p_extend]/(900*100), label='Extended Pylon, AOA-4$^\circ$')
            ##ax2.plot(probeDis[navy_subset][p_navy],press[navy_subset][p_navy]/(900*100), label='Navy Pylon, AOA-4$^\circ$')
            ##ax2.plot(probeDis[extend_subset][p_extend],press[extend_subset][p_extend]/(900*100), label='Extended Pylon, AOA 0$^\circ$')
            ##ax2.plot(probeDis[navy_subset][p_navy],press[navy_subset][p_navy]/(900*100), label='Navy Pylon, AOA 0$^\circ$')
            ##ax2.plot(probeDis[extend_subset][p_extend],press[extend_subset][p_extend]/(900*100), label='Extended Pylon, AOA+4$^\circ$')
            ##ax2.plot(probeDis[navy_subset][p_navy],press[navy_subset][p_navy]/(900*100), label='Navy Pylon, AOA+4$^\circ$')
            #--------------------------------------------------------------------------------------------------------------------------#
            ax2.plot(probeDis[extend_subset][p_extend],press[extend_subset][p_extend]/(800*100), label='Extended Pylon, AOA-4$^\circ$')
            ax2.plot(probeDis[navy_subset][p_navy],press[navy_subset][p_navy]/(800*100), label='Navy Pylon, AOA-4$^\circ$')
            ##ax2.plot(probeDis[extend_subset][p_extend],press[extend_subset][p_extend]/(800*100), label='Extended Pylon, AOA 0$^\circ$')
            ##ax2.plot(probeDis[navy_subset][p_navy],press[navy_subset][p_navy]/(800*100), label='Navy Pylon, AOA 0$^\circ$')
            ##ax2.plot(probeDis[extend_subset][p_extend],press[extend_subset][p_extend]/(800*100), label='Extended Pylon, AOA+4$^\circ$')
            ##ax2.plot(probeDis[navy_subset][p_navy],press[navy_subset][p_navy]/(800*100), label='Navy Pylon, AOA+4$^\circ$')
        #if j == 1:
        #    ax2.plot(probeDis[extend_subset][p_extend],press[extend_subset][p_extend]/(900*100), label='Extended, AOA +4')
        #    ax2.plot(probeDis[navy_subset][p_navy],press[navy_subset][p_navy]/(900*100), label='Navy, AOA +4')
        #else:
        #if j == 2:
        #    ax2.plot(probeDis[extend_subset][p_extend],press[extend_subset][p_extend]/(900*100), label='Extended, AOA -4')
        #    ax2.plot(probeDis[navy_subset][p_navy],press[navy_subset][p_navy]/(900*100), label='Navy, AOA -4')
    ax2.set_xlabel('Distance from Canister [m]')
    ax2.set_ylabel(r'Pressure (Pt/Po)', color='b')
    legend = ax2.legend(loc='lower left')
    ax2.set_xlim([-0.5,0.5])
    ax2.grid()

    #---------------------------------------------------------------------------------------------------
    # Temperature
    #----------------------------------------------------------------------------------------------------

    #-------------------------------
    # 900 mb, 30C Environmental 
    #-------------------------------
    # Create the differential velocity plot.  
    ax3 = fig.add_subplot(313)
    for j in range(len(perPress9)):
        extend_x      = np.where(case == 'extendedPylon')
        extend_subset = np.intersect1d(perPress9[j],extend_x)
        navy_x      = np.where(case == 'navyPylon')
        navy_subset = np.intersect1d(perPress9[j],navy_x)
        # Find the velocity subset. 
        t_extend = np.where(tempC[extend_subset].mask != True)
        t_navy   = np.where(tempC[navy_subset].mask != True)
        # Plot the velocity data.
        if j == 0:
            ##ax3.plot(probeDis[extend_subset][t_extend],tempC[extend_subset][t_extend]/(303), label='Extended Pylon, AOA-4$^\circ$')
            ##ax3.plot(probeDis[navy_subset][t_navy],tempC[navy_subset][t_navy]/(303), label='Navy Pylon, AOA-4$^\circ$')
            ##ax3.plot(probeDis[extend_subset][t_extend],tempC[extend_subset][t_extend]/(303), label='Extended Pylon, AOA 0$^\circ$')
            ##ax3.plot(probeDis[navy_subset][t_navy],tempC[navy_subset][t_navy]/(303), label='Navy Pylon, AOA 0$^\circ$')
            ##ax3.plot(probeDis[extend_subset][t_extend],tempC[extend_subset][t_extend]/(303), label='Extended Pylon, AOA+4$^\circ$')
            ##ax3.plot(probeDis[navy_subset][t_navy],tempC[navy_subset][t_navy]/(303), label='Navy Pylon, AOA+4$^\circ$')
            #----------------------------------------------------------------------------------------------------------------------#
            ax3.plot(probeDis[extend_subset][t_extend],tempC[extend_subset][t_extend]/(293), label='Extended Pylon, AOA-4$^\circ$')
            ax3.plot(probeDis[navy_subset][t_navy],tempC[navy_subset][t_navy]/(293), label='Navy Pylon, AOA-4$^\circ$')
            ##ax3.plot(probeDis[extend_subset][t_extend],tempC[extend_subset][t_extend]/(293), label='Extended Pylon, AOA 0$^\circ$')
            ##ax3.plot(probeDis[navy_subset][t_navy],tempC[navy_subset][t_navy]/(293), label='Navy Pylon, AOA 0$^\circ$')
            ##ax3.plot(probeDis[extend_subset][t_extend],tempC[extend_subset][t_extend]/(293), label='Extended Pylon, AOA+4$^\circ$')
            #3ax3.plot(probeDis[navy_subset][t_navy],tempC[navy_subset][t_navy]/(293), label='Navy Pylon, AOA+4$^\circ$')
        #if j == 1:
        #    ax3.plot(probeDis[extend_subset][t_extend],tempC[extend_subset][t_extend]/(303), label='Extended, AOA +4')
        #    ax3.plot(probeDis[navy_subset][t_navy],tempC[navy_subset][t_navy]/(303), label='Navy, AOA +4')
        #else:
        #if j == 2:
        #    ax3.plot(probeDis[extend_subset][t_extend],tempC[extend_subset][t_extend]/(303), label='Extended, AOA -4')
        #    ax3.plot(probeDis[navy_subset][t_navy],tempC[navy_subset][t_navy]/(303), label='Navy, AOA -4')
    ax3.set_xlabel('Distance from Canister [m]')
    ax3.set_ylabel(r'Temperature (T/To)', color='b')
    legend = ax3.legend(loc='upper left')
    ax3.set_xlim([-0.5,0.5])
    ax3.grid()

    # save the plot. 
    pltTitle = 'OBrien_'+case[0]+'_'+version[0]+'_'+parms[i]
    plt.savefig(pltTitle)
"""
"""
parms     = ['AOA_0'] 
perPress9 = [tas120_aoa0, tas140_aoa0]
#perPress9 = [tas120_aoa0]
#perPress9 = [tas140_aoa0]
    
# Define the figure
fig3 = plt.figure()

#------------------------------------------------------------------------------------------------------
# Calculated OPENFOAM TAS Comparison 
#-------------------------------------------------------------------------------------------------------

# Find the pylon subset.
extend_x      = np.where(case == 'extendedPylon')
extend_subset = np.intersect1d(perPress9[0],extend_x)
navy_x        = np.where(case == 'navyPylon')
navy_subset   = np.intersect1d(perPress9[0],navy_x)
    
# Create the TAS plot. 
ax1 = fig3.add_subplot(221)
# Plot the data
#ax1.scatter(tas_calcMach[extend_subset][1:], tas_calcMach[navy_subset])
ax1.scatter(tas_calcMach[extend_subset], tas_calcPres[extend_subset])
# Label the axes
ax1.set_xlabel(r'TAS-extendedPylon (Mach-based) [m/s]')
ax1.set_ylabel(r'TAS-extendedPylon (Pressure-based) [m/s]')
# Label the Plot
ax1.set_title("FreeSteam Velocity - 120m/s, AOA 0")
# Set limits on the axes
ax1.set_xlim([10,130])
ax1.set_ylim([10,130])
# Define a 1:1 ratio line. 
newline = np.arange(1,150,1)
ax1.plot(newline,newline,'g')
# Set grid lines
ax1.grid()

# Create the TAS plot. 
ax2 = fig3.add_subplot(222)
# Plot the data
ax2.scatter(tas_calcMach[navy_subset], tas_calcPres[navy_subset])
# Label the axes
ax2.set_xlabel(r'TAS-navyPylon (Mach-based) [m/s]')
ax2.set_ylabel(r'TAS-navyPylon (Pressure-based) [m/s]')
# Label the Plot
ax2.set_title("FreeSteam Velocity - 120 m/s AOA 0")
# Set limits on the axes
ax2.set_xlim([10,130])
ax2.set_ylim([10,130])
# Define a 1:1 ratio line. 
newline = np.arange(1,150,1)
ax2.plot(newline,newline,'g')
# Set grid lines
ax2.grid()

# Find the pylon subset.
extend2_subset = np.intersect1d(perPress9[1],extend_x)
navy2_subset   = np.intersect1d(perPress9[1],navy_x)

# Create the TAS plot. 
ax3 = fig3.add_subplot(223)
# Plot the data
ax3.scatter(tas_calcMach[extend2_subset], tas_calcPres[extend2_subset])
# Label the axes
ax3.set_xlabel(r'TAS-extendedPylon (Mach-based) [m/s]')
ax3.set_ylabel(r'TAS-extendedPylon (Pressure-based) [m/s]')
# Label the Plot
ax3.set_title("FreeSteam Velocity - 140 m/s AOA 0")
# Set limits on the axes
ax3.set_xlim([50,150])
ax3.set_ylim([50,150])
# Define a 1:1 ratio line. 
newline = np.arange(1,170,1)
ax3.plot(newline,newline,'g')
# Set grid lines
ax3.grid()

# Create the TAS plot. 
ax4 = fig3.add_subplot(224)
# Plot the data
ax4.scatter(tas_calcMach[navy2_subset], tas_calcPres[navy2_subset])
# Label the axes
ax4.set_xlabel(r'TAS-navyPylon (Mach-based) [m/s]')
ax4.set_ylabel(r'TAS-navyPylon (Pressure-based) [m/s]')
# Label the Plot
ax4.set_title("FreeSteam Velocity - 140 m/s AOA 0")
# Set limits on the axes
ax4.set_xlim([50,150])
ax4.set_ylim([50,150])
# Define a 1:1 ratio line. 
newline = np.arange(1,170,1)
ax4.plot(newline,newline,'g')
# Set grid lines
ax4.grid()
"""

"""
parms     = ['AOA_0'] 
perPress9 = [tas140_aoa0,tas140_aoaPos4,tas140_aoaNeg4]

# Iterate over the three AOA cases. 
for i in range(len(parms)):
    # Define the figure
    fig = plt.figure()
    # Define the title of the plot.
    #fig.suptitle('Ratio of Freestream Values through the Instrument Location')

    #---------------------------------------------------------------------------------------------------
    # Dynamic Pressure
    #----------------------------------------------------------------------------------------------------

    #-------------------------------
    # 900 mb, 30C Environmental 
    #-------------------------------
    # Create the differential velocity plot.  
    ax1 = fig.add_subplot(311)
    extend_x      = np.where(case == 'extendedPylon')
    extend_subset = np.intersect1d(perPress9[0],extend_x)
    navy_x      = np.where(case == 'navyPylon')
    navy_subset = np.intersect1d(perPress9[0],navy_x)
    # Find the velocity subset. 
    t_extend = np.where(pressDynam[extend_subset].mask != True)
    t_navy   = np.where(pressDynam[navy_subset].mask != True)
    # Plot the velocity data.
    ax1.plot(probeDis[extend_subset][t_extend],pressDynam[extend_subset][t_extend]/100, label='Extended Pylon, AOA 0$^\circ$')
    ax1.plot(probeDis[navy_subset][t_navy],pressDynam[navy_subset][t_navy]/100, label='Navy Pylon, AOA 0$^\circ$')
    ax1.set_xlabel('Distance from Canister [m]')
    ax1.set_ylabel(r'Dynamic Pressure [hPa]', color='b')
    legend = ax1.legend(loc='lower left')
    ax1.set_xlim([-0.5,0.5])
    ax1.set_title("800 mb, 20C, 140 m/s Freestream Environmental Values")
    ax1.grid()
    
    #-------------------------------
    # 900 mb, 30C Environmental 
    #-------------------------------
    # Create the differential velocity plot.  
    ax2 = fig.add_subplot(312)
    extend_x      = np.where(case == 'extendedPylon')
    extend_subset = np.intersect1d(perPress9[1],extend_x)
    navy_x      = np.where(case == 'navyPylon')
    navy_subset = np.intersect1d(perPress9[1],navy_x)
    # Find the velocity subset. 
    t_extend = np.where(pressDynam[extend_subset].mask != True)
    t_navy   = np.where(pressDynam[navy_subset].mask != True)
    # Plot the velocity data.
    ax2.plot(probeDis[extend_subset][t_extend],pressDynam[extend_subset][t_extend]/(100), label='Extended Pylon, AOA+4$^\circ$')
    ax2.plot(probeDis[navy_subset][t_navy],pressDynam[navy_subset][t_navy]/(100), label='Navy Pylon, AOA+4$^\circ$')
    ax2.set_xlabel('Distance from Canister [m]')
    ax2.set_ylabel(r'Dynamic Pressure [hPa]', color='b')
    legend = ax2.legend(loc='lower left')
    ax2.set_xlim([-0.5,0.5])
    ax2.grid()
    
    #-------------------------------
    # 900 mb, 30C Environmental 
    #-------------------------------
    # Create the differential velocity plot.  
    ax3 = fig.add_subplot(313)
    extend_x      = np.where(case == 'extendedPylon')
    extend_subset = np.intersect1d(perPress9[1],extend_x)
    navy_x      = np.where(case == 'navyPylon')
    navy_subset = np.intersect1d(perPress9[2],navy_x)
    # Find the velocity subset. 
    t_extend = np.where(pressDynam[extend_subset].mask != True)
    t_navy   = np.where(pressDynam[navy_subset].mask != True)
    # Plot the velocity data.
    ax3.plot(probeDis[extend_subset][t_extend],pressDynam[extend_subset][t_extend]/(100), label='Extended Pylon, AOA-4$^\circ$')
    ax3.plot(probeDis[navy_subset][t_navy],pressDynam[navy_subset][t_navy]/(100), label='Navy Pylon, AOA-4$^\circ$')
    ax3.set_xlabel('Distance from Canister [m]')
    ax3.set_ylabel(r'Dynamic Pressure [hPa]', color='b')
    legend = ax3.legend(loc='lower left')
    ax3.set_xlim([-0.5,0.5])
    ax3.grid()

    # save the plot. 
    pltTitle = 'OBrien_'+case[0]+'_'+version[0]+'_'+parms[i]
    plt.savefig(pltTitle)
"""
parms     = ['AOA_0'] 
perPress9 = [tas120_aoa0,tas120_aoaPos4,tas120_aoaNeg4]

# Iterate over the three AOA cases. 
for i in range(len(parms)):
    # Define the figure
    fig = plt.figure()
    # Define the title of the plot.
    #fig.suptitle('Ratio of Freestream Values through the Instrument Location')

    #---------------------------------------------------------------------------------------------------
    # Dynamic Pressure
    #----------------------------------------------------------------------------------------------------

    #-------------------------------
    # 900 mb, 30C Environmental 
    #-------------------------------
    # Create the differential velocity plot.  
    ax1 = fig.add_subplot(311)
    extend_x      = np.where(case == 'extendedPylon')
    extend_subset = np.intersect1d(perPress9[0],extend_x)
    navy_x      = np.where(case == 'navyPylon')
    navy_subset = np.intersect1d(perPress9[0],navy_x)
    # Find the velocity subset. 
    t_extend = np.where(pressDynam[extend_subset].mask != True)
    t_navy   = np.where(pressDynam[navy_subset].mask != True)
    # Plot the velocity data.
    ax1.plot(probeDis[extend_subset][t_extend],tas_calcMach[extend_subset][t_extend], label='Extended Pylon, AOA 0$^\circ$')
    ax1.plot(probeDis[navy_subset][t_navy],tas_calcMach[navy_subset][t_navy], label='Navy Pylon, AOA 0$^\circ$')
    ax1.set_xlabel('Distance from Canister [m]')
    ax1.set_ylabel(r'True Airspeed [m/s]', color='b')
    legend = ax1.legend(loc='lower left')
    ax1.set_xlim([-0.5,0.5])
    ax1.set_title("900 mb, 30C, 120 m/s Freestream Environmental Values")
    ax1.grid()
    
    #-------------------------------
    # 900 mb, 30C Environmental 
    #-------------------------------
    # Create the differential velocity plot.  
    ax2 = fig.add_subplot(312)
    extend_x      = np.where(case == 'extendedPylon')
    extend_subset = np.intersect1d(perPress9[1],extend_x)
    navy_x      = np.where(case == 'navyPylon')
    navy_subset = np.intersect1d(perPress9[1],navy_x)
    # Find the velocity subset. 
    t_extend = np.where(tas_calcMach[extend_subset].mask != True)
    t_navy   = np.where(tas_calcMach[navy_subset].mask != True)
    # Plot the velocity data.
    ax2.plot(probeDis[extend_subset][t_extend],tas_calcMach[extend_subset][t_extend], label='Extended Pylon, AOA+4$^\circ$')
    ax2.plot(probeDis[navy_subset][t_navy],tas_calcMach[navy_subset][t_navy], label='Navy Pylon, AOA+4$^\circ$')
    ax2.set_xlabel('Distance from Canister [m]')
    ax2.set_ylabel(r'True Airspeed [m/s]', color='b')
    legend = ax2.legend(loc='lower left')
    ax2.set_xlim([-0.5,0.5])
    ax2.grid()
    
    #-------------------------------
    # 900 mb, 30C Environmental 
    #-------------------------------
    # Create the differential velocity plot.  
    ax3 = fig.add_subplot(313)
    extend_x      = np.where(case == 'extendedPylon')
    extend_subset = np.intersect1d(perPress9[1],extend_x)
    navy_x      = np.where(case == 'navyPylon')
    navy_subset = np.intersect1d(perPress9[2],navy_x)
    # Find the velocity subset. 
    t_extend = np.where(tas_calcMach[extend_subset].mask != True)
    t_navy   = np.where(tas_calcMach[navy_subset].mask != True)
    # Plot the velocity data.
    ax3.plot(probeDis[extend_subset][t_extend],tas_calcMach[extend_subset][t_extend], label='Extended Pylon, AOA-4$^\circ$')
    ax3.plot(probeDis[navy_subset][t_navy],tas_calcMach[navy_subset][t_navy], label='Navy Pylon, AOA-4$^\circ$')
    ax3.set_xlabel('Distance from Canister [m]')
    ax3.set_ylabel(r'True Airspeed [m/s]', color='b')
    legend = ax3.legend(loc='lower left')
    ax3.set_xlim([-0.5,0.5])
    ax3.grid()

    # save the plot. 
    pltTitle = 'OBrien_'+case[0]+'_'+version[0]+'_'+parms[i]
    plt.savefig(pltTitle)

