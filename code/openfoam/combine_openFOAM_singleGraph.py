#!/usr/bin/env python 
"""
# Name:
#   combine_openFOAM_singleGraph.py
#
# Purpose:
#   To collect the singleGraph files and combine them into a single, comprehensive file.  
#   File will be CSV format. 
#
# Syntax:
#   python combine_openFOAM_singleGraph.py 
#
# Modification History:
#   2021/07/26 - Joe O'Brien <joseph.r.obrien@und.edu>:  Created.
#   2021/08/25 - Joe O'Brien <joseph.r.obrien@und.edu>:  Updated for new singleGraph files (rho, Ma, total(p))
#
"""
import datetime 
import sys 
import numpy as np
import glob 
import calendar 

#--------------------------------------------------------------------------------------------------------------
# A) Input Parameters
#--------------------------------------------------------------------------------------------------------------

# Define syntax. 
def help_message():
    print('\n')
    print('Syntax: combine_openFOAM_singleGraph.py  <file1,file2,file3...> syntax\n\n')
    print(' Input:                                    ')
    print('   syntax: extentsion of the files that want to be combined (e.g. *.xy)')
    print(' Keywords:                                 ')
    print('   file1,file2,file3 - Specific files to combine. ')
    print('  Example: python combine_openFOAM_singleGraph.py .xy\n')

# Check input parameters. 
for param in sys.argv:
    if param.startswith('-h') | param.startswith('-help') | param.startswith('--help') | param.startswith('--h'):
        help_message()
        exit()

# Check to see if there are inputed files.
# If not, search for all Summary files  
if (len(sys.argv) >= 3):   # Theres input files 
    files=[]
    for file in sys.argv[1:]:
        files.append(file) 
else:     		   # No input file. Search for all *.oracles within the cwd
    if (len(sys.argv) < 3):
        syntax = sys.argv[-1]
        files = glob.glob('*'+syntax)

# Sort the files. Sorts by TAS.  
files.sort()

# Define a count for initialization of the output file. 
i = 0 

# Status Output
print("Combining singleGraph files .... > \n")

#-------------------------------------------------------------------------------------------------------------------
# B) Initalize Variables
#-------------------------------------------------------------------------------------------------------------------

# Define blank lists to hold all the data. 
case         = []                             # Geometry case 
version      = []                             # Version # for specific case
tas_init     = []                             # Freestream TAS [m/s]
aoa_init     = []                             # Angle of Attack of the Freestream Velocity 
press_init   = []                             # Freestream Static Pressure [Pa]
tempC_init   = []                             # Freestream Static Temperature [C] 
#-----------------------------------------------------------------------------------------
vel_xx_init  = []                             # Freestream Velocity (x-component) [m/s]
vel_xy_init  = []                             # Freestream Velocity (y-component) [m/s]
vel_xz_init  = []                             # Freestream Velocity (z-component) [m/s]
#-----------------------------------------------------------------------------------------
distance     = []                             # Case dimension (y-z plane; along x-axis) [m]
vel_xx       = []                             # Instanteous Velocity (x-component) [m/s]
vel_xy       = []                             # Instanteous Velocity (y-component) [m/s]
vel_xz       = []                             # Instanteous Velocity (z-component) [m/s]
press        = []                             # Dynamic Pressure [Pa] 
tempC        = []                             # Total Temperature [C] 
#-----------------------------------------------------------------------------------------
Ma           = []                             # Mach Number
rho          = []                             # Freestream Air Density
pressTot     = []                             # Total Pressure [Pa] 
#-----------------------------------------------------------------------------------------
tas_calcMach = []                             # Calculated True Air Speed (using Mach Number) [m/s]
tas_calcPres = []                             # Calculated True Air Speed (using total pressure) [m/s]
pressDy      = []                             # Dynamic Pressure [Pa]
probeDis     = []                             # Distance from Geometry (case specific) [m]
magU         = []                             # Velocity Magnitude [m/s]

#---------------------------------------------------------------------------------------------------------------------
# C) Read in the Data
#---------------------------------------------------------------------------------------------------------------------

# Define some comments. 
soundSpd = 340.29                            # Speed of Sound [m/s]
refTemp  = 288.15                            # Reference Temperature [K]

# Loop through each file.  
for nfile in files:
    # Let user know where we are at. 
    print(nfile) 
    # Define the date of the file.
    ncase   = nfile.split('_')[1]
    nversion = nfile.split('_')[2]
    ntas     = nfile.split('_')[3][3:]
    naoa     = nfile.split('_')[4][3:]
    npress   = nfile.split('_')[5].split('T')[0]
    ntempC   = nfile.split('_')[5].split('T')[1]
    
    # If this is the first time through, initialize the output file. 
    if i == 0:
        # Define output file name.  
        output = ncase+'_'+nversion+'_combined_singleGraph'+syntax 
        # Open the output file. Set up to APPEND. 
        fileout = open(output, "a")
        # Write the author
        fileout.write("OBrien, Joseph R.\n")
        # Write the date created.
        date = datetime.datetime.now()
        fileout.write(str(date.day)+'/'+str(date.month)+'/'+str(date.year)+'\n')
        # Write the softare. 
        fileout.write('openFOAM - version 7 - singleGraph postProcessing ouput\n')
        # Write the geometry. 
        fileout.write(ncase+'\n')
        # Write the Version of the geometry. 
        fileout.write(nversion+'\n')
        # Define the column header. 
        header = ['Geometry','Version#','AOA','TAS','Temp_0','Press_0','Vel0_xx','Vel0_xy','Vel0_xz',\
                  'Distance','Vel_xx','Vel_xy','Vel_xz','Press','\tTemp','Mach#','Density','Press_Tot',\
                  'MagU','probeDistance','staticPress','TAS_Mach','TAS_Press\n']
        fileout.write('\t'.join(header))
        # Define the units. 
        units  = ['N/A', '\t#', 'Degree', 'm/s', 'Celsius', 'Pa', 'm/s', 'm/s', 'm/s', 'meters', \
                  '\tm/s', 'm/s', 'm/s', 'Pa', '\tCelsius', '#', 'kg/m3', 'Pa', 'm/s', 'meters', \
                  'Pa', 'm/s', 'm/s\n']
        ##fileout.write('\t'.join(units))
    
    # Determine the number of primary variables.
    if len(nfile.split('_')) >= 8:                        #pressure and temperature file.
        print('yes', '\t->', nfile)
        # Read/Write the selected header information. 
        filein = open(nfile)
        for j, line in enumerate(filein):
            # Send all the data to blank lists. 
            case.append(ncase)
            version.append(nversion)
            tas_init.append(ntas)
            aoa_init.append(naoa)
            press_init.append(npress)
            tempC_init.append(ntempC)
            # Add 2 for UNIX time and Date. 
            # Split the line on tab, strip out newline and define parameters. 
            distance.append(line.strip().split('\t')[0])        # velocity file.
            vel_xx.append('-999')                               # not in this file - set MVC
            vel_xy.append('-999')                               # not in this file - set MVC
            vel_xz.append('-999')                               # not in this file - set MVC
            press.append(line.strip().split('\t')[3])           # static pressure [Pa]
            tempC.append(line.strip().split('\t')[4])           # temperature [celcius]
            # Check to see if the file has Ma, rho or total pressure. 
            # add mvc if it does not. 
            if 'rho' in nfile:
                rho.append(line.strip().split('\t')[5])             # Density of Air [kg/m3]
            else:
                rho.append('-999')
            if 'Ma' in nfile:
                Ma.append(line.strip().split('\t')[6])              # Mach Number [#]
            else:
                Ma.append('-999')
            if 'totalP' in nfile:
                pressTot.append(line.strip().split('\t')[7])        # total pressure [Pa]
            else:
                pressTot.append('-999')
            # Define initial velocity components based off of the file header.
            if ntas == '100':
                if naoa == '0':
                    vel_xx_init.append(99.76)
                    vel_xy_init.append(0.00)
                    vel_xz_init.append(0.00)
                elif naoa == 'Pos4':
                    vel_xx_init.append(99.76)
                    vel_xy_init.append(6.98)
                    vel_xz_init.append(0.00)
                else:
                    vel_xx_init.append(99.76)
                    vel_xy_init.append(-6.98)
                    vel_xz_init.append(0.00)
            elif ntas == '120':
                if naoa == '0':
                    vel_xx_init.append(120.00)
                    vel_xy_init.append(0.00)
                    vel_xz_init.append(0.00)
                elif naoa == 'Pos4':
                    vel_xx_init.append(120.293)
                    vel_xy_init.append(8.3912)
                    vel_xz_init.append(0.00)
                else:
                    vel_xx_init.append(120.293)
                    vel_xy_init.append(-8.3912)
                    vel_xz_init.append(0.00)
            elif ntas == '140':
                if naoa == '0':
                    vel_xx_init.append(140.00)
                    vel_xy_init.append(0.00)
                    vel_xz_init.append(0.00)
                elif naoa == 'Pos4':
                    vel_xx_init.append(139.66)
                    vel_xy_init.append(9.77)
                    vel_xz_init.append(0.00)
                else:
                    vel_xx_init.append(139.66)
                    vel_xy_init.append(-9.77)
                    vel_xz_init.append(0.00)
            else:
                print("ERROR: Initalized TAS not found.")
                sys.exit()
    else:
        # Read/Write the selected header information. 
        filein = open(nfile)
        for j, line in enumerate(filein):
            # Send all the data to blank lists. 
            case.append(ncase)
            version.append(nversion)
            tas_init.append(ntas)
            aoa_init.append(naoa)
            press_init.append(npress)
            tempC_init.append(ntempC)
            # Add 2 for UNIX time and Date. 
            # Split the line on tab, strip out newline and define parameters. 
            distance.append(line.strip().split('\t')[0])        # velocity file. 
            vel_xx.append(line.strip().split('\t')[3])          # x-component velocity file.
            vel_xy.append(line.strip().split('\t')[4])          # y-component velocity file. 
            vel_xz.append(line.strip().split('\t')[5])          # z-component velocity file. 
            press.append('-999')                                # not in this file - set MVC
            tempC.append('-999')                                # not in this file - set MVC
            Ma.append('-999')                                   # Mach Number [#]
            rho.append('-999')                                  # Density of Air [kg/m3]
            pressTot.append('-999')                             # total pressure [Pa]
            # Calculate additional parameters.
            # Define initial velocity components based off of the file header.
            if ntas == '100':
                if naoa == '0':
                    vel_xx_init.append(99.76)
                    vel_xy_init.append(0.00)
                    vel_xz_init.append(0.00)
                elif naoa == 'Pos4':
                    vel_xx_init.append(99.76)
                    vel_xy_init.append(6.98)
                    vel_xz_init.append(0.00)
                else:
                    vel_xx_init.append(99.76)
                    vel_xy_init.append(-6.98)
                    vel_xz_init.append(0.00)
            elif ntas == '120':
                if naoa == '0':
                    vel_xx_init.append(120.00)
                    vel_xy_init.append(0.00)
                    vel_xz_init.append(0.00)
                elif naoa == 'Pos4':
                    vel_xx_init.append(120.293)
                    vel_xy_init.append(8.3912)
                    vel_xz_init.append(0.00)
                else:
                    vel_xx_init.append(120.293)
                    vel_xy_init.append(-8.3912)
                    vel_xz_init.append(0.00)
            elif ntas == '140':
                if naoa == '0':
                    vel_xx_init.append(140.00)
                    vel_xy_init.append(0.00)
                    vel_xz_init.append(0.00)
                elif naoa == 'Pos4':
                    vel_xx_init.append(139.66)
                    vel_xy_init.append(9.77)
                    vel_xz_init.append(0.00)
                else:
                    vel_xx_init.append(139.66)
                    vel_xy_init.append(-9.77)
                    vel_xz_init.append(0.00)
            else:
                print("ERROR: Initalized TAS not found.")
    # Close the file. 
    filein.close()
    # Add to the counter. 
    i += 1
                
# Calculate additional parameters. 
for i in range(len(case)):
    # Define the case in order to determine geometry location. 
    if   case[i] == 'pmsCanister':
        probe_location = [-0.3926, 0.3926] # along the x-dimension [meters]
        # Distance from probe.
        if float(distance[i]) <= probe_location[0]:
            probeDis.append(float(distance[i])-probe_location[0])    # meters
        else:
            probeDis.append(float(distance[i])-probe_location[1])    # meters
    elif case[i] == 'navyPylon':
        probe_location = [14.63, 14.63] # along the x-dimension [meters]
        # Distance from probe. 
        if float(distance[i]) <= probe_location[0]:
            probeDis.append(float(distance[i])-probe_location[0])    # meters
        else:
            if float(distance[i]) >= probe_location[1]:
                probeDis.append(float(distance[i])-probe_location[1])    # meters
            else:
                probeDis.append(0.0)    # meters
    elif case[i] == 'extendedPylon':
        probe_location = [13.95, 13.95] # along the x-dimension [meters]
        # Distance from probe. 
        if float(distance[i]) <= probe_location[0]:
            probeDis.append(float(distance[i])-probe_location[0])    # meters
        else:
            if (float(distance[i]) >= probe_location[1]):
                probeDis.append(float(distance[i])-probe_location[1])    # meters
            else:
                probeDis.append(0.0)    # meters
    else:
        # Distance from probe. 
        probeDis.append(float(distance[i]))    # meters
    
    #---------------------------------
    # Calculate additional parameters.
    #---------------------------------
    # Calculate magnitude of the velocity. 
    if (float(vel_xx[i]) != -999) and (float(vel_xy[i]) != -999) and (float(vel_xz[i]) != -999):
        magU.append(np.sqrt((float(vel_xx[i])**2)+(float(vel_xy[i])**2)+(float(vel_xz[i])**2)))
    else:
        magU.append(-999)
    
    # Calculate dynamic pressure [Pa] (dynamic = total - static). 
    if (float(press[i]) != -999) and (float(pressTot[i]) != -999):
        pressDy.append(float(pressTot[i])-float(press[i]))
    else:
        pressDy.append(-999)

    # Calculate true air speed from the mach number. 
    if (float(Ma[i]) != -999) and (float(tempC[i]) != -999):
        tas_calcMach.append(soundSpd*float(Ma[i])*np.sqrt((float(tempC[i]))/refTemp))
    else:
        tas_calcMach.append(-999)

    # Calculate true air speed using the dynamic and static pressure. 
    if (float(tempC[i]) != -999) and (float(press[i]) != -999) and (float(pressDy[i]) != -999):
        tas_calcPres.append(soundSpd*np.sqrt(((5*(float(tempC[i])))/refTemp)*((((float(pressDy[i])/float(press[i]))+1)**(2/7))-1)))
    else:
        tas_calcPres.append(-999)

#-----------------------------------------------------------------------------------------
# D) Write out the data
#-----------------------------------------------------------------------------------------

for i in range(len(case)):
    fileout.write(case[i]+'\t'+version[i]+'\t\t'+aoa_init[i]+'\t'+tas_init[i]+'\t'+tempC_init[i]+'\t'+press_init[i]+'\t'\
            +str(vel_xx_init[i])+'\t'+str(vel_xy_init[i])+'\t'+str(vel_xz_init[i])+'\t'+distance[i]+'\t'+vel_xx[i]+'\t'\
            +vel_xy[i]+'\t'+vel_xz[i]+'\t'+press[i]+'\t'+tempC[i]+'\t'+str(Ma[i])+'\t'+str(rho[i])+'\t'+str(pressTot[i])+\
            '\t'+str(magU[i])+'\t'+str(probeDis[i])+'\t'+str(pressDy[i])+'\t'+str(tas_calcMach[i])+'\t'+str(tas_calcPres[i])+'\n')

# Print Status
print("\n")
print("DONE\n")

# END OF PROGRAM 
