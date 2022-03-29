##########################################################
#                                                        # 
# Script to transform a cluster file with the format     #
# needed by the Diffusion module into a xyz file         #
#                                                        #
#                                                        #
# python cluster_to_xyz.py file.dat                      #
#                                                        #
##########################################################

import sys

# input file
inputfile = str(sys.argv[1])
outputfile = 'cluster.xyz'
cfile = open( inputfile,'r' )
xyzfile = open(outputfile,'w')

# number of atoms
line = cfile.readline()
column = line.split()
nmolecules = int(column[0])
natoms = 3*nmolecules        # water molecules
line = str(nmolecules)+'\n'
xyzfile.write(line)
line = '\n'
xyzfile.write(line)

for imol in range(nmolecules):

    line = cfile.readline()
    column = line.split()
    line = 'O     '+str(column[0])+'   '+str(column[1])+'   '+str(column[2])+'\n'
    xyzfile.write(line)

    line = cfile.readline()
    column = line.split()
    line = 'H     '+str(column[0])+'   '+str(column[1])+'   '+str(column[2])+'\n'
    xyzfile.write(line)

    line = cfile.readline()
    column = line.split()
    line = 'H     '+str(column[0])+'   '+str(column[1])+'   '+str(column[2])+'\n'
    xyzfile.write(line)    
    
xyzfile.close()
cfile.close()
