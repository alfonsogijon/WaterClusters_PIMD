##########################################################
#                                                        # 
# Script to transform a xyz file in a file with the      #
# format needed by the Diffusion module                  #
#                                                        #
#                                                        #
# python xyt_to_cluster.py file.xyz                      #
#                                                        #
##########################################################

import sys

# input file
inputfile = str(sys.argv[1])
outputfile = 'cluster.dat'
xyzfile = open( inputfile,'r' )
cfile = open(outputfile,'w')

# number of atoms
line = xyzfile.readline()
column = line.split()
natoms = int(column[0])
nmolecules = int(natoms/3)
next(xyzfile)
line = str(nmolecules)+'\n'
cfile.write(line)

for iat in range(natoms):
    line = xyzfile.readline()
    column = line.split()
    line = str(column[1])+'   '+str(column[2])+'   '+str(column[3])+'\n'
    cfile.write(line)
    
xyzfile.close()
cfile.close()
