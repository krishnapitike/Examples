import numpy as np
import glob
import re
from collections import Counter

natoms=5
ry2ev=13.605662285137
rypbr2evpang=25.71104309541616
br2ang=ry2ev/rypbr2evpang

def qe2poscar(finp):
    lattice=[]
    coords=[]
    species=[]
    
    f=open(finp, 'r')
    for line in f:
        if re.search('nat', line):
            natoms=int(line.split()[2])
            if line == None:
                nat=NAN
    f.close()
    
    f=open(finp, 'r')
    for line in f:
        if re.search("CELL_PARAMETERS", line):
            if "angstrom" in line:
                factor=1
            if "bohr" in line:
                factor=br2ang
            it=0
            for line in f:
                if (it<3):
                    line=line.split()
                    lattice.append(line)
                it=it+1
            if line == None:
                lattice=NAN
    f.close()
    f=open(finp, 'r')
    for line in f:
        if re.search("ATOMIC_POSITIONS", line):
            if "crystal" in line:
                vaspName="DIRECT"
            if "angstrom" in line:
                vaspName="CARTESIAN"
            it=0
            for line in f:
                if(it<natoms):
                    line=line.split()
                    species.append(line[0])
                    coords.append(line[1:4])
                it=it+1
            if line == None:
                coords=NAN
                species=NAN
    f.close()
    
    elements=list(Counter(species).keys()) # equals to list(set(words))
    nElements=list(Counter(species).values()) # counts the elements' frequency
    lattice=np.array(lattice).astype('float')
    coords=np.array(coords).astype('float')
    
    fvaspname=finp[:-4]+'.POSCAR.vasp'
    fvasp=open(fvaspname,'w')
    for iline in range(8+natoms):
        if(iline==0):
            outline="Written by Krishna\n"
        if(iline==1):
            outline="1.0\n"
        if(iline>=2 and iline<5):
            outline="{:23.16f} {:21.16f} {:21.16f}\n".format(lattice[iline-2,0],lattice[iline-2,1],lattice[iline-2,2])
        if(iline==5):
            outline="{}   {}   {}\n".format(elements[0],elements[1],elements[2])
        if(iline==6):
            outline="{}   {}   {}\n".format(nElements[0],nElements[1],nElements[2])
        if(iline==7):
            outline=vaspName+"\n"
        if(iline>=8):
            outline="{:21.16f} {:21.16f} {:21.16f}\n".format(coords[iline-8,0],coords[iline-8,1],coords[iline-8,2])
        fvasp.write(outline)
    return(natoms)

def qe2contcar(fout,natoms):
    ##yet to write the script
    return()

topDir='./'
for inp in glob.glob(topDir+"*.inp"):
    out=inp[:-4]+'.out'
    natom=qe2poscar(inp)
    qe2contcar(out,natoms)
print("done")
