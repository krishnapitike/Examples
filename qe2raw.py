import numpy as np
import csv
import glob
import re

natoms=5
ry2ev=13.605662285137
rypbr2evpang=25.71104309541616

def getNatoms(file):
    f=open(file, 'r')
    for line in f:
        if re.search('nat', line):
            nat=int(line.split()[2])
            if line == None:
                nat=NAN
    f.close()
    return(nat)

def getLattice(file):
    lattice=[]
    coords=[]
    species=[]
    f=open(file, 'r')
    for line in f:
        if re.search("angstrom", line):
            it=0
            for line in f:
                if (it<3):
                    line=line.split()
                    lattice.append(line)
                elif(it>=4 and it<4+natoms):
                    line=line.split()
                    species.append(line[0])
                    coords.append(line[1:4])
                it=it+1
            lattice=np.array(lattice).astype('float')
            coords=np.array(coords).astype('float')
            if line == None:
                lattice=NAN
                coords=NAN
                species=NAN
    f.close()
    return(species,lattice,coords)

def getEnergy(file):
    f=open(file, 'r')
    for line in f:
        if re.search('! ', line):
            toten=float(line.split()[4])*ry2ev
            if line == None:
                toten=NAN
    f.close()
    return(toten)

def getStress(file):
    stress=[]
    f=open(file, 'r')
    for line in f:
        if re.search("P=", line):
            it=0
            for line in f:
                if (it<3):
                    line=line.split()
                    stress.append(line[3:6])
                it=it+1
            stress=np.array(stress).astype('float')
            if line == None:
                stress=NAN
    f.close()
    return(stress)

def getForces(file):
    forces=[]
    f=open(file, 'r')
    for line in f:
        if re.search("Ry/au", line):
            it=0
            for line in f:
                if (it>0 and it<=natoms):
                    line=line.split()
                    forces.append(line[6:9])
                it=it+1
            forces=np.array(forces).astype('float')*rypbr2evpang
            if line == None:
                forces=NAN
    f.close()
    return(forces)

box=open('box.raw','w')
boxWriter=csv.writer(box, delimiter =' ')
coord=open('coord.raw','w')
coordWriter=csv.writer(coord, delimiter =' ')
energy=open('energy.raw','w')
energyWriter=csv.writer(energy, delimiter =' ')
force=open('force.raw','w')
forceWriter=csv.writer(force, delimiter =' ')
virial=open('virial.raw','w')
virialWriter=csv.writer(virial, delimiter =' ')

topDir='./'
for out in glob.glob(topDir+"*.out"):
    inp=out[:-4]+'.inp'
    print(out)
    #inp=topDir+str(i)+'.inp'
    #out=topDir+str(i)+'.out'
    natoms=getNatoms(inp)
    species,lattice,coords=getLattice(inp)
    toten=[getEnergy(out)]
    stress=getStress(out)
    forces=getForces(out)
    energyWriter.writerow( '{:0.8e}'.format(x) for x in toten)
    boxWriter.writerow( '{:0.8e}'.format(x) for x in lattice.reshape(9))
    virialWriter.writerow( '{:0.8e}'.format(x) for x in stress.reshape(9))
    coordWriter.writerow( '{:0.8e}'.format(x) for x in coords.reshape(natoms*3))
    forceWriter.writerow( '{:0.8e}'.format(x) for x in forces.reshape(natoms*3))
box.close()
coord.close()
energy.close()
force.close()
virial.close()
