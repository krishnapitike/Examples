import numpy as np
import glob

def readvasp(vasp):
    fvasp=open(vasp,'r')
    lat=[]
    coord=[]
    elementArray=[]
    iline=0
    for line in fvasp:
        if iline==1: multFactor=float(line.split()[0])
        if iline>=2 and iline<5:
            line = line.split()
            for it,val in enumerate(line): line[it]=float(val)*multFactor
            lat.append(line)
        if iline==5:
            elements=line.split()
        if iline==6:
            line=line.split()
            for it,val in enumerate(line): line[it]=int(val)
            nelements=line
        if iline==7:
            line=line.split()
            if(line[0]!='Direct'): sys.exit("poscar doesnt have direct coordinates")
            for it,val in enumerate(nelements):
                for dummy in range(val):
                    elementArray.append(elements[it])
        if iline>=8 and iline<8+sum(nelements):
            line=line.split()
            for it,val in enumerate(line): line[it]=float(val)
            coord.append(line)
        iline += 1
    fvasp.close()
    coord=np.array(coord)
    lat=np.array(lat)
    return(lat,coord,elementArray)

def writeqe(qe,lat,coord,elementArray):
    fqe=open(qe,'w')
    for iline in range(5+len(elementArray)):
        if(iline==0):
            outline="CELL_PARAMETERS (angstrom)\n"
        if(iline>=1 and iline<4):
            outline="{:23.16f} {:21.16f} {:21.16f}\n".format(lat[iline-1,0],lat[iline-1,1],lat[iline-1,2])
        if(iline==4):
            outline="ATOMIC_POSITIONS (crystal)\n"
        if(iline>=5):
            outline="{:2s}{:21.16f} {:21.16f} {:21.16f}\n".format(elementArray[iline-5],coord[iline-5,0],coord[iline-5,1],coord[iline-5,2])
        fqe.write(outline)
        iline += 1
    return()

def concatenateFiles(inpfiles,out):
    with open(out, 'w') as outfile:
        for fname in inpfiles:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)

topDir='./'
qeHeader=topDir+'qeHeader'
for vasp in glob.glob(topDir+"*.vasp"):
    qe=vasp[:-4]+'qe'
    qeinp=vasp[:-4]+'inp'
    lat,coord,elementArray = readvasp(vasp)
    writeqe(qe,lat,coord,elementArray)
    concatenateFiles([qeHeader, qe],qeinp)
