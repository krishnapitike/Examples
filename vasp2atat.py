import ase.io
import numpy as np
topDir='./'
strin='POSCAR.vasp'
a=[[3.83,0.00,0.00],    #a lat constant
   [0.00,3.83,0.00],    #b
   [0.00,0.00,12.6],    #c
   [20,   0,   0   ],   #supercell size in a direction
   [0,   20,   0   ],   #supercell size in b direction
   [0,    0,   10  ]]   #supercell size in c direction

system=ase.io.read(topDir+strin,format='vasp')
coords=system.get_scaled_positions()
specList=system.get_chemical_symbols()
for it,valit in enumerate(coords):
    for jt,valjt in enumerate(valit):
        coords[it,jt] = valjt*a[jt+3][jt]
print(len(coords))
print(coords[0:10])
fstr=open(topDir+'str.out','w')
for it,val in enumerate(a):
    line=" {:6.3f} {:6.3f} {:6.3f}\n".format(val[0],val[1],val[2])
    fstr.write(line)
for it,val in enumerate(coords):
    line=" {:6.3f} {:6.3f} {:6.3f}  {:s}\n".format(float(val[0]),float(val[1]),float(val[2]),specList[it])
    fstr.write(line)
fstr.close()
