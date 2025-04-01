from ase.io.lammpsdata import read_lammps_data, write_lammps_data
from ase.build import make_supercell

#################################################################
# reading and writing lammps data
#################################################################
# reading lammps data
unit_cell=read_lammps_data('init.data',style='atomic')

# making supercell
transformation = [[2,0,0],[0,2,0],[0,0,2]]
super_cell=make_supercell(unit_cell, transformation)

# writing lammps data
write_lammps_data('supercell.data',super_cell,atom_style='atomic')

