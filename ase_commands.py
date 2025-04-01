from ase.io.lammpsdata import read_lammps_data, write_lammps_data

#################################################################
# reading and writing lammps data
#################################################################
# reading lammps data
sys=read_lammps_data('init.data',style='atomic')

# writing lammps data
write_lammps_data('test.data',sys,atom_style='atomic')

