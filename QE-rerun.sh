#!/bin/bash

#PBS -W group_list=cades-theory  ##group_list=cades-birthright
#PBS -A theory                   ##birthright
#PBS -q skylake
#PBS -N test
#PBS -M pitikek@ornl.gov
#PBS -l nodes=2:ppn=36
#PBS -l walltime=48:00:00
#PBS -l qos=std                ##burst, std, long
#PBS -V

export OMP_NUM_THREADS=1
module load PE-intel/1.0
module load QE/6.1

cd  $PBS_O_WORKDIR
pwd

natm=      #number of atoms
max=       #maximum number of iterations
system=    #name of the system
phase=     #cubic, tetragonal etc
strain=    #bulk, st_-1, etc
calc=      #pw.vc-relax, pw.scf etc
fname=     #dont worry if not specified
fname=`echo ${system}.${natm}atom.${phase}.${strain}.${calc}`

mkdir input
cp    *     input/

for i in `seq 1 9` ; do if [ -a lattice.$i ] ; then m=`echo $((i+1))` ; fi ; done

for n in `seq $m $max` ; do

  if [ $m -gt 1 ] ; then
    grep -A`echo $((natm+5))` "CELL_PARAMETERS (angstrom)" ${fname}.`echo $((n-1))`.out > temp.allcelldata
    tail -`echo $((natm+6))`  temp.allcelldata                                                                                  > temp.celldata
    head -4                   temp.celldata                                                                                                     > lattice.${n}
    tail -`echo $((natm+1))`  temp.celldata                                                                                     > coord.${n}
    rm                        temp.*
  fi

  if [ -s lattice.$n ] ; then echo "lattice.${n} is good" ; else cp lattice.`echo $((n-1))` lattice.${n} ; fi
  if [ -s coord.$n ]   ; then echo "coord.${n}   is good" ; else cp coord.`echo $((n-1))`   coord.${n}   ; fi
  
  cat       incar potcar lattice.${n} kpoints coord.${n} > ${fname}.${n}.inp
  mpirun -n 72 pw.x -nk 2 -nb 6 -input ${fname}.${n}.inp > ${fname}.${n}.out

  grep -A`echo $((natm+5))` "CELL_PARAMETERS (angstrom)" ${fname}.${n}.out > temp.allcelldata
  tail -`echo $((natm+6))`  temp.allcelldata                               > temp.celldata
  head -4                   temp.celldata                                  > lattice.`echo $((n+1))`
  tail -`echo $((natm+1))`  temp.celldata                                  > coord.`echo $((n+1))`
  rm                        temp.*

done
rm -r scratch
