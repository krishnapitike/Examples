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

max=
system=
mkdir input
cp    *     input/

for i in `seq 1 9` ; do if [ -a lattice.$i ] ; then m=`echo $((i+1))` ; fi ; done

for n in `seq $m $max` ; do

  if [ $m -gt 1 ] ; then
    grep -A69 "CELL_PARAMETERS (angstrom)" ${system}.64atom.pseudocubic.bulk.pw.vc-relax.`echo $((n-1))`.out > temp.allcelldata
    tail -70  temp.allcelldata                                                                               > temp.celldata
    head -4   temp.celldata                                                                                  > lattice.${n}
    tail -65  temp.celldata                                                                                  > coord.${n}
    rm        temp.*
  fi

  if [ `wc -l lattice.$n` -eq 0 ] ; then cp lattice.`echo $((n-1))` lattice.${n} ; fi
  if [ `wc -l coord.$n` -eq 0 ]   ; then cp coord.`echo $((n-1))`   coord.${n}   ; fi
  
  cat       incar potcar lattice.${n} kpoints coord.${n}                                        > ${system}.64atom.pseudocubic.bulk.pw.vc-relax.${n}.inp
  mpirun -n 72 pw.x -nk 2 -nb 6 -input ${system}.64atom.pseudocubic.bulk.pw.vc-relax.${n}.inp   > ${system}.64atom.pseudocubic.bulk.pw.vc-relax.${n}.out

  grep -A69 "CELL_PARAMETERS (angstrom)" ${system}.64atom.pseudocubic.bulk.pw.vc-relax.${n}.out > temp.allcelldata
  tail -70  temp.allcelldata                                                                    > temp.celldata
  head -4   temp.celldata                                                                       > lattice.`echo $((n+1))`
  tail -65  temp.celldata                                                                       > coord.`echo $((n+1))`
  rm        temp.*

done
rm -r scratch
