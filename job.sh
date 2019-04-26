#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -N sf
#PBS -j oe
#PBS -e error.log
#PBS -o out.log
cd $PBS_O_WORKDIR
for NX in {0..11}
do
    for NY in {0..11}
    do
        qsub ./sf2.out 8 8 12 12 $NX $NY -10.0 10.0 10000 10.0 0.2
    done
done