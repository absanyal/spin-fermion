for NX in {0..31}
do
    for NY in {0..31}
    do
      	sub="amit_${NX}_${NY}"
        echo "#!/bin/bash" > $sub
        echo "#PBS -l nodes=1:ppn=1" >> $sub
        echo "#PBS -N sf" >> $sub
        echo "#PBS -j oe" >> $sub
        echo "#PBS -e error.log" >> $sub
        echo "#PBS -o out.log" >> $sub
        echo 'cd $PBS_O_WORKDIR' >> $sub
        echo "time ./sf2.out 8 8 32 32 $NX $NY -10.0 10.0 10000 10.0 0.2 > outp$
        qsub $sub
    done
done
