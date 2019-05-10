for J in 0.5 1.0 10.0
do
    for fill in 0.5 0.8 0.1
    do
        mkdir thetascan_${J}_${fill}
        cp sf_thetascan.out thetascan_${J}_${fill}/.
        cd thetascan_${J}_${fill}
        sub="amit_${J}_${fill}"
        echo "#!/bin/bash" > $sub
        echo "#PBS -l nodes=1:ppn=1" >> $sub
        echo "#PBS -N sf" >> $sub
        echo "#PBS -j oe" >> $sub
        echo "#PBS -e error.log" >> $sub
        echo "#PBS -o out.log" >> $sub
        echo 'cd $PBS_O_WORKDIR' >> $sub
        echo "time ./sf_thetascan.out 16 16 ${J} 100 ${fill}" >> outp$
        qsub $sub
        cd ../
    done
done