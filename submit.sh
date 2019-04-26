for NX in {0..11}
do
    for NY in {0..11}
    do
        echo $NX $NY
        time ./sf2.out 8 8 12 12 $NX $NY -10.0 10.0 10000 10.0 0.2
    done
done