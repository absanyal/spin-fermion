./compile.sh
for NX in {0..3}
do
    for NY in {0..3}
    do
        echo $NX $NY
        time ./sf2.out 4 4 4 4 $NX $NY -10.0 10.0 10000 10.0 0.2
    done
done