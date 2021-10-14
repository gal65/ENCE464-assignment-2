#!/bin/bash

ns=(101 201 301 401 501 601 701 801 901)
t=4
i=100
exe=main

for n in ${ns[@]}; do
    echo "running for n=$n"
    echo "$exe -n $n -t $i -t $t"
    time ( ./$exe -n $n -t $i -t $t > /dev/null 2>&1 ) 2>&1
done