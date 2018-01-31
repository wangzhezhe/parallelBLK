#!/bin/bash


# run bmf for different configuration five times for each

for i in {1..20}

do

echo "---------(256,$i)---------"
./sssp --input ./input/roadNet-CA.txt --bsize 256 --bcount $i --output output.txt --method bmf

done