#!/bin/bash


# run bmf for different configuration five times for each

for i in {1..5}

do

echo "---------(256,8)----------"
./sssp --input ./input/roadNet-CA.txt --bsize 256 --bcount 8 --output output.txt --method tpe

echo "---------(384,5)----------"
./sssp --input ./input/roadNet-CA.txt --bsize 384 --bcount 5 --output output.txt --method tpe

echo "---------(512,4)----------"
./sssp --input ./input/roadNet-CA.txt --bsize 512 --bcount 4 --output output.txt --method tpe

echo "---------(768,2)----------"
./sssp --input ./input/roadNet-CA.txt --bsize 768 --bcount 2 --output output.txt --method tpe

echo "---------(1024,2)----------"
./sssp --input ./input/roadNet-CA.txt --bsize 1024 --bcount 2 --output output.txt --method tpe

done