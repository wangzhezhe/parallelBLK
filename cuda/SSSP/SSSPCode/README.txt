compile:
$ make

run code:

run sssp by balman-fold algorithm and sort edge list by dst:
$ ./sssp --input <input file path> --bsize 512 --bcount 2 --output output.txt --method bmf


run sssp by balman-fold algorithm + to-process optimization and sort edge list by src:
$ ./sssp --input <input file path> --bsize 512 --bcount 2 --output output.txt --method tpe

run ssp by balman fold algorithm and sort edge list by src:
$ ./sssp --input <input file path> --bsize 512 --bcount 2 --output output.txt --method opt