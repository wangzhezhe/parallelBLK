#!/bin/bash

# test for src and dst

# sort by dst
#./sssp --input ./input/amazon0312.txt --bsize 1024 --bcount 2 --output output.txt --method bmf
#./sssp --input ./input/roadNet-CA.txt --bsize 1024 --bcount 2 --output output.txt --method bmf
#./sssp --input ./input/soc-LiveJournal1.txt --bsize 1024 --bcount 2 --output output.txt --method bmf
#./sssp --input ./input/soc-pokec-relationships.txt --bsize 1024 --bcount 2 --output output.txt --method bmf
#./sssp --input ./input/web-Google.txt --bsize 1024 --bcount 2 --output output.txt --method bmf

# sort by src

#./sssp --input ./input/amazon0312.txt --bsize 1024 --bcount 2 --output output.txt --method opt
#./sssp --input ./input/roadNet-CA.txt --bsize 1024 --bcount 2 --output output.txt --method opt
#./sssp --input ./input/soc-LiveJournal1.txt --bsize 1024 --bcount 2 --output output.txt --method opt
#./sssp --input ./input/soc-pokec-relationships.txt --bsize 1024 --bcount 2 --output output.txt --method opt
#./sssp --input ./input/web-Google.txt --bsize 1024 --bcount 2 --output output.txt --method opt

#tpe sort by dst
#./sssptpedst --input ./input/amazon0312.txt --bsize 1024 --bcount 2 --output output.txt --method tpe
#./sssptpedst --input ./input/roadNet-CA.txt --bsize 1024 --bcount 2 --output output.txt --method tpe
./sssptpedst --input ./input/soc-LiveJournal1.txt --bsize 1024 --bcount 2 --output output.txt --method tpe
#./sssptpedst --input ./input/soc-pokec-relationships.txt --bsize 1024 --bcount 2 --output output.txt --method tpe
#./sssptpedst --input ./input/web-Google.txt --bsize 1024 --bcount 2 --output output.txt --method tpe



#tpe sort by src
#./sssp --input ./input/amazon0312.txt --bsize 1024 --bcount 2 --output output.txt --method tpe
#./sssp --input ./input/roadNet-CA.txt --bsize 1024 --bcount 2 --output output.txt --method tpe
./sssp --input ./input/soc-LiveJournal1.txt --bsize 1024 --bcount 2 --output output.txt --method tpe
#./sssp --input ./input/soc-pokec-relationships.txt --bsize 1024 --bcount 2 --output output.txt --method tpe
#./sssp --input ./input/web-Google.txt --bsize 1024 --bcount 2 --output output.txt --method tpe