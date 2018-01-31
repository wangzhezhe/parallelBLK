
Name:ZheWang
netid:zw241

Please use following command to run source code:

```
make clean && make

# [prject1]original spmv for project1
./spmv -mat <matrix path> -ivec <vector path> -alg segment -blksize  512 -blknum 4

# [prject1]optimization the thread partition based on original implementation
./spmv -mat <matrix path> -ivec <vector path> -alg design -blksize  512 -blknum 4

# [prject3A] optimization for project3 based on fist touch packaging algorithm
./spmv -mat <matrix path> -ivec <vector path> -alg optfirst -blksize  512 -blknum 4

# [prject3A] optimization for project3 based on graph partition algorithm
# [test ok for sample matrix for ./matrix/test8.mtx ./matrix/vector8.txt, only work well when blocknumber*blocksize > number of non-zero element]

./spmv -mat ./matrix/test8.mtx -ivec ./matrix/vector8.txt -alg optgraph -blksize 4 -blknum 2

#long running time for large matrix for searching the largest edge
./spmv -mat <matrix path> -ivec <vector path> -alg optgraph -blksize 512 -blknum 4

# /bin/bash evaluationb.sh could test all the matrix in requirements of project3A by using -alg optfirst
```