#/bin/bash

#MLIST="il2 kne nc2 pds roac roau webg webs hu1 hu2"
MLIST="roac"
for mname in ${MLIST}

do
#mname="pwt"
echo "matrix name ${mname}"
mpath=/.freespace/zw241/opt/matrix/${mname}.mtx
vpath=/.freespace/zw241/opt/vector/${mname}.txt

echo ${mpath}
echo ${vpath}

#echo "------(1,256)------"

#./spmv -mat ${mpath} -ivec ${vpath} -alg optfirst -blksize 256 -blknum 1

#echo "------(2,256)------"

#./spmv -mat ${mpath} -ivec ${vpath} -alg optfirst -blksize 256 -blknum 2


echo "------(4,512) optfirst ------"

./spmv -mat ${mpath} -ivec ${vpath} -alg optfirst -blksize 512 -blknum 4

sleep 1

echo "------(4,512) design ------"

./spmv -mat ${mpath} -ivec ${vpath} -alg design -blksize 512 -blknum 4

#echo "------(6,256)------"

#./spmv -mat ${mpath} -ivec ${vpath} -alg optfirst -blksize 256 -blknum 6


#echo "------(8,256)------"

#./spmv -mat ${mpath} -ivec ${vpath} -alg optfirst -blksize 256 -blknum 8


#echo "------(10,256)------"

#./spmv -mat ${mpath} -ivec ${vpath} -alg optfirst -blksize 256 -blknum 10


#echo "------(12,256)------"

#./spmv -mat ${mpath} -ivec ${vpath} -alg optfirst -blksize 256 -blknum 12

#echo "------(14,256)------"

#./spmv -mat ${mpath} -ivec ${vpath} -alg optfirst -blksize 256 -blknum 14

#echo "------(16,256)------"

#./spmv -mat ${mpath} -ivec ${vpath} -alg optfirst -blksize 256 -blknum 16

done