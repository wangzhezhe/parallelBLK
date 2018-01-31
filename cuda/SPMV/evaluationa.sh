#/bin/bash


#MLIST="can con mac pdb web shi wat cir ful mc2"

#MLIST="il2 kne nc2 pds roac roau webg webs hu1 hu2"

#for mname in ${MLIST}

#do

mname="hu1"
echo "matrix name ${mname}"
mpath=/.freespace/zw241/opt/matrix/${mname}.mtx
vpath=/.freespace/zw241/opt/vector/${mname}.txt

#mpath=/.freespace/zw241/matrix/${mname}.mtx
#vpath=/.freespace/zw241/vector/${mname}.txt

echo ${mpath}
echo ${vpath}

#echo "------(10,192)------"

#./spmv -mat ${mpath} -ivec ${vpath} -alg segment -blksize 192 -blknum 10

#echo "------(8,256)------"

#./spmv -mat ${mpath} -ivec ${vpath} -alg segment -blksize 256 -blknum 8


#echo "------(5,384)------"

#./spmv -mat ${mpath} -ivec ${vpath} -alg segment -blksize 384 -blknum 5


#echo "------(4,512)------"

#./spmv -mat ${mpath} -ivec ${vpath} -alg segment  -blksize 512 -blknum 4


#echo "------(2,768)------"

#./spmv -mat ${mpath} -ivec ${vpath} -alg segment -blksize 768 -blknum 2


#echo "------(2,1024)------"

#./spmv -mat ${mpath} -ivec ${vpath} -alg segment -blksize 1024 -blknum 2

#done


#for mname in ${MLIST}

#do


#mpath=/.freespace/zw241/opt/matrix/${mname}.mtx
#vpath=/.freespace/zw241/opt/vector/${mname}.txt

#echo ${mpath}
#echo ${vpath}

echo "------test opt design------"

echo "------(10,192)------"

./spmv -mat ${mpath} -ivec ${vpath} -alg design -blksize 192 -blknum 10

echo "------(8,256)------"

./spmv -mat ${mpath} -ivec ${vpath} -alg design -blksize 256 -blknum 8


echo "------(5,384)------"

./spmv -mat ${mpath} -ivec ${vpath} -alg design -blksize 384 -blknum 5


echo "------(4,512)------"

./spmv -mat ${mpath} -ivec ${vpath} -alg design  -blksize 512 -blknum 4


echo "------(2,768)------"

./spmv -mat ${mpath} -ivec ${vpath} -alg design -blksize 768 -blknum 2


echo "------(2,1024)------"

./spmv -mat ${mpath} -ivec ${vpath} -alg design -blksize 1024 -blknum 2




#mpath=/.freespace/zw241/opt/matrix/${mname}.mtx
#vpath=/.freespace/zw241/opt/vector/${mname}.txt

#echo ${mpath}
#echo ${vpath}

sleep 1

echo "------test opt first------"

echo "------(10,192)------"

./spmv -mat ${mpath} -ivec ${vpath} -alg optfirst -blksize 192 -blknum 10

echo "------(8,256)------"

./spmv -mat ${mpath} -ivec ${vpath} -alg optfirst -blksize 256 -blknum 8


echo "------(5,384)------"

./spmv -mat ${mpath} -ivec ${vpath} -alg optfirst -blksize 384 -blknum 5


echo "------(4,512)------"

./spmv -mat ${mpath} -ivec ${vpath} -alg optfirst  -blksize 512 -blknum 4


echo "------(2,768)------"

./spmv -mat ${mpath} -ivec ${vpath} -alg optfirst -blksize 768 -blknum 2


echo "------(2,1024)------"

./spmv -mat ${mpath} -ivec ${vpath} -alg optfirst -blksize 1024 -blknum 2

#done