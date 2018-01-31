CC = /usr/local/cuda/bin/nvcc

segscan: *.cu *.c
	$(CC) -std=c++11 segscan_naive.cu utils.c -O3 -arch=sm_30 -o segscan_naive
#parse_graph.cpp impl1.cu impl2.cu opt.cu 

clean:
	rm -f *.o segscan
