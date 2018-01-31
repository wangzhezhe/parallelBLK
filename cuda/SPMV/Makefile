#For this file, I owe my soul to:
#https://devblogs.nvidia.com/parallelforall/separate-compilation-linking-cuda-device-code

CC = /usr/local/cuda/bin/nvcc
GCC = g++
objects = mmio.o main.o genresult.o

all: $(objects) 
	$(CC) -arch=sm_30 -g -O0 $(objects) -o spmv

mmio.o: mmio.c
	$(GCC) -w mmio.c -c 

main.o: main.c
	$(CC) -x cu -arch=sm_30 -I. -dc main.c -o $@

genresult.o: genresult.c spmv_segment.cu spmv_design_opt.cu
	$(CC) -x cu -arch=sm_30 -g -I. -dc genresult.c -o $@

%.o: %.c
	$(CC) -x cu -arch=sm_30 -I. -dc $< -o $@

clean:
	rm -f *.o spmv

#Quick conveniences for submission etc.
tar:
	tar -cvf spmv_proj.tar *.c *.cuh *.pdf makefile

untar:
	tar -xvf spmv_proj.tar
