CC=g++
CFLAGS=-Wall

main: main.cpp main_naive main_warp main_shared
	$(CC) -o main main.cpp

main_naive: main_naive.cu
	nvcc -o main_naive main_naive.cu

main_warp: main_warp.cu
	nvcc -o main_warp main_warp.cu

main_shared: main_shared.cu
	nvcc -o main_shared main_shared.cu

clean:
	rm -f main main_naive main_warp main.o 



