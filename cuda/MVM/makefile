CC=g++
CFLAGS=-Wall

main: main.cpp main_naive main_warp
	$(CC) -o main main.cpp

main_naive: main_naive.cu
	nvcc -o main_naive main_naive.cu

main_warp: main_warp.cu
	nvcc -o main_warp main_warp.cu

clean:
	rm -f main main_naive main_warp main.o 



