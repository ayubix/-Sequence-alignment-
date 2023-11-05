build:
	mpicxx -c finalParllel.c -o finalParllel.o -fopenmp
	nvcc -c cuda.cu -o cuda.o 
	mpicxx -o exec finalParllel.o cuda.o /usr/local/cuda-11.0/targets/x86_64-linux/lib/libcudart_static.a -ldl -lrt -fopenmp

clean:
	rm -f *.o exec

run:
	mpiexec -n 2 ./exec
