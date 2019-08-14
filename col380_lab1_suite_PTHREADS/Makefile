seq: lab1_io.o lab1_sequential.o main_sequential.o
	g++ lab1_io.o lab1_sequential.o main_sequential.o -o seq

p_th: lab1_io.o lab1_pthread.o main_pthread.o
	g++ lab1_io.o lab1_pthread.o main_pthread.o -lpthread -o pth

omp: lab1_io.o lab1_omp.o main_omp.o
	g++ lab1_io.o lab1_omp.o main_omp.o -fopenmp -o omp

check: omp_training.o 
	g++ omp_training.o -fopenmp -o om
	

lab1_io.o: lab1_io.c
	g++ -c lab1_io.c

lab1_sequential.o: lab1_sequential.cpp
	g++ -c lab1_sequential.cpp

main_sequential.o: main_sequential.c
	g++ -c main_sequential.c

lab1_pthread.o: lab1_pthread.cpp
	g++ -c lab1_pthread.cpp -lpthread

main_pthread.o: main_pthread.c
	g++ -c main_pthread.c
	
lab1_omp.o: lab1_omp.cpp
	g++ -c lab1_omp.cpp -fopenmp

main_omp.o: main_omp.c
	g++ -c main_omp.c

omp_training.o: omp_training.cpp
	g++ -c omp_training.cpp -fopenmp
