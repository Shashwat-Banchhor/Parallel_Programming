omp: lab2.o
	g++ lab2.o -fopenmp -o omp

lab2.o: lab2.cpp
	g++ -c lab2.cpp -fopenmp
