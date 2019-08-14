#include <stdlib.h>
#include <iostream>
#include <omp.h>
#include <pthread.h>
#include <string.h>
#include <future>

using namespace std;

double pi =0.0;
static long num_steps = 100000;
double step = 1.0/num_steps;
int num_threads = 10;
int thread_part = num_steps/num_threads;


void integrate(int id){
	int i;double private_pi = 0.0;
	if (id==0){num_threads = omp_get_num_threads();}
	int start = id*thread_part;
	for(i =start;i<start + thread_part ;i++ ){
		double x = (i+0.5)*step;
		private_pi = private_pi + (4.0/(1+x*x))*step;
	}
	#pragma omp critical
		{			
			pi = pi + private_pi;
		}
	return;
	
}

// int main (int argc, char *argv[])
// {	
// 	double start,end;
// 	start = omp_get_wtime();	
// 	double A[1000];	
// 	omp_set_num_threads(num_threads);
// 	#pragma omp parallel 
// 	{
// 		int ID = omp_get_thread_num();
// 		integrate(ID);	


















// 		8#pragma omp critical
// 		{
			
			
// 			/*printf("hello(%d)",ID);
// 			printf("world(%d)\n",ID);
// 		//}*/	
// 		//int nt = omp_get_num_threads();
// 		//printf("%d",nt);
// 	}
// 	end = omp_get_wtime();
// 	std::cerr << "Time Taken: "<< end-start << " Value: "<<pi<<"\n";

// 	return 0;

// }




#define NN 10000000
// int * a = (int *) malloc(sizeof(int)*n)
// for (int i = 0; i < 4; ++i)
// {
// 	*(a+i) = i;
//  }

// #define NUM_THREADS 2

int main(int argc, char  *argv[])
{
	// omp_set_num_threads(1);
	// srand(1);
	// int * a = (int *) malloc(sizeof(int)*NN);
	// int * tmp = (int *) malloc(sizeof(int)*NN);
	// int * tmp1 = (int *) malloc(sizeof(int)*NN);


	// int N = atoi(argv[1]);
	// // int NUM_THREADS = atoi(argv[2]);
	// cout << NUM_THREADS;
	// for (int i = 0; i < N; ++i)
	// {
	// 	a[i] = rand()%3*N;
	// }

	// cout << "The  array is : \n";
	// for (int i = 0; i < N; ++i)
	// {
	// 	cout << a[i] <<" ";
	// }
	// cout << "\n";

	// int count;
	// double start,end,startp,endp;

	// start  = omp_get_wtime();
	// for (int i = 0; i < N; ++i)
	// 		{
	// 			int count = 0;
	// 			for (int j = 0; j < N; ++j)
	// 			{
	// 				if ((a[j]<a[i])|| (a[j]==a[i] && j<i))
	// 				{
	// 					#pragma omp critical
	// 					{
	// 						count++;
	// 					}
	// 				}
	// 			}
	// 			tmp1[count] = a[i];
	// 		}
		
	// end = omp_get_wtime();

	 
	// startp = omp_get_wtime();
	#pragma omp parallel 
	{	
			
			// #pragma omp for
			// for (int i = 0; i < N; ++i)
			// {
			// 	int count = 0;
			// 	for (int j = 0; j < N; ++j)
			// 	{
			// 		if ((a[j]<a[i])|| (a[j]==a[i] && j<i))
			// 		{
			// 			#pragma omp critical
			// 			{
			// 				count++;
			// 			}
			// 		}
			// 	}
			// 	tmp[count] = a[i];
			// }
		
		cout << "Hey I am one !!! \n"; 

	}
	// endp = omp_get_wtime();

	// cout << "The sorted array is : \n";
	// for (int i = 0; i < N; ++i)
	// {
	// 	cout << "tmp1 "<< tmp1[i] <<" ";
	// }
	// cout << "\n";

	// cout << "The sorted array is : \n";
	// for (int i = 0; i < N; ++i)
	// {
	// 	cout << "tmp "<< tmp[i] <<" ";
	// }
	// cout << "\n";

	// cout << "Time Taken: \n Sequential : "<< end - start << "\n Parallel : " << endp - startp <<"\n" 
	// << "Speedup : " << (end - start)/(endp - startp) 
	// << "Efficiency : " << (end - start)/(endp - startp) /NUM_THREADS ;  
	return 0;
}


// #include <future>
// #include <iostream>
// #include <string.h>




///////// QUESTION 4 ///////////////////
// void sleep(int n){
// 	for (int i = 0; i < n*10000000; ++i)
// 	{
// 		int a = 100000000;
// 	}
// }

// using namespace std;



// string GoGetEducation(future<string> & s){

//     string child_name = s.get();

//   	cout  << "\nYour child "<< child_name <<"has been admitted. School is going on ...\n";

//   	return "Schooling goes  on and on.." ;
// }


// int main(int argc, char const *argv[])
// {

//     promise<string> child ;
//     future<string> baby = child.get_future();


//     cout << "My Child is Currently : yound  I must send him to School when he is 4 \n" ; 
//     // call 
//     future<string> time_went_past  = async(launch::async,GoGetEducation,ref(baby));
//     cout << "He is 1 Year old ...\n He is 2 Year old ... \n He is 3 Year old ... \n He is 4 Year old I must send him to school now."; 

//     child.set_value("Shashwat");

//    	sleep(3);

//     return 0;
// }
