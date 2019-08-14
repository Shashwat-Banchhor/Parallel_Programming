#include <stdlib.h>
#include <cmath>
#include <iostream>
// #include <omp.h>
#include "lab1_omp.h"


/*
	Arguments:

	N 				   : no. of data points (input)
	K 				   : no. of clusters (input)
	data_points 	   : 1D array of data_points (input)
						 format - 1 data point uses 3 consecutive indices, ie
	     				 -----------------------------------------------
						| pt1_x | pt1_y | pt1_z | pt2_x | pt2_y | pt2_z | ...
		 				 -----------------------------------------------
	data_point_cluster : 1D array of data_points and their cluster id (to be computed)
						 format - 1 data point uses 4 consecutive indices, ie
						  -------------------------------------------------------------------------------
						 | pt1_x | pt1_y | pt1_z | pt1_clusterid | pt2_x | pt2_y | pt2_z | pt2_clusterid | ...
	 					  -------------------------------------------------------------------------------
	 					  cluster_id - range between 0 to K-1
	centroids		   : 1D array of centroids of each iteration (to be computed)
						 format - 1 iteration uses K consecutive indices
						 *include initial centroids in the array (the 1st K entries of the array), ie
						  -------------------------------------------------------------------------------------
						 | cent0_init_x | cent0_init_y | cent0_init_z | ... | centK_init_x centK_init_y centK-1_init_z |
						  -------------------------------------------------------------------------------------
						 | cent0_it1_x | cent0_it1_y | cent0_it1_z | ... | centK_it1_x centK_it1_y centK-1_it1_z |
						  -------------------------------------------------------------------------------------
						 | cent0_it2_x | cent0_it2_y | cent0_it2_z | ... | centK_it2_x centK_it2_y centK-1_it2_z | ...
						  -------------------------------------------------------------------------------------
	num_iterations     : no of iterations of k-means algo = (length of centroids array / K) - 1 (to be computed)

*/


// double start,end;

long changed_points ;
int * new_cluster_sum;
int * new_cluster_count;
float * temp_centroids;
int * data_points_global;
int* data_points_cluster_global ;
float * centroids_global;
int * cluster;

int K_G ;
int N_G;
int t;


// calculates the distance b/w the point and the given cetroid
float distance(int * point , float * centroid){
	float x = *(point) - *(centroid);
	float y = *(point+1) - *(centroid+1);
	float z = *(point+2) - *(centroid+2);
	return sqrt(x*x + y*y + z*z);
}


// write in  the  data_points_cluster
void* logging(){//int K,int N,int * data_points, int ** data_points_cluster,int * cluster){
	
	
	int length_per_thread = N_G;
	int start = 0;
	
	for (int i = start; i < start + N_G; ++i)
	{
		*(data_points_cluster_global+4*i+0) = *(data_points_global + 3*i +0);
		*(data_points_cluster_global+4*i+1) = *(data_points_global + 3*i + 1);
		*(data_points_cluster_global+4*i+2) = *(data_points_global + 3*i + 2);
		*(data_points_cluster_global+4*i+3) = *(cluster + i);
	}
	
	return NULL ;

}

// returns the cluster  in which this pt belongs
int find_cluster_data_point(int K , int  N,float * temporary_centroids,int* point){
	
	float min_distance = 4000.0; //  the distance b/w two points cant be greater than 4000
	int min_cluster = -1;  // initialize cluster with -1
	

	// Check for all centroids 0 1 ... K-1 
	for (int i = 0; i < K; ++i)
	{
		// the distance b/w the datapoint and the centroid
		float d = distance(point,temp_centroids+3*i);

					//std::cerr<<"Printing: "<< d<<" "<< point[0]<<" "<< temp_centroids[3*i]<<std::endl;
		
		// If distanc is less than min _distance then update the min_cluster ;
		
		if ( d< min_distance){
			min_distance = d;
			min_cluster  = i;

		}
	}
	
	// Error Check
	if(min_cluster==-1){std::cerr << "Error in find_cluster_data_point";exit(0);}

	//return the cluster of the datapoint 
	return min_cluster;
}


// updates the cluster and returns the number of points whose assignmen tgot changed
void*  find_cluster_data_points(){//int K , int  N,int * temp_centroids,int ** cluster,int* data_points){
	

	int start = 0;
	

 	

	for(int i =start; i < start+N_G; i++)
	{
		

		int old_cluster = *(cluster + i);
				//std::cerr << "pts-1"<<std::endl;

		//Update the cluster vector for each point
		*(cluster + i) = find_cluster_data_point(K_G,N_G,temp_centroids,data_points_global+3*i);
		
				//std::cerr << "pts0"<<std::endl;
		int new_cluster = *(cluster + i);
		
		new_cluster_count[new_cluster] += 1;
		new_cluster_sum[3*new_cluster] += data_points_global[3*i] ;
		new_cluster_sum[3*new_cluster+1] += data_points_global[3*i+1] ;
		new_cluster_sum[3*new_cluster+2] += data_points_global[3*i+2] ;
		
		
		if (new_cluster!=old_cluster){
			changed_points += 1;
		}

	}

  

  	return NULL;
	

	
}

void kmeans_sequential(int N,int K,int* data_points,int** data_point_cluster,float** centroids,int* num_iterations){

	// initialization of output #####################################


	// start = omp_get_wtime();
	K_G = K;
	N_G = N;
	temp_centroids = (float*)malloc(3*K*sizeof(float));
	data_points_global = data_points;
	data_points_cluster_global = (int *) malloc(N*4*sizeof(int));
	centroids_global = (float *) malloc((N*3+100)*sizeof(float));
	cluster = (int *) malloc(N*sizeof(int)); // will store the cluster in which ith point is present curretnly
	*num_iterations = 0;
	
	//will store the sum of the co-ordinates of the points in a cluster
	new_cluster_sum = (int *)calloc(3*K,sizeof(int)); 
	//will store the number of points in a cluster
	new_cluster_count = (int *)calloc(K,sizeof(int));

	changed_points = 0;
	int convergence = 0;//0.0001*N; // when should the algo converge

	

	// initialization ends #######################################



	
	int *tid = (int *) malloc (sizeof (int) * t);



	///***** initial logging ******//  ?????????????(check not same if time permits)
		//centroids at beginning
	int centroid_index = 0; // index_to_write in int** centroids
    for(int j=0;j<K;j++){	
        int random_selection = rand()%N;
    	for (int i = 0; i <3; ++i)
    	{
    		
    		temp_centroids[3*j + i] = data_points[random_selection*3 + i]; // initialize cetroids with first k datapoints
    		centroids_global[3*j + i] = data_points[random_selection*3 + i];
    					
    	}
    }
    //******************************//

	
	

	while(1){ 	
	  	
		
	
		find_cluster_data_points();
		

		
		 
		

	 	// one iteration is done
	 	*num_iterations += 1;



	 	// write the points on centroids_global  and update temp_centroids
	 	// re-initialize the new_cluster_count & new_cluster_sum for next iteration
	 	for (int index = 0; index < K ; ++index)
	  	{

	  		for (int index_cord = 0; index_cord < 3; ++index_cord)
	  		{
	  		
	  			if (new_cluster_count[index]!=0){
	  			    
	  			    
	  			
    				 	
	  			    temp_centroids[3*index+index_cord] = float(new_cluster_sum[3*index+index_cord]) /new_cluster_count[index];
	  			
	  			    
    	  			centroids_global[*num_iterations*3*K + 3*index+index_cord] = float(new_cluster_sum[3*index+index_cord]) /new_cluster_count[index];
    	  			
    	  			
    	  			new_cluster_sum[3*index+index_cord] =  0;
    				
    	  		}
	  			else{
				
    				temp_centroids[3*index+index_cord] =centroids_global[*num_iterations*3*K + 3*index+index_cord - 3*K] ;
    	  			centroids_global[*num_iterations*3*K + 3*index+index_cord] = centroids_global[*num_iterations*3*K + 3*index+index_cord - 3*K] ;
    	  			new_cluster_sum[3*index+index_cord] =  0;
	  			    
	  			}
	  		}
	  		new_cluster_count[index] = 0;
		
	  	}


		if(changed_points <= convergence){
			break;
		}
		changed_points = 0;

	}

	
	//  update the output data_points_cluster_global
		
	logging();
		
	

	/// Assigning the outputs
	*centroids = centroids_global;
	*data_point_cluster = data_points_cluster_global;

	// end = omp_get_wtime();
	// std::cerr << end - start <<"\n";

	return ;

}


