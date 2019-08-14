#include "lab1_sequential.h"
#include <iostream>
#include <cmath>
//---------------------------------------------------------------------
	// int N;					//no. of data points (input)
	// int K;					//no. of clusters to be formed (input)
	// int* data_points;		//data points (input)
	// int* cluster_points;	//clustered data points (to be computed)
	// int* centroids;			//centroids of each iteration (to be computed)
	// int num_iterations;    //no of iterations performed by algo (to be computed)
	//---------------------------------------------------------------------



// calculates the distance b/w the point and the given cetroid
float distance(int * point , float * centroid){
	float x = *(point) - *(centroid);
	float y = *(point+1) - *(centroid+1);
	float z = *(point+2) - *(centroid+2);
	return sqrt(x*x + y*y + z*z);
}


// updates the new centroid/means in temp_centroids
void update_temp_centroids(int K,int N,int * data_points, float ** temp_centroids,int * cluster){
	
	int  * cluster_count = (int *) calloc(K,sizeof(int));
	float * new_means = (float *) calloc(3*K,sizeof(float));
	

	// for (int i = 0; i < N; ++i)
	// {
	//     if (*(cluster_count+i/3)!= 0 ){
	// 	    *(*temp_centroids + i) = *(cluster_sum + i)/ *(cluster_count+i/3);
	//     }
	// }

	for (int i = 0; i < N; ++i)
	{
		
		int c = cluster[i];
		cluster_count[c] += 1;
		new_means[3*c] = (new_means[3*c]*(cluster_count[c]-1)+data_points[3*i])/cluster_count[c]; 
		new_means[3*c+1] = (new_means[3*c+1]*(cluster_count[c]-1)+data_points[3*i+1])/cluster_count[c];
		new_means[3*c+2] = (new_means[3*c+2]*(cluster_count[c]-1)+data_points[3*i+2])/cluster_count[c];
		/*for(int j=0;j<3*K;j++){
		    std::cerr << new_means[j]<<" ";
		}
		std::cerr<<std::endl;
		for(int j=0;j<N;j++){
		    std::cerr << cluster[j]<<" ";
		}
		std::cerr<<std::endl;*/
	}


	for (int i = 0; i < 3*K; ++i)
	{
		temp_centroids[0][i] = (float) new_means[i];
	}

	return ;

}

// returns the cluster  in which this pt belongs
int find_cluster_data_point(int K , int  N,float * temp_centroids,int* point){
	
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
int  find_cluster_data_points(int K , int  N,float * temp_centroids,int ** cluster,int* data_points){
	int changed_points = 0;

	for (int i = 0; i < N; ++i)
	{
		int old_cluster = *(*cluster + i);
				//std::cerr << "pts-1"<<std::endl;

		//Update the cluster vector for each point
		*(*cluster + i) = find_cluster_data_point(K,N,temp_centroids,data_points+3*i);
		
				//std::cerr << "pts0"<<std::endl;
		int new_cluster = *(*cluster + i);
				//std::cerr << i<<" "<<new_cluster<<std::endl;
		
		// Check whether its assignment changed or not.
		if (new_cluster!=old_cluster){
			changed_points += 1;
		}

				// 		else{
				// //std::cerr << "pts2"<<std::endl;
				// 			cluster_count[0][new_cluster] += 1;
				// 			cluster_count[0][old_cluster] -= 1;
				// //std::cerr << "pts3"<<std::endl;
				// 			*(*cluster_sum + 3*new_cluster+ 1) += *(data_points + 3*i +1);
				// 			*(*cluster_sum + 3*new_cluster+ 2) += *(data_points + 3*i +2);
				// 			*(*cluster_sum + 3*new_cluster+ 3) += *(data_points + 3*i +3);

				// 			*(*cluster_sum + 3*old_cluster+ 1) -= *(data_points + 3*i +1);
				// 			*(*cluster_sum + 3*old_cluster+ 2) -= *(data_points + 3*i +2);
				// 			*(*cluster_sum + 3*old_cluster+ 3) -= *(data_points + 3*i +3);
				// 		}
	}

	return changed_points;
}

void kmeans_sequential(int N,int K,int* data_points,int** data_point_cluster,float** centroids,int* num_iterations){
	
	
	// initialization #######################################

	*data_point_cluster = (int *) malloc(N*4*sizeof(int));
	*centroids = (float *) malloc((N*3+10)*sizeof(float));
	float * temp_centroids = (float *) malloc(K*3*sizeof(float));
	int * cluster = (int *) malloc(N*sizeof(int)); // will store the cluster in which ith point is present curretnly

	srand(1); // initialize the random seed

	*num_iterations = 0;

	// initialization ends #######################################
	
	int convergence = 0;//0.0001*N; // when should the algo converge
			//std::cerr << "S-1"<<std::endl;
	
	//centroids at beginning
	int centroid_index = 0; // index_to_write in int** centroids
    for(int j=0;j<K;j++){	
        int random_selection = rand()%N/3;
    	for (int i = 0; i <3; ++i)
    	{
    		
    		temp_centroids[3*j + i] = data_points[random_selection*3 + i]; // initialize cetroids with first k datapoints
    		centroids[0][3*j + i] = data_points[random_selection*3 + i];
    					
    	}
    }
    	
	centroid_index += 3*K;

	while(true){
		

				//std::cerr << "S01"<<std::endl;

		// number of points whose cluster assignment changed
		int pts_changed = find_cluster_data_points(K,N,temp_centroids,&cluster,data_points);
				//std::cerr << "S02"<<std::endl;		

		// new means of the K clusters
		update_temp_centroids(K,N,data_points,&temp_centroids,cluster);
		
		// one iteration completed
		*num_iterations = *num_iterations + 1;
		
		
		//*** log the centroids ***//
		for (int i = centroid_index; i < centroid_index + 3*K; ++i)
		{
			
			*(*centroids+i) = *(temp_centroids+i);
			
		}
		centroid_index += 3*K;
	    //*** done logging the centroids ***//
		
		
				//std::cerr << centroid_index<<std::endl ;;
		
		
		//check convergence
			//if convergence then exit;
		if ( pts_changed<= convergence){
		
				//std::cerr << pts_changed <<std::endl;			
			break;
		}

	}


			//std::cerr << "S1"<<std::endl;



	
	// make the data_point_cluster array  
	int k = 0;
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			*(*data_point_cluster+k) = *(data_points+i*3+j);
			k = k+1;
		}

		//which cluster is belongs
		*(*data_point_cluster+ k) = *(cluster + k/4);
		k=k+1;
	}

		

	return ;



}
