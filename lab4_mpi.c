#include "lab4_mpi.h"

#include <malloc.h>
// #include "mpi.h"


// Utility function
int
min (int a, int b){
	return a<b?a:b;
}

int 
Ones(int * S, int n, int index, int sz){
	
	for (int i = 0; i < sz; ++i)
	{	
		if(i+index>n){
			return 0;
		}
		if (S[index+i]!=1){
			return 0;
		}
			
	}
	return 1;
		
}

int 
Equal(char * X , char *Y, int sz ){
	for (int i = 0; i < sz; ++i)
	{
		if(X[i]!= Y[i]){
			return 0;
		}
	}
	return 1;
}

int
getPeriod(char *Y, int m){
	int curr_period = 0;
	int shift  = 1;
	for (int shift = 1; shift < m  ; ++shift)
	{
		if(Equal(Y,Y+shift,m-shift)){
			return shift;
		}

	}

	return m-1;
}



int *
WITNESS(char *Y, int m, int size){
	int  * WitnessArr = (int *)malloc(sizeof(int)*size);
	WitnessArr[0] = 0;
	for (int i = 1; i < size; ++i)
	{
		int k = 0;
		for (k = 0; k < m; ++k)
		{
			if (Y[k]!=Y[i+k]){
				break;
			}
		}
		WitnessArr[i] = k;
	}
	return WitnessArr;
}


int 
duel(char * Z, int n, char * Y, int m, int *WtnsArr, int i, int j)
{
	int k = WtnsArr[j-i];
	if ((j+k) >= n){
		return i;
	}
	else{
		if (Z[j+k]!=Y[k]){
			return i;
		}
		return j;
	}

}

int
Occurs(char * T, int index, char *P, int m){
	// printf("indx: %d\n",index );
	int j =0;
	for (int i = 0; i < m; ++i , index++)
	{
		// printf("T: %c P: %c i: %d \n",T[index],P[i],i );
		if (T[index]!=P[i]){
			return 0;
		}
		
	}
	return 1;
} 


char *
makeu2v(char * u, int m_u, char *v, int m_v){

	int size = (2*m_u+m_v);
	char *u2v = (char *)malloc(sizeof(char)*size);

	for (int i = 0; i < size; ++i)
	{
		if (i<m_u){
			u2v[i] = u[i];
			// printf("a: %d %c\n",i,u2v[i] );
		}
		else if(i<2*m_u){
			u2v[i] = u[i-m_u];
			// printf("b: %d %c\n",i,u2v[i] );

		}
		else{
			// printf("vim %d\n", i-2*m_u);
			u2v[i] = v[i-2*m_u];
			// printf("c: %d %c\n",i,u2v[i] );

		}
	}
	return u2v;

}

// P is non-periodic
int *
NP_periodic_pattern_matching(char *T, int n, char *P, int m,int *WtnsArr, int * total_match)
{
	


	// Algorithm
	int blk_sz = ((m+1)/2);
	int num_blocks = n/blk_sz;
	int *T_partition_indx = (int *) malloc(sizeof(int)*num_blocks);
	int *T_partition_sz   = (int *) malloc(sizeof(int)*num_blocks);
	int blk_start_indx = 0;
	int blk_id = 0;
	while(blk_start_indx < n){
		T_partition_indx[blk_id] = blk_start_indx;	
		blk_start_indx = min(blk_start_indx + blk_sz , n );
		T_partition_sz[blk_id] = min(blk_sz, n- blk_sz);
		blk_id ++;

	}

	
	// T_partition_indx and T_partition_sz check
	// for (int i = 0; i < num_blocks; ++i)
	// {
	// 	printf("%d %d\n",T_partition_indx[i] , T_partition_sz[i] );
	// }

	// Initialize
	//--------------------------
	// to_compute
	int * match_positions; 
	// intermidiate storage
	int * potential_positions;
	int * matched_positions;
	int match_count = 0;
	potential_positions = (int *)malloc(sizeof(int)*num_blocks);
	matched_positions = (int *)malloc(sizeof(int)*num_blocks);

	//---------------------------


	for (int b = 0; b< num_blocks;b++){
		int i = 0;
		for (int blks = 0; blks < b; blks++)
		{
			i = i + T_partition_sz[blks];
		}
		int b_loop_sz = 0;
		for (int blks = 0; blks < b+1; blks++)
		{
			b_loop_sz = b_loop_sz + T_partition_sz[blks];
		}


		for(int j = i+1 ; j<b_loop_sz;j++){
			i = duel(T, n, P, m, WtnsArr, i, j);
		}
		potential_positions[b] = i;

	}

	for (int i = 0; i < num_blocks; ++i)
	{
		if(Occurs(T,potential_positions[i],P,m)){
			matched_positions[i] = 1;
			match_count += 1;
		}
		else{
			matched_positions[i] = 0;
		}

	}

	match_positions = (int *) malloc(sizeof(int)*match_count);
	int m_indx = 0;
	for (int i = 0; i < num_blocks ; ++i)
	{	
		if (matched_positions[i]==1){
			match_positions[m_indx] = potential_positions[i];
			m_indx +=1;
		}
	}

	// for (int i = 0; i < num_blocks ; ++i)
	// {
	// 	printf("%d ",potential_positions[i] );
	// }

	// printf("%d\n",match_count );
	// for (int i = 0; i < match_count ; ++i)
	// {
	// 	printf("%d\n",match_positions[i] );
	// }

	total_match[0] = match_count;
	return match_positions;

}


int 
P_periodic_pattern_matching(char *T, int n, char *P, int m, int period, int ** matching_positions, int *match_counts)
{
	int p = period;
	int m_prime = (2*p-1);
	char * P_prime = (char *)malloc(sizeof(char)*m_prime);
	for (int i = 0; i < m_prime; ++i)
	{
		P_prime[i] = P[i];
	}

	for (int i = 0; i < m_prime  ; ++i)
	{	
		printf("%c ",P_prime[i] );
	}

	int p_prime = getPeriod(P_prime,m_prime);
	printf("P_prime: %d\n",p_prime );
	int wa_sz = min((m_prime+1)/2,p_prime);
	int * wa =  WITNESS(P_prime, m_prime, wa_sz);
	for (int i = 0; i < (wa_sz); ++i)
	{
		printf("%d\n",wa[i] );
	}
	int pos_count=0;
	int * pos = NP_periodic_pattern_matching(T, n, P_prime, m_prime, wa, &pos_count);
	for (int i = 0; i < (pos_count); ++i)
	{
		printf("%d ",pos[i] );
	}
	printf("\n");
	int m_u = p;
	char * u = (char *)malloc(sizeof(char)*m_u);
	for (int i = 0; i < m_u; ++i)
	{
		u[i] = P[i];
		// printf("%c ",u[i] );
	}
	// printf("\n");


	int k = m/p;

	int m_v = 0;
	char * v = (char *)malloc(sizeof(char)*m_v);
	for (int i = k*p; i < m; ++i)
	{
		v[i-k*p] = P[i];
		printf("%c ",v[i] );

		m_v += 1;
	}

	int * M = (int *)malloc(sizeof(int)*n);



	for (int i = 0; i < n; ++i)
	{
		M[i] = 0;
		int flag = 0;
		for(int j = 0; j<pos_count;j++){
			if (pos[j]==i){
				flag = 1;
				break;
			}
		}

		if (flag==1){
			int m_u2v = 2*m_u + m_v;
			char *u2v = makeu2v(u, m_u, v, m_v);
			// printf("u: %d v: %d\n",m_u,m_v);
			// for (int df = 0; df < m_u2v; ++df)
			// {
			// 	printf("%c ",u2v[df] );
			// }
		
			if(n - i >= m_u2v){  
				
				if (Occurs(T, i, u2v,m_u2v)){
					// printf("%s\n","Oh Yeah " );
					M[i] = 1;
				}
			}
			// return 0;
		}
	}
 	// check M
	// for (int i = 0; i < n; ++i)
	// {
	// 	printf("%d ",M[i] );
	// }
	
// ???????????????????????????????????CHECK n/p ////////////
	int ** S = (int **)malloc(sizeof(int *)*p);
	int ** C = (int **)malloc(sizeof(int *)*p);
	int * SC_sz = (int *)malloc(sizeof(int)*p);

	// for (int i = 0; i < p; ++i)
	// {
	// 	S[i] = (int *)malloc(sizeof(int)*n/p);
	// }

	for (int i = 0; i < p; ++i)
	{
		int j_s = 0;
		for (int jump = 0; jump < n; jump = jump+p)
		{
			
			j_s += 1;
		}
		S[i] = (int *)malloc(sizeof(int)*j_s);
		C[i] = (int *)malloc(sizeof(int)*j_s);

		j_s = 0;
		for (int jump = 0; jump < n; jump = jump+p)
		{
			S[i][j_s] = M[i + jump];
			j_s += 1;
		}
		SC_sz[i] = j_s;
		for (int j = 0; j < j_s; ++j)
		{
			C[i][j] = 0;
			////CHECK j ////////////////
			if(Ones(S[i],j_s,j,k-1)){
				C[i][j] = 1;
			}
		}
		

	}

	int * matched_positions = (int *)malloc(sizeof(int)*(n-m));
	for (int j = 0; j < n-m+1; ++j)
	{
		for (int i = 0; i < p; ++i)
		{
			for (int l = 0; l < SC_sz[i]; ++l)
			{
				if (j==i+l*p){
					matched_positions[j] = C[i][l];
				}
			}
		}
	}
	int match_count = 0;
	for (int i = 0; i < n-m; ++i)
	{
		if(matched_positions[i]==1){
			match_count += 1;
		}
	}

	int * match_positions = (int *)malloc(sizeof(int)*match_count);
	int index = 0;
	for (int i = 0; i < n-m, index < match_count; ++i)
	{
		if(matched_positions[i]==1){
			
			match_positions[index] = i;
			index+=1;
		}

	}


	// S check
	printf("\n");
	for (int i = 0; i < p; ++i)
	{
		for (int j = 0; j < n/p; ++j)
		{
			printf("%d ",C[i][j] );
		}
		printf("\nSC_sz: %d\n",SC_sz[i]);
	}

	for (int i = 0; i < match_count; ++i)
	{
		printf("%d ",match_positions[i] );
	}


	match_counts[0] = match_count;
	*matching_positions = match_positions;


	return 0;

}





// /*
// 	*****************************************************
// 		TODO -- You must implement this function
// 	*****************************************************
// */
void periodic_pattern_matching (
		int n, 
		char *text, 
		int num_patterns, 
		int *m_set, 
		int *p_set, 
		char **pattern_set, 
		int **match_counts, 
		int **matches)
{

	int numOfProc , rank;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numOfProc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank)
	match_counts[0] = (int *)malloc(sizeof(int)*num_patterns);

	for (int i = 0; i < num_patterns; ++i)
	{
		

		// int * wa = (int *)malloc(sizeof(int)*m_set[0]);
		int wa_sz = getPeriod(pattern_set[i],m_set[i]);//(m_set[0]+1)/2;

		int * wa = WITNESS(pattern_set[i], m_set[i], wa_sz);
		// for (int i = 0; i < (wa_sz); ++i)
		// {
		// 	printf("%d\n",wa[i] );
		// }

		int total_match= 0;
		int * matching_positions;
		if (wa_sz <= m_set[0]/2){
			P_periodic_pattern_matching(text,n,pattern_set[0],m_set[0],p_set[0],&matching_positions,&total_match);

		}
		else{ 
			NP_periodic_pattern_matching(text,n,pattern_set[0],m_set[0],wa,&total_match);
		}

		match_counts[0][i] = total_match;
 		// matches[i] = matching_positions;

 		printf("%s\n","Final Match" );
 		for (int j = 0; j < match_counts[0][i]; ++j)
		{
			printf("%d ",matching_positions[j] );
		}
	}





	MPI_Finalize();


	
	return;
}
