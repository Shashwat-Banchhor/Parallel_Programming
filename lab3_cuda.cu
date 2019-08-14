#include "lab3_cuda.h"





#include <iostream>
#include <cmath>
#include <malloc.h>
#include <fstream>

using namespace std;

#define TOLERANCE 0.001
#define JACOBI_UPDATE_TOLERANCE 0.001
// #define FILENAME "testcase_1000_300"
// #define samples 1000
// #define features 300

double **S; //Symmetric matrix (input)
double  *e; //eigenvalues
double **E; //eigenvectors
int  *ind;
bool *changed;
int  state;
int  N;

// void read_file(char* filename, int num_samples, int num_features, double** A) {
//     ifstream ifile;
//     ifile.open(filename, ios::in);

//     double tmp;
//     for (int i=0; i<num_samples; i++) {
//         for (int j=0; j<num_features; j++){
//             ifile >> tmp;
//             A[i][j] = tmp;
//         }
//     }

//     ifile.close();
// }

// #define BLKSIZ 4
#define BLOCK_SIZE 16

#define Blk_sz 16


struct Matrix
{
   int  width;
   int height;
   int stride;
   double * elements;
};

__device__ Matrix GetSubMatrix(Matrix A, int row, int col )
{
    Matrix A_sub;
    A_sub.width = Blk_sz;
    A_sub.height = Blk_sz;
    A_sub.stride = A.stride;
    A_sub.elements = & A.elements[A.stride*row*Blk_sz + col*Blk_sz];
    return A_sub;
}

__device__ double GetElement(Matrix A, int row, int col )
{
   

    return A.elements[A.stride*row*Blk_sz + col*Blk_sz];
}

__device__ void SetElement(Matrix A, int row, int col,double value )
{
   

    A.elements[A.stride*row*Blk_sz + col*Blk_sz] = value;
  
}


__global__ void MatMulKernel(Matrix A, Matrix B, Matrix C)
{
    // Block row and column
    int blockRow = blockIdx.y;
    int blockCol = blockIdx.x;
    // Each thread block computes one sub-matrix Csub of C
    Matrix Csub = GetSubMatrix(C, blockRow, blockCol);
    // Each thread computes one element of Csub
    // by accumulating results into Cvalue
    double Cvalue = 0;
    // Thread row and column within Csub
    int row = threadIdx.y;
    int col = threadIdx.x;


    int row_index = blockIdx.y * blockDim.y + threadIdx.y; 
    int col_index = blockIdx.x * blockDim.x + threadIdx.x;
    double sum = 0;
    

    if(blockIdx.y==0 and blockIdx.x==0 and threadIdx.x==1 and threadIdx.y==1){
        //  printf("A2:::%d::%d\n",blockRow,row);
        // for (int i = 0; i < A.height; ++i)
        // {
        //     for (int j = 0; j < A.width; ++j)
        //     {
        //         printf("%lf ",A.elements[i*A.width+j] );
        //     }
        //     printf("\n");
        // }

        // printf("B2:::::\n");
        // for (int i = 0; i < B.height; ++i)
        // {
        //     for (int j = 0; j < B.width; ++j)
        //     {
        //         printf("%lf ",B.elements[i*B.width+j] );
        //     }
        //     printf("\n");
        // }
    }



    if( col_index < B.width && row_index < A.height) 
    {
        for(int i = 0; i < A.width; i++) 
        {
            sum += A.elements[row_index * A.width + i] * B.elements[i * B.width + col_index];
            // printf("R: %d C: %d a: %lf(%d,%d) b : %lf(%d,%d) SUM: %lf\n",row_index,col_index,  A.elements[row_index * A.width + i],row_index,i, B.elements[i * B.width + col_index],i,col_index, sum);
        }
        // C will have same number of columns as B
        // printf("R%d C%d S%d\n",row_index, col_index,sum );
        C.elements[B.width * row_index + col_index] = sum;
    }
    __syncthreads();
   // Loop over all the sub-matrices of A and B that are
    // required to compute Csub
    // Multiply each pair of sub-matrices together
    // and accumulate the results
    return;
    // Loop over all the sub-matrices of A and B that are
    // required to compute Csub
    // Multiply each pair of sub-matrices together
    // and accumulate the results

    for (int m = 0; m < (A.width / BLOCK_SIZE); ++m) {
        // Get sub-matrix Asub of A
        Matrix Asub = GetSubMatrix(A, blockRow, m);
        // Get sub-matrix Bsub of B
        Matrix Bsub = GetSubMatrix(B, m, blockCol);
        // Shared memory used to store Asub and Bsub respectively
        __shared__ double As[BLOCK_SIZE][BLOCK_SIZE];
        __shared__ double Bs[BLOCK_SIZE][BLOCK_SIZE];
        // Load Asub and Bsub from device memory to shared memory
       
        // Each thread loads one element of each sub-matrix
        As[row][col] = GetElement(Asub, row, col);
        Bs[row][col] = GetElement(Bsub, row, col);
        // Synchronize to make sure the sub-matrices are loaded
        // before starting the computation
        __syncthreads();

            // Multiply Asub and Bsub together
        for (int e = 0; e < BLOCK_SIZE; ++e)
        Cvalue += As[row][e] * Bs[e][col];
        // Synchronize to make sure that the preceding
        // computation is done before loading two new
        // sub-matrices of A and B in the next iteration
        __syncthreads();
    }
    
    // Write Csub to device memory
    // Each thread writes one element
    SetElement(Csub, row, col, Cvalue);
 }





// Give transpose of A
__host__ void MatMulCu( Matrix A, Matrix B, Matrix C)
{
// Load A and B to device memory
    Matrix d_A;
    d_A.width = d_A.stride = A.width; d_A.height = A.height;
    size_t size = A.width * A.height * sizeof(double);
    cudaMalloc(&d_A.elements, size);
    cudaMemcpy(d_A.elements, A.elements, size,
    cudaMemcpyHostToDevice);
    Matrix d_B;
    d_B.width = d_B.stride = B.width; d_B.height = B.height;
    size = B.width * B.height * sizeof(double);
    cudaMalloc(&d_B.elements, size);
    cudaMemcpy(d_B.elements, B.elements, size,
    cudaMemcpyHostToDevice);



    Matrix d_C;
    d_C.width = d_C.stride = C.width; d_C.height = C.height;
    size = C.width * C.height * sizeof(double);
    cudaMalloc(&d_C.elements, size);
    // Invoke kernel
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
    dim3 dimGrid((B.width+BLOCK_SIZE-1) / BLOCK_SIZE, (A.height+BLOCK_SIZE -1) / BLOCK_SIZE);

    // printf("A:::::\n");
    // for (int i = 0; i < A.height; ++i)
    // {
    //     for (int j = 0; j < A.width; ++j)
    //     {
    //         printf("%lf ",A.elements[i*A.width+j] );
    //     }
    //     printf("\n");
    // }

    // printf("B:::::\n");
    // for (int i = 0; i < B.height; ++i)
    // {
    //     for (int j = 0; j < B.width; ++j)
    //     {
    //         printf("%lf ",B.elements[i*B.width+j] );
    //     }
    //     printf("\n");
    // }



    MatMulKernel<<<dimGrid, dimBlock>>>(d_A, d_B, d_C);
    // Read C from device memory
    cudaMemcpy(C.elements, d_C.elements, size,
    cudaMemcpyDeviceToHost);
    // Free device memory
    cudaFree(d_A.elements);
    cudaFree(d_B.elements);
    cudaFree(d_C.elements);


}

__host__ double** mat_transpose(double** A, int Am, int An) {
    double **B;
    B = (double**)malloc(__SIZEOF_POINTER__*An);
    for (int i=0; i<An; i++)
        B[i] = (double*)malloc(__SIZEOF_DOUBLE__*Am);

    for (int i=0; i<Am; i++){
        for (int j=0; j<An; j++){
            B[j][i] = A[i][j];
        }
    }

    return B;
}


// MUltiply matrices A and B C = A*B
__host__ double** mat_mul(double** A, int Am, int An, 
                 double** B, int Bm, int Bn)
{
    // double **C;
    // C = (double**)malloc(__SIZEOF_POINTER__*Am);
    // for (int i=0; i<Am; i++)
    //     C[i] = (double*)malloc(__SIZEOF_DOUBLE__*Bn);

    // for (int i=0; i<Am; i++){
    //     for (int j=0; j<Bn; j++){
    //         C[i][j] = 0;
    //         for (int k=0; k<An; k++){
    //             C[i][j] += A[i][k] * B[k][j];
    //         }
    //     }
    // }

    // return C;

     double **C;
    C = (double**)malloc(__SIZEOF_POINTER__*Am);
    for (int i=0; i<Am; i++)
        C[i] = (double*)malloc(__SIZEOF_DOUBLE__*Bn);

    double *A_cu = (double *)malloc(sizeof(double)*Am*An);
    double *B_cu = (double *)malloc(sizeof(double)*Bm*Bn);    
    double *C_cu = (double *)malloc(sizeof(double)*Am*Bn);    

    // printf("A in mat_mul \n");
    for (int i = 0; i < Am; ++i)
    {
        for (int j = 0; j < An; ++j)
        {
            A_cu[i*An + j] = A[i][j];
            // printf("%lf ", A_cu[i*An + j]);
        }
        // printf("\n");
    }

    // printf("B in mat_mul \n");

    for (int i = 0; i < Bm; ++i)
    {   
        for (int j = 0; j < Bn; ++j)
        {
           B_cu[i*Bn+j] = B[i][j];
            // printf("%lf ", B_cu[i*Bn + j]);

        }
        // printf("\n");
    }

    Matrix CudaA ;
    CudaA.width = An;
    CudaA.height  = Am;
    CudaA.stride = An;
    // CudaA.actual_ht = Am;
    CudaA.elements = A_cu;//A[0];

    Matrix CudaB ;
    CudaB.width = Bn;
    CudaB.height  = Bm;
    CudaB.stride = Bn;
    // CudaB.actual_ht = Bm;
    CudaB.elements =  B_cu;//B[0];

    Matrix CudaC ;
    CudaC.width = Bn;
    CudaC.height  = Am;
    CudaC.stride = Bn;
    // CudaC.actual_ht = Am;
    CudaC.elements = C_cu;


    MatMulCu(CudaA,CudaB,CudaC);



    // C = CUdaC.elements;
    // printf("CUda out :\n");
    for (int i=0; i<Am; i++){
        for (int j=0; j<Bn; j++){
            C[i][j] = CudaC.elements[i*Bn+j];
            // printf("%lf ",C[i][j] );
        }
        // printf("\n");
    }

    return C;
}


// MUltiply matrices A and B C = A*B
__host__ double** mat_mul_o(double** A, int Am, int An, 
                 double** B, int Bm, int Bn)
{
    

    double **C;
    C = (double**)malloc(__SIZEOF_POINTER__*Am);
    for (int i=0; i<Am; i++)
        C[i] = (double*)malloc(__SIZEOF_DOUBLE__*Bn);

    for (int i=0; i<Am; i++){
        for (int j=0; j<Bn; j++){
            C[i][j] = 0;
            for (int k=0; k<An; k++){
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    return C;



}


int maxind(int k) {
    int m = k+1;

    for (int i = k+2; i < N; i++){
        if (fabs(S[k][i]) > fabs(S[k][m])){
            m = i;
        }
    }

    return m;
}

__host__ void update(int k, double t) {
    double ek_prev = e[k];
    e[k] = ek_prev + t;

    if (e[k] < 0) e[k] = 0;

    if (changed[k] && fabs(ek_prev - e[k]) < JACOBI_UPDATE_TOLERANCE) {
        changed[k] = false;
        state = state - 1;
    }
    else if ((! changed[k]) && fabs(ek_prev - e[k]) > JACOBI_UPDATE_TOLERANCE) {
        changed[k] = true;
        state = state + 1;
    }
}

__host__ void rotate(int k, int l, int i, int j, double c, double s,
            bool eigenvectors)
{
    
//   We must change the rotate function in order to reduce the 
    // overhead of creating 2-d arrays multiple times



    double* mat1;//double** mat1;
    double* mat2;//double** mat2;
    double* mat3;//double** mat3;

    mat1 = (double*)malloc(sizeof(double)*4);
    // mat1 = (double*)malloc(__SIZEOF_POINTER__*4);
    
    // mat1[0] = (double*)malloc(__SIZEOF_DOUBLE__*2);
    // mat1[1] = (double*)malloc(__SIZEOF_DOUBLE__*2);
    mat1[0]=c;// mat1[0][0] = c; 
    mat1[1]=-s;//mat1[0][1] = -s;
    mat1[2]=s;// mat1[1][0] = s;
    mat1[3]=c;//mat1[1][1] = c;



    mat2 = (double*)malloc(sizeof(double)*2);
    // mat2[0] = (double*)malloc(__SIZEOF_DOUBLE__*1);
    // mat2[1] = (double*)malloc(__SIZEOF_DOUBLE__*1);
    if (eigenvectors){
        mat2[0] = E[i][k];//mat2[0][0] = E[i][k];
        mat2[1] = E[i][l];//mat2[1][0] = E[i][l];
    }
    else {
        mat2[0] = S[k][l];//mat2[0][0] = S[k][l];
        mat2[1] = S[i][j];//mat2[1][0] = S[i][j];
    }

    // mat3 = mat_mul_o(mat1, 2, 2, mat2, 2, 1);

    // Its only a 2x2 matrix multiply it directly
    mat3 = (double*)malloc(sizeof(double)*2);
    mat3[0] = mat1[0]*mat2[0] + mat1[1]*mat2[1];
    mat3[1] = mat1[2]*mat2[0] + mat1[3]*mat2[1];



    if (eigenvectors){
        E[i][k] = mat3[0];//mat3[0][0];
        E[i][l] = mat3[1];//mat3[1][0];
    }
    else{
        S[k][l] = mat3[0];//mat3[0][0];
        S[i][j] = mat3[1];//mat3[1][0];
    }

    // free(mat1[0]);
    // free(mat1[1]);
    free(mat1);
    // free(mat2[0]);
    // free(mat2[1]);
    free(mat2);
    // free(mat3[0]);
    // free(mat3[1]);
    free(mat3);
}


// print the matrix
__host__ void print_matrix(double** A, int Am, int An) {
    cout << "[";
    for (int i=0; i<Am; i++){
        if (i>0)
            cout<<" ";
        cout<<"[";
        for (int j=0; j<An-1; j++){
            cout << A[i][j] << ", ";
        }
        if (i < Am-1)
            cout << A[i][An-1] << "]" << endl;
    }
    cout << A[Am-1][An-1] << "]]" << endl;
}


// print the vector
__host__ void print_vector(double* A, int An) {
    cout << "[";
    for(int i=0; i<An-1; i++)
        cout << A[i] << ",";
    cout << A[An-1] << "]" << endl;
}


// initialize Jacobi 
    // --  malloc e and E and changed
__host__ void init_jacobi() {
    E = (double**)malloc(__SIZEOF_POINTER__*N);
    for (int i=0; i<N; i++){
        E[i] = (double*)malloc(__SIZEOF_DOUBLE__*N);
        for (int j=0; j<N; j++){
            E[i][j] = 0;
        }
        E[i][i] = 1;
    }

    state = N;

    e = (double*)malloc(__SIZEOF_DOUBLE__*N);
    ind = (int*)malloc(__SIZEOF_INT__*N);
    changed = (bool*)malloc(sizeof(bool)*N);

    for (int k=0; k<N; k++){
        ind[k]     = maxind(k);
        e[k]       = S[k][k];
        changed[k] = true;
    }
}



// Given a Matrix gives its eigenvalues and eigenvectors
__host__ void Jacobi(double **input_matrix, int n, 
            double **eigenvalues, double ***eigenvectors) 
{
    N = n;
    S = input_matrix;

    init_jacobi();

    while(state != 0){
        int m = 0;

        for (int k=1; k<N-1; k++){
            if (fabs(S[k][ind[k]]) > fabs(S[m][ind[m]])){
                m = k;
            }
        }

        int k = m;
        int l = ind[m];
        double p = S[k][l];
        double y = (e[l] - e[k]) / 2.0;
        double d = fabs(y) + sqrt(p*p + y*y);
        double r = sqrt(p*p + d*d);
        double c = d / r;
        double s = p / r;
        double t = (p*p) / d;

        if (y < 0.0) { s = -s; t = -t; }

        S[k][l] = 0.0;
        update(k, -t);
        update(l, t);

        for (int i=0; i<k; i++)  { rotate(i, k, i, l, c, s, false); }
        for (int i=k+1; i<l; i++){ rotate(k, i, i, l, c, s, false); }
        for (int i=l+1; i<N; i++)  { rotate(k, i, l, i, c, s, false); }

        for (int i=0; i<N; i++){
            rotate(k, l, i, i, c, s, true);
        }

        ind[k] = maxind(k);
        ind[l] = maxind(l);
    }

    *eigenvalues = e;
    *eigenvectors = E;
}



















///////////////////////////   LAB 2 COPY ////////////////////////////
#include <iostream>
#include "stdlib.h"
// #include "lab2_omp.h"
// #include <omp.h>
#include <math.h>

using namespace std;


struct EigenValIdx{
    double value;
    int index;
};

// Transpose M
__host__ void transpose(double * M, int rows, int cols,double ** M_T){

    // #pragma omp parallel
    // {   
        // #pragma omp for collapse(2)
        for (int i = 0; i < cols; ++i)
        {
            for(int j=0 ; j< rows ; j++){

                // D_T[i][j] = D[j][i]
                M_T[0][i*rows + j] = M[j*cols+i]; 

            }
        }
    // }

    return;
}





__host__ void tiling_multilply_matrices(double * M1 ,int rows1,int cols1,double * M2,int rows2,int cols2, double ** M3){
    
    int BLKSIZ = 16;

    double* M2_T = (double *) malloc(sizeof(double) * cols2 * rows2);
    transpose(M2,rows2,cols2,&M2_T);
    
    // int blkYA = rows1/BLKSIZ;
    // int blkYB = cols2/BLKSIZ;
    // int blkX = cols1/BLKSIZ;
    
    for (int i = 0; i < rows1; ++i)
    {
        for (int j= 0; j < cols2; ++j)
        {
            M3[0][i* cols2 + j]  = 0.0;
            for(int t = 0; t < cols1/BLKSIZ+1;t++)
            {
                for (int k = 0; k < 16 && t*BLKSIZ+k < cols1; ++k)
                {
                    // cerr <<                D[i*cols + k] << " "<<D_T[k*tcols + j]<<" "<<D[i*cols + k] *D_T[k*tcols + j] <<endl;
                    // cerr << DD_T[i*tcols + j] <<endl;
                    // cerr <<"Is : "<< 1*1 ;
                    // cerr <<"Multiplication is " << D[i*cols + k] * D_T[k*tcols + j] <<endl;
                    // M3[0][i*cols2 + j]  += M1[i*cols1 + k] * M2[k*cols2 + j];
                    M3[0][i*cols2 + j]  += M1[i*cols1 + t*BLKSIZ + k] * M2_T[j*rows2 + t*BLKSIZ + k];

                    // cerr << DD_T[i*tcols + j] <<endl;
                }
            }
        }
    }


    free(M2_T);
    return;
    
}

// Multiply M1 and M2 and stores them in M3
__host__ void multilply_matrices(double * M1 ,int rows1,int cols1,double * M2,int rows2,int cols2, double ** M3){
    

    int BLKSIZ = 16;

    double* M2_T = (double *) malloc(sizeof(double) * cols2 * rows2);
    transpose(M2,rows2,cols2,&M2_T);
    
    // int blkYA = rows1/BLKSIZ;
    // int blkYB = cols2/BLKSIZ;
    // int blkX = cols1/BLKSIZ;
    #pragma omp parallel
    {
        #pragma omp for collapse(2)
        for (int i = 0; i < rows1; ++i)
        {
            for (int j= 0; j < cols2; ++j)
            {
                M3[0][i* cols2 + j]  = 0.0;
                for(int t = 0; t < cols1/BLKSIZ+1;t++)
                {
                    for (int k = 0; k < 16 && t*BLKSIZ+k < cols1; ++k)
                    {
                        // cerr <<                D[i*cols + k] << " "<<D_T[k*tcols + j]<<" "<<D[i*cols + k] *D_T[k*tcols + j] <<endl;
                        // cerr << DD_T[i*tcols + j] <<endl;
                        // cerr <<"Is : "<< 1*1 ;
                        // cerr <<"Multiplication is " << D[i*cols + k] * D_T[k*tcols + j] <<endl;
                        // M3[0][i*cols2 + j]  += M1[i*cols1 + k] * M2[k*cols2 + j];
                        M3[0][i*cols2 + j]  += M1[i*cols1 + t*BLKSIZ + k] * M2_T[j*rows2 + t*BLKSIZ + k];

                        // cerr << DD_T[i*tcols + j] <<endl;
                    }
                }
            }
        }
    }

    free(M2_T);
    return;

}

__host__ void multilply_matrices_diagonal(double * M1 ,int rows1,int cols1,double * M2,int rows2,int cols2, double ** M3,double ** E){
    


    int BLKSIZ = 16;

    double* M2_T = (double *) malloc(sizeof(double) * cols2 * rows2);
    transpose(M2,rows2,cols2,&M2_T);
    
    // int blkYA = rows1/BLKSIZ;
    // int blkYB = cols2/BLKSIZ;
    // int blkX = cols1/BLKSIZ;
    #pragma omp parallel
    {
        #pragma omp for collapse(2)
        for (int i = 0; i < rows1; ++i)
        {
            for (int j= 0; j < cols2; ++j)
            {
                M3[0][i* cols2 + j]  = 0.0;
                if (i==j){E[0][i*cols2+j] = 1.0;}
                
                else{E[0][i*cols2+j] = 0.0;}
                for(int t = 0; t < cols1/BLKSIZ+1;t++)
                {
                    for (int k = 0; k < 16 && t*BLKSIZ+k < cols1; ++k)
                    {
                        // cerr <<                D[i*cols + k] << " "<<D_T[k*tcols + j]<<" "<<D[i*cols + k] *D_T[k*tcols + j] <<endl;
                        // cerr << DD_T[i*tcols + j] <<endl;
                        // cerr <<"Is : "<< 1*1 ;
                        // cerr <<"Multiplication is " << D[i*cols + k] * D_T[k*tcols + j] <<endl;
                        // M3[0][i*cols2 + j]  += M1[i*cols1 + k] * M2[k*cols2 + j];
                        M3[0][i*cols2 + j]  += M1[i*cols1 + t*BLKSIZ + k] * M2_T[j*rows2 + t*BLKSIZ + k];

                        // cerr << DD_T[i*tcols + j] <<endl;
                    }
                }
            }
        }
    }

    free(M2_T);
    return;

}
__host__ void multilply_matrices_together(double * M1 ,int rows1,int cols1,double * M2,int rows2,int cols2, double ** M3, double * M4 , double * M5, double ** M6){
    


    int BLKSIZ = 16;

    double* M2_T = (double *) malloc(sizeof(double) * cols2 * rows2);
    transpose(M2,rows2,cols2,&M2_T);
    double* M5_T = (double *) malloc(sizeof(double) * cols2 * rows2);
    transpose(M5,rows2,cols2,&M5_T);
    // int blkYA = rows1/BLKSIZ;
    // int blkYB = cols2/BLKSIZ;
    // int blkX = cols1/BLKSIZ;
    #pragma omp parallel
    {
        #pragma omp for collapse(2)
        for (int i = 0; i < rows1; ++i)
        {
            for (int j= 0; j < cols2; ++j)
            {
                M3[0][i* cols2 + j]  = 0.0;
                M6[0][i* cols2 + j]  = 0.0;
                for(int t = 0; t < cols1/BLKSIZ+1;t++)
                {
                    for (int k = 0; k < 16 && t*BLKSIZ+k < cols1; ++k)
                    {
                        // cerr <<                D[i*cols + k] << " "<<D_T[k*tcols + j]<<" "<<D[i*cols + k] *D_T[k*tcols + j] <<endl;
                        // cerr << DD_T[i*tcols + j] <<endl;
                        // cerr <<"Is : "<< 1*1 ;
                        // cerr <<"Multiplication is " << D[i*cols + k] * D_T[k*tcols + j] <<endl;
                        // M3[0][i*cols2 + j]  += M1[i*cols1 + k] * M2[k*cols2 + j];
                        M3[0][i*cols2 + j]  += M1[i*cols1 + t*BLKSIZ + k] * M2_T[j*rows2 + t*BLKSIZ + k];
                        M6[0][i*cols2 + j]  += M4[i*cols1 + t*BLKSIZ + k] * M5_T[j*rows2 + t*BLKSIZ + k];


                        // cerr << DD_T[i*tcols + j] <<endl;
                    }
                }
            }
        }
    }

    free(M2_T);
    return;

}




// Print Matrix M
__host__ void display_matrix(double * M,int rows,int cols){

    for (int i = 0; i < rows*cols; ++i)
    {
        if (i%cols==0){cerr << endl;}
        cerr << M[i] <<" ";
    }

    cerr << endl;
}

//    /
//  \/
__host__ void inner_product(double * A , double * B ,int n, double * ans){

    ans[0] = 0.0;
    // #pragma omp parallel for reduction(+:ans[0])
    for (int i = 0; i < n; ++i)
    {
        ans[0] += A[i]*B[i];        
    }
}

//                                         /
// Projection of Vector A on Vector U    \/
__host__ void projection(double * U , double * A ,int n, double ** P){

    double num,den;
    inner_product(A,U,n,&num);
    inner_product(U,U,n,&den);

    // #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i)
    {
        P[0][i] = (num/den)*U[i];
    }

    return;
}
//                                                   /
// Normalize ColNo row of Matrix U with n columns  \/
__host__ void normalize(double ** U,int ColNo , int n){
    
    double IIeII = 0.0;
    
    // Norm of E
    // cerr << "Row "<< ColNo <<endl;
    for (int i = 0; i < n; ++i)
    {
        // cerr << U[0][ColNo*n + i] <<" "<<U[0][ColNo*n + i] <<endl; 
        IIeII += U[0][ColNo*n + i] * U[0][ColNo*n + i];
    }
    IIeII = sqrt(IIeII);
    // cerr << IIeII;
    //

    // Update E
    // #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i)
    {
        U[0][ColNo*n + i] = U[0][ColNo*n + i]/IIeII;    
    }
    // 

    return;
}

// Normalize ColNo row of Matrix U with n columns and paste in A  \/
__host__ void normalize2(double ** U,double ** A,int ColNo , int n){
    
    double IIeII = 0.0;
    
    // Norm of E
    // cerr << "Row "<< ColNo <<endl;
    for (int i = 0; i < n; ++i)
    {
        // cerr << U[0][ColNo*n + i] <<" "<<U[0][ColNo*n + i] <<endl; 
        IIeII += U[0][ColNo*n + i] * U[0][ColNo*n + i];
    }
    IIeII = sqrt(IIeII);
    // cerr << IIeII;
    //

    // Update E
    for (int i = 0; i < n; ++i)
    {
        A[0][ColNo*n + i] = U[0][ColNo*n + i]/IIeII;    
    }
    // 

    return;
}
//   /                                                                  
// \/ Assign ColNo Row of A to ColNo row of U : Name is Problem Specific 
__host__ void assign(double ** U , int ColNo, double * A , int n){

    
    for (int i = 0; i < n; ++i)
    {
        U[0][ColNo*n+i] = A[ColNo*n+i];
    }
    return;
}

  
// \/ Assign ColNo Row of A to ColNo row of U : Name is Problem Specific
// n_U is < n_A  
__host__ void assign2(double ** U , int ColNo, double * A , int n_U,int n_A){

    
    for (int i = 0; i < n_U; ++i)
    {
        U[0][ColNo*n_U+i] = A[ColNo*n_A+i];
    }
    return;
}
//                                                              /
// U is Matrix From ColNo row vector of U subtract P_i vector \/
__host__ void subtract(double ** U,int ColNo , int n ,  double * P_i){


    for (int i = 0; i < n; ++i)
    {
        U[0][ColNo*n + i] = U[0][ColNo*n + i] - P_i[i];
    }
    return;
}


//   /
// \/  Make E a diagonal Matrix with all diagonal entries 1.0
__host__ void makeDiagonal(double ** E , int rows , int cols, double n){

    // #pragma omp parallel
    // {    
        #pragma omp for collapse(2)
        for (int i = 0; i < rows; ++i)
        {
            for(int j =0;j<cols;j++){
                
                if (i==j){E[0][i*cols+j] = n;}
                
                else{E[0][i*cols+j] = 0.0;}
            }
        }
    // }
}


//   /
// \/
__host__ void copyMatrixTogether(double ** E,double * D,int rows, int cols,double ** F, double * G){
    
    #pragma omp parallel
    {
        #pragma omp for collapse(2)
        for (int i = 0; i < rows; ++i)
        {
            for (int j = 0; j < cols; ++j)
            {
                E[0][i*cols+j] = D[i*cols+j];
                F[0][i*cols+j] = G[i*cols+j];

            }
        
        }
    }
    return;
}

__host__ void copyMatrix2(double ** E,double * D,int rows, int cols1, int cols2){

        for (int i = 0; i < rows; ++i)
        {
            for (int j = 0; j < cols1; ++j)
            {
                E[0][i*cols1+j] = D[i*cols2+j];
            }
        }

    // }    
    return;
}


__host__ void GramSchmidt(double * A,int rows , int cols, double ** Q , double ** R){

    // This is the transpose of matrix M as we had to deal in column major form in M we are
    // goig to deal in row major form in A

    // This is the transpose of A      [ --- a1----]
                                    // | --- a2--- |
                                    // | --- . --- |
                                    // | --- . --- |
                                    // | --- . --- |
                                    // [ --- an ---]
    // cerr << "    A "<<endl;
    
    // display_matrix(A,rows,cols);
    double* A_T = (double *) malloc(sizeof(double) * cols * rows);
    transpose(A,rows,cols,&A_T);


    // This is the transpose of U U_T =    [ --- u1----]   
                                        // | --- u2--- |
                                        // | --- . --- |
                                        // | --- . --- |
                                        // | --- . --- |
                                        // [ --- un ---]

    double* U_T = (double *) malloc(sizeof(double) * cols * rows);
    double* Q_T = (double *) malloc(sizeof(double) * cols * rows);
    
    //WIKIPEDIA 
    for (int i = 0; i < cols; ++i)
    {
        
        // Assign the ith row of U_T to A_T
        // u_i = a_i
        assign(&U_T,i,A_T,rows);
        
        
        
        double* P_i = (double *) malloc(sizeof(double) * rows);
        
        for(int j=0;j<i;j++){
            
            
            projection(U_T+j*rows, A_T+i*rows ,rows, &P_i);

            subtract(&U_T,i,rows,P_i);
            
        }
        
    }


    //  This is the transpose of Q Q_T =   [ --- e1----]    where ei = ui/||ui||
                                        // | --- e2--- |
                                        // | --- . --- |
                                        // | --- . --- |
                                        // | --- . --- |
                                    // [ --- en ---]
    // cerr << "U_T : "<<endl; 
    // display_matrix(U_T,cols,rows);
    for (int i = 0; i < cols; ++i)
    {
        // update Q_T with ei = ui/||ui||
        normalize2(&U_T,&Q_T,i,rows);
        
    }
    // cerr << "Q_T : "<<endl; 

    // display_matrix(Q_T,cols,rows);
    transpose(Q_T,cols,rows,Q);

    // R = Q_TA
    makeDiagonal(R,cols,cols,0.0);
    for (int i = 0; i < cols; ++i)
    {
        double inr_pdt = 0.0; 
        for(int j=i;j<cols;j++){
            
            inner_product(Q_T+i*rows,A_T + j*rows,rows,&inr_pdt);
            R[0][i*cols + j] = inr_pdt;
            inr_pdt = 0.0;
        }
    }
    // multilply_matrices(U_T,cols,rows,A,rows,cols,R);

    

}

__host__ void ModifiedGramSchmidt(double * A,int rows , int cols, double ** Q , double ** R){

    // This is the transpose of matrix M as we had to deal in column major form in M we are
    // goig to deal in row major form in A

    // This is the transpose of A      [ --- a1----]
                                    // | --- a2--- |
                                    // | --- . --- |
                                    // | --- . --- |
                                    // | --- . --- |
                                    // [ --- an ---]
    // cerr << "    A "<<endl;
    
    // display_matrix(A,rows,cols);
    double* A_T = (double *) malloc(sizeof(double) * cols * rows);
    transpose(A,rows,cols,&A_T);


    // This is the transpose of U U_T =    [ --- u1----]   
                                        // | --- u2--- |
                                        // | --- . --- |
                                        // | --- . --- |
                                        // | --- . --- |
                                        // [ --- un ---]

    double* U_T = (double *) malloc(sizeof(double) * cols * rows);
    
    //WIKIPEDIA 
                            // for (int i = 0; i < cols; ++i)
                            // {
                                
                            //  // Assign the ith row of U_T to A_T
                            //  // u_i = a_i
                            //  assign(&U_T,i,A_T,rows);
                                
                                
                                
                            //  double* P_i = (double *) malloc(sizeof(double) * rows);
                                
                            //  for(int j=0;j<i;j++){
                                    
                            //      // Projection of A_i on U_j 
                            //      // cerr << "jth sub-new: "<<j<<endl;
                            //          // for (int k = 0; k < rows; ++k)
                            //          // {
                            //          //  cerr << "       "<<U_T[i*rows + k]<<" ";
                            //          // }
                            //          // cerr << endl;
                            //      projection(U_T+j*rows, U_T+i*rows ,rows, &P_i);

                            //      subtract(&U_T,i,rows,P_i);
                            //      // if(i==3){
                            //          // cerr << "jth sub-row: "<<j<<endl;
                            //          // for (int k = 0; k < rows; ++k)
                            //          // {
                            //          //  cerr << "       "<<U_T[j*rows + k]<<" ";
                            //          // }
                            //          // cerr << "jth sub-proj: "<<j<<endl;
                            //          // for (int k = 0; k < rows; ++k)
                            //          // {
                            //          //  cerr << "       "<<P_i[k]<<" ";
                            //          // }
                            //          // cerr << endl;
                            //          // cerr << "jth sub-subt: "<<j<<endl;
                            //          // for (int k = 0; k < rows; ++k)
                            //          // {
                            //          //  cerr << "       "<<U_T[i*rows + k]<<" ";
                            //          // }
                            //          // cerr << endl;
                            //      // }
                            //  }
                            //  normalize(&U_T,i,rows);
                            //  // cerr << "Ith COL: "<<i<<endl;
                            //  //      for (int k = 0; k < rows; ++k)
                            //  //      {
                            //  //          cerr << "   "<<U_T[i*rows + k]<<" ";
                            //  //      }
                            //  //      cerr << endl;

                            // }
                            // //   This is the transpose of Q Q_T =   [ --- e1----]    where ei = ui/||ui||
                            //                                  // | --- e2--- |
                            //                                  // | --- . --- |
                            //                                  // | --- . --- |
                            //                                  // | --- . --- |
                            //                                  // [ --- en ---]
                            


                            // transpose(U_T,cols,rows,Q);
                            // // R = Q_TA
                            // multilply_matrices(U_T,cols,rows,A,rows,cols,R);

                                // // NEW UPDT
    // for (int i = 0; i < cols; ++i)
    // {
    //  //  // Assign the ith row of U_T to A_T
    //      //  // u_i = a_i
    //  assign(&U_T,i,A_T,rows);
    // }


    // REMOVED THIS AND ADDED OMP
    // copyMatrix(&U_T, A_T,rows, cols);
    // makeDiagonal(R,cols,cols,0.0);

    #pragma omp parallel for collapse(2)
    for (int i = 0; i < rows; ++i)
    {
        for(int j=0;j<rows;++j)
        {
            U_T[i*rows + j] = A_T[i*rows+j];
            R[0][i*rows + j ] = 0.0;
        }
    }

    double* P_i = (double *) malloc(sizeof(double) * rows);
    double* Q_T = (double *) malloc(sizeof(double) * cols * rows);
    
    // 
    for (int i = 0; i < cols; ++i)
    {
        double norm = 0.0;
        inner_product(U_T+i*rows, U_T+i*rows, rows, &norm);
        norm = sqrt(norm);
        R[0][i*cols + i] = norm;
        for(int k=0;k<rows;k++){
            Q_T[i*rows+k] = (1.0/norm)*U_T[i*rows + k];
        }
        
        for (int j = i+1; j < rows; ++j)
        {
            // projection(Q[0]+i*rows, U_T+k*rows ,rows, &P_i);
            double dot_p = 0.0;
            inner_product(Q_T+i*rows,U_T+j*rows,rows,&dot_p);
            R[0][i*cols + j] = dot_p;
            for(int k=0;k<rows;k++){
                P_i[k] = dot_p*Q_T[i*rows+k];
            }
            subtract(&U_T,j,rows,P_i);
        }
    }
    transpose(Q_T,cols,rows,Q);
    // cerr << "    Q"<<endl;

    // display_matrix(Q[0],rows,cols);
    // cerr << "    R "<<endl;
    // display_matrix(R[0],rows,cols);
}






// Check whether A is upper triangular or not
__host__ bool convergence(double * A , double * A_old, int rows, int cols){
    double epsilon = 0.001;
    double max_diff = 0.0;
    for (int i = 0; i < rows; ++i)
    {
        // for(int j =0;j<rows;j++){
        int j = 1;
            // if (abs(A[i*cols+j] - A_old[i*cols+j])>epsilon){ cerr << abs(A[i*cols+j] - A_old[i*cols+j]) << endl; return false;}
            if (fabs(A[i*cols+j] - A_old[i*cols+j])>max_diff){  max_diff = fabs(A[i*cols+j] - A_old[i*cols+j]);}

        // }
    }
    if (max_diff > epsilon){//cerr << max_diff<<" "<<epsilon << " Diff: " << max_diff - epsilon<<endl; 
    return false;}

    return true;
}

__host__ int  compare(const void * A , const void * B){
    double a = *((double *) A);
    double b = *((double *) B);
    if (a<b){return 1;}
    else{ return -1;}
}

__host__ int  compare_struct(const void * A , const void * B){
    EigenValIdx a = *((EigenValIdx *) A);
    EigenValIdx b = *((EigenValIdx *) B);
    if (a.value<b.value){return 1;}
    else{ return -1;}
}

__host__ void check_product(double * U , double * SIGMA, double * V_T, double * D_T, int M, int N){

    double* SIGMA_MAT= (double *) calloc( N * M, sizeof(double));
    for (int i = 0; i < N; ++i)
    {
        SIGMA_MAT[i*M + i] = SIGMA[i];
    }

    double* USIGMA_MAT= (double *) calloc( N * M, sizeof(double));
    double* OUTPUT= (double *) calloc( N * M, sizeof(double));
    double* OUTPUT_T= (double *) calloc( N * M, sizeof(double));

    multilply_matrices(U, N, N, SIGMA_MAT, N, M, & USIGMA_MAT);
    multilply_matrices(USIGMA_MAT, N, M, V_T, M, M,&OUTPUT);
    // display_matrix(U, N, N);
    // display_matrix(SIGMA_MAT, N, M );
    // display_matrix(V_T, M, M);
    transpose(OUTPUT,N,M,&OUTPUT_T);
    display_matrix(OUTPUT_T, M, N);

}


// void read_matrix (const char* input_filename, int* M, int* N, double** D){
//  FILE *fin = fopen(input_filename, "r");

//  fscanf(fin, "%d%d", M, N);
    
//  int num_elements = (*M) * (*N);
//  *D = (double*) malloc(sizeof(double)*(num_elements));
    
//  for (int i = 0; i < num_elements; i++){
//      fscanf(fin, "%f", (*D + i));
//  }
//  fclose(fin);
// }



__host__ void copytheOutput(double ** D,double ** D_HAT, int rows, int cols){
    for (int i = 0; i < rows; ++i)
    {
       for(int j=0;j<cols;j++){
         // printf("i %d j %d\n",i,j );
         D_HAT[0][i*cols+j] = D[i][j];
       }
    }
    return;
}


__host__ double ** make2D(double * W , int rows , int cols){
    double ** W2D ; 


    W2D = (double**)malloc(__SIZEOF_POINTER__*rows);
    for (int i=0; i<rows; i++)
        W2D[i] = (double*)malloc(__SIZEOF_DOUBLE__*cols);

    for (int i=0; i<rows; i++){

        for (int j=0; j<cols; j++){

            W2D[i][j] = W[i*cols+j];
            
        }
    }

    return W2D ;
}

// /*
//  *****************************************************
//      TODO -- You must implement this function
//  *****************************************************
// */
// Calculates the SVD of Transpose of D i.e. D_T
__host__ void SVD_T(int M, int N, double* D, double** U, double** SIGMA, double** V_T){
    int BLKSIZ = 16;
    // Display Input
    // cerr << "Matrix D : ";
    // cerr << M<<"x"<< N ;
    // display_matrix(D,M,N);
    // exit(0);
    
    // Calculate transpose of D
    
    // double* D_T = (double *) malloc(sizeof(double) * N * M);
    // transpose(D,M,N,&D_T);

    double ** D_2d = make2D(D,M,N);
    double ** D_T = mat_transpose(D_2d,M,N);


    // cerr << "Calculating SVD of D_T : ";
    // display_matrix(D_T, N, M);
    
    // Calculate D multiplied by D_T
    // as we have to calculate SVD  of D of D_T


    // double* DD_T = (double *) malloc(sizeof(double) * M * M) ;
    // multilply_matrices(D, M, N, D_T, N, M, &DD_T);
    double ** DD_T = mat_mul(D_2d,M,N,D_T,N,M);    



    // display_matrix(DD_T, M, M);
    
    // Memory allocation for matrices
    // double* Q = (double *) malloc(sizeof(double) * M * M) ;
    // double* R = (double *) malloc(sizeof(double) * M * M) ;
    // double* Q_T = (double *) malloc(sizeof(double) * M * M) ;
    
    // // D0
    // double* EValues = (double *) malloc(sizeof(double) * M * M)  ;
    // double* EValues_old = (double *) malloc(sizeof(double) * M * M)  ;

    // // E0
    // double* EVectors = (double *) malloc(sizeof(double) * M * M) ;
    // double* EVectors_old = (double *) malloc(sizeof(double) * M * M)  ;
    // double* Check = (double *) malloc(sizeof(double) * M * M) ;
    // double* Check_inp = (double *) malloc(sizeof(double) * M * M) ;

    // Initialize 
    // D0 = D
    // copyMatrix(&EValues , DD_T, M, M);
// multilply_matrices_diagonal(D, M, N, D_T, N, M, &EValues,&EVectors);

    // E0 = I
    // makeDiagonal(&EVectors, M, M,1.0);
    // cerr << "This is DD_T: "<<endl;
    // display_matrix(DD_T, M, M);
   

    ///////////////////////////////////// REPLACE BY JACOBI ///////////
// int iter = 0;
// do{
    
//     // Store old values for convergence check
//     // copyMatrix(&EValues_old, EValues, M, M);
//     copyMatrixTogether(&EVectors_old, EVectors, M, M,&EValues_old,EValues);

//     // Q-R Decomposition
//     // We get Qi,Ri
//     ModifiedGramSchmidt(EValues, M, M, &Q , &R);
//     // GramSchmidt(EValues, M, M, &Q , &R);
    
//     // cerr << "Q"<<endl;
//     // display_matrix(Q,M,M);
//     // cerr << "R"<<endl;
//     // display_matrix(R,M,M);
//     // cerr << "Check "<<iter<<endl;
//     //  multilply_matrices(Q, M, M, R, M, M,&Check);
//     //  display_matrix(Check,M,M);
//     //  cerr << "Next_INPUT:" <<endl;
//     //  multilply_matrices(R, M, M, Q, M, M,&Check_inp);
//     //  display_matrix(Check_inp,M,M);
//     // exit(0);
//     // Multiply Ri,Qi
//     // Di+1 = Ri * Qi
//     // transpose(Q,M,M,&Q_T);
//     multilply_matrices_together(R, M, M, Q, M, M,&EValues,EVectors_old,Q,&EVectors);

//     // Multiply Ei,Qi
//     // Ei+1 = Ei*Qi
//     // multilply_matrices(EVectors_old, M, M, Q, M, M,&EVectors);
    

    
//     iter+= 1;
    
// }while(!(convergence(EValues, EValues_old, M, M) || iter >30));// || convergence(EVectors, EVectors_old, M, M)) );

//     cudaEvent_t start, stop;

// cudaEventCreate(&start);
//     cudaEventCreate(&stop);

//     float computation_time;

    double * eigenvalues , **eigenvectors;
    // printf("ENtered Jacobi\n"); 
    Jacobi(DD_T, M, &eigenvalues, &eigenvectors);
    // printf("Exited Jacobi\n"); 
//      cudaEventRecord(start);
// cudaEventRecord(stop);
//     cudaEventSynchronize(stop);
//     cudaEventElapsedTime(&computation_time, start, stop);
//     printf("Time Jacobi: %f\n",computation_time );

    ///////////////////////////////////////////////////////////////////








    // iter<100);//
    // cout << "Broke From Loop EValues: \n";
    // display_matrix(EValues, M, M);
    // cout << "Broke From Loop EVectors: \n";

    // display_matrix(EVectors, M, M);
    // for (int i = 0; i < 4; ++i)
    // {
    //  cerr << EValues[i*M+i]<<" ";
    // }
    // exit(0);
    // 
    // Maintain A structure to arrange eigenvectors
    // in order of sorted eigenvalues
    EigenValIdx * EigenValIdxs = (EigenValIdx *) malloc(sizeof(EigenValIdx)*M);
    // U is NxN
    // SIGMA is NxM but a 1D array of size N rest all are 0 entries
    // V_T is MxM

    // Assuming M >= N
    for (int i = 0; i < M; ++i)
    {
        // SIGMA[0][i] = sqrt(EValues[i*M + i]);
        SIGMA[0][i] = sqrt(eigenvalues[i]);

        
        EigenValIdx egnval;
        // egnval.value = sqrt(EValues[i*M + i]);
        egnval.value = sqrt(eigenvalues[i]);

        egnval.index = i;

        EigenValIdxs[i] = egnval;
        // cerr << SIGMA[0][i] << " ";
    }
    // display_matrix(EVectors, M, M);
    // cout << "Broke From Loop EValues agian: \n";

    // for (int i = 0; i < M; ++i)
    // {
    //  cerr << EValues[i*M+i]<<" ";
    // }
    // cout << "Broke From Loop SIGMA agian: \n";

    // for (int i = 0; i < M; ++i)
    // {
    //  cout << endl;
    //  cerr << i<<" " ;
    //  cerr << SIGMA[0][i]<<" ";
    // }

    // store singular value in descending order
        // qsort((void*)arr, size, sizeof(arr[0]), comparator); 
    qsort((void *)SIGMA[0],M,sizeof(SIGMA[0][0]),compare);
    qsort((void *)EigenValIdxs, M, sizeof(EigenValIdxs[0]),compare_struct);
    // cout << "Broke From Loop SIGMA soreted: \n";

    // for (int i = 0; i < M; ++i)
    // {
    //  cout << endl;
    //  cerr << i<<" " ;
    //  cerr << SIGMA[0][i]<<" ";
    // }
    // SIGMA_INV is MxN
    


    // double* SIGMA_INV = (double *) calloc( M * N, sizeof(double))  ;
    // for (int i = 0; i < M; ++i)
    // {
    //     SIGMA_INV[i*N + i] = 1.0/SIGMA[0][i];
    // }




    // double** mat_transpose(double** A, int Am, int An) {
    double **SIGMA_INV;
    SIGMA_INV = (double**)malloc(__SIZEOF_POINTER__*M);
    for (int i=0; i<M; i++)
        SIGMA_INV[i] = (double*)malloc(__SIZEOF_DOUBLE__*N);

    for (int i=0; i<M; i++){
        for (int j=0; j<N; j++){
            SIGMA_INV[i][j] = 0;
            if (i==j){
                 SIGMA_INV[i][j] = 1.0/SIGMA[0][i];
            }
        }
    }

    // return B;













    // cout << "SIGMA_INV::::::::::::::::::::::::::::::::::::";
    // display_matrix(SIGMA_INV, M, N);


    double* EVectors_T = (double *) malloc(sizeof(double) * M * M)  ;
    // transpose(EVectors, M, M, &EVectors_T);
    
    // Fill the eigenvectors in V_T row-wise

    double ** V_T_temp ;
    V_T_temp = (double**)malloc(__SIZEOF_POINTER__*M);
    for (int i=0; i<M; i++)
        V_T_temp[i] = (double*)malloc(__SIZEOF_DOUBLE__*M);

    // #pragma omp parallel for collapse(2)
    for (int i = 0; i < M; ++i)
    {
        // int idx = EigenValIdxs[i].index;
        // REMOVED ADDED OMP
        // assign(V_T, idx, EVectors_T, M);
        for (int j = 0; j < M; ++j)
        {
            if (i<M){
                // V_T[0][EigenValIdxs[i].index*M+j] = EVectors[j*M+EigenValIdxs[i].index];
                V_T_temp[EigenValIdxs[i].index][j] = eigenvectors[j][EigenValIdxs[i].index];

            }
            else{
                // V_T[0][i*M + j] = 0;
                V_T_temp[i][j] = 0;

            }
        }
    }
    // (double ** U , int ColNo, double * A , int n){

    
    


    // Fill rows with index >= N of V_T with zeros 
    // for (int i = N; i < M; ++i)
    // {
    //  for (int j = 0; j < M; ++j)
    //  {
    //      V_T[0][i*M + j] = 0;
    //  }
    // }

    // double* V = (double *) malloc(sizeof(double) * M * M) ;
    // transpose(V_T[0], M, M, &V);/

    double ** V = mat_transpose(V_T_temp, M, M);
    

// double* MV = (double *) malloc(sizeof(double) * N * M) ;
// multilply_matrices(D_T, N, M, V, M, M, &MV);
// multilply_matrices(MV, N, M, SIGMA_INV, M, N, U);

    
    double ** MV  = mat_mul(D_T, N, M, V, M, M); 
    double ** U_temp   = mat_mul(MV, N, M, SIGMA_INV, M, N);


    copytheOutput(U_temp,U,N,N);
    copytheOutput(V_T_temp,V_T,M,M);

    // cerr << "This is Matrix U  inside SVD_T :"<<endl;
    // display_matrix(U[0], N, N);
    // cerr << "These are the singularvalues"<<endl;
    // for (int i = 0; i < N; ++i){cerr << SIGMA[0][i]<<" ";}cerr << endl;
    // cerr << "This is Matrix V :"<<endl;
    // display_matrix(V, M, M);
    // cout << "Product Check \n";
    // check_product(U[0], SIGMA[0], V_T[0], D_T, M, N);

    return;
}








__host__ void SVD(int M, int N, double* D, double** U, double** SIGMA, double** V_T){
    
    double* D_T = (double *) malloc(sizeof(double) * N * M);
    transpose(D,M,N,&D_T);
    double* U_T = (double *) malloc(sizeof(double) * N * N);
    double* V = (double *) malloc(sizeof(double) * M * M);
    double* SIGMA_HAT = (double *) malloc(sizeof(double) * N);


    // printf("ENtered SVD_T \n");

    SVD_T(N, M, D_T, &V, &SIGMA_HAT, &U_T);
    // printf("Exited SVD_T\n");
    // check_product(V, SIGMA_HAT, U_T, D, N, M);
    // cout << "This is V output: \n";
    // display_matrix(V,M,M);
    // display_matrix(U_T,)
    // display_matrix()
    // copyMatrix2(V_T,V_T,N,N,N);
    

////////////////////////////////////////////////////////
    
    U[0] = (double*) malloc(sizeof(double) * N*N);
    SIGMA[0] = (double*) malloc(sizeof(double) * N);
    V_T[0] = (double*) malloc(sizeof(double) * M*M);
    // transpose(U_T,N,N,U);
    // transpose(V,M,M,V_T);
    double ** V_2d = make2D(V,M,M);
    double ** U_T_2d = make2D(U_T,N,N);


    double ** U_temp = mat_transpose(U_T_2d,N,N);
    double ** V_T_temp    = mat_transpose(V_2d,M,M);




    // Alloc U and V_T

    copytheOutput(U_temp,U,N,N);
    copytheOutput(V_T_temp,V_T,M,M);
    

    // copyMatrix2(U,V,M,M,M);
    for (int i = 0; i < N; ++i)
    {
        SIGMA[0][i] = SIGMA_HAT[i];
    }
    

////////////////////////////////////////////////////////



    // cerr << "This is Matrix U :"<<endl;
    // display_matrix(U[0], N, N);
    // cerr << "These are the singularvalues"<<endl;
    // for (int i = 0; i < N; ++i){cerr << SIGMA[0][i]<<" ";}cerr << endl;
    // cerr << "This is Matrix V_T :"<<endl;
    // display_matrix(V_T[0], M, M);
    // check_product(U[0], SIGMA[0], V_T[0], D_T, M, N);
    // check_product(U, SIGMA_HAT, U_T, D, N, M);

    return;

    // int BLKSIZ = 16;
    // // Display Input
    // cerr << "Matrix D : ";
    // cerr << M<<"x"<< N ;
    // // display_matrix(D,M,N);
    // // exit(0);
    
    // // Calculate transpose of D
    // double* D_T = (double *) malloc(sizeof(double) * N * M);
    // transpose(D,M,N,&D_T);
    // cerr << "Calculating SVD of D_T : ";
    // // display_matrix(D_T, N, M);
    
    // // Calculate D multiplied by D_T
    // // as we have to calculate SVD  of D of D_T
    // double* DD_T = (double *) malloc(sizeof(double) * M * M) ;
    // // multilply_matrices(D, M, N, D_T, N, M, &DD_T);
    // // display_matrix(DD_T, M, M);
    
    // // Memory allocation for matrices
    // double* Q = (double *) malloc(sizeof(double) * M * M) ;
    // double* R = (double *) malloc(sizeof(double) * M * M) ;
    // double* Q_T = (double *) malloc(sizeof(double) * M * M) ;
    
    // // D0
    // double* EValues = (double *) malloc(sizeof(double) * M * M)  ;
    // double* EValues_old = (double *) malloc(sizeof(double) * M * M)  ;

    // // E0
    // double* EVectors = (double *) malloc(sizeof(double) * M * M) ;
    // double* EVectors_old = (double *) malloc(sizeof(double) * M * M)  ;
    // double* Check = (double *) malloc(sizeof(double) * M * M) ;
    // double* Check_inp = (double *) malloc(sizeof(double) * M * M) ;

    // // Initialize 
    // // D0 = D
    // // copyMatrix(&EValues , DD_T, M, M);
    // multilply_matrices_diagonal(D, M, N, D_T, N, M, &EValues,&EVectors);

    // // E0 = I
    // // makeDiagonal(&EVectors, M, M,1.0);
    // // cerr << "This is DD_T: "<<endl;
    // // display_matrix(DD_T, M, M);
    // int iter = 0;
    // do{
        
    //  // Store old values for convergence check
    //  // copyMatrix(&EValues_old, EValues, M, M);
    //  copyMatrixTogether(&EVectors_old, EVectors, M, M,&EValues_old,EValues);

    //  // Q-R Decomposition
    //  // We get Qi,Ri
    //  ModifiedGramSchmidt(EValues, M, M, &Q , &R);
    //  // GramSchmidt(EValues, M, M, &Q , &R);
        
    //  // cerr << "Q"<<endl;
    //  // display_matrix(Q,M,M);
    //  // cerr << "R"<<endl;
    //  // display_matrix(R,M,M);
    //  // cerr << "Check "<<iter<<endl;
    //  //  multilply_matrices(Q, M, M, R, M, M,&Check);
    //  //  display_matrix(Check,M,M);
    //  //  cerr << "Next_INPUT:" <<endl;
    //  //  multilply_matrices(R, M, M, Q, M, M,&Check_inp);
    //  //  display_matrix(Check_inp,M,M);
    //  // exit(0);
    //  // Multiply Ri,Qi
    //  // Di+1 = Ri * Qi
    //  // transpose(Q,M,M,&Q_T);
    //  multilply_matrices_together(R, M, M, Q, M, M,&EValues,EVectors_old,Q,&EVectors);

    //  // Multiply Ei,Qi
    //  // Ei+1 = Ei*Qi
    //  // multilply_matrices(EVectors_old, M, M, Q, M, M,&EVectors);
        

        
    //  iter+= 1;

    // }while(!(convergence(EValues, EValues_old, N, M)) || iter >100);// || convergence(EVectors, EVectors_old, M, M)) );
    // // iter<100);//
    // // display_matrix(EValues, M, M);
    // // display_matrix(EVectors, M, M);
    // for (int i = 0; i < 4; ++i)
    // {
    //  cerr << EValues[i*M+i]<<" ";
    // }
    // // exit(0);
    // // 
    // // Maintain A structure to arrange eigenvectors
    // // in order of sorted eigenvalues
    // EigenValIdx * EigenValIdxs = (EigenValIdx *) malloc(sizeof(EigenValIdx)*N);
    // // U is NxN
    // // SIGMA is NxM but a 1D array of size N rest all are 0 entries
    // // V_T is MxM

    // // Assuming M >= N
    // for (int i = 0; i < N; ++i)
    // {
    //  SIGMA[0][i] = sqrt(EValues[i*M + i]);
        
    //  EigenValIdx egnval;
    //  egnval.value = sqrt(EValues[i*M + i]);
    //  egnval.index = i;

    //  EigenValIdxs[i] = egnval;
    //  // cerr << SIGMA[0][i] << " ";
    // }


    // // store singular value in descending order
    //     // qsort((void*)arr, size, sizeof(arr[0]), comparator); 
    // qsort((void *)SIGMA[0],N,sizeof(SIGMA[0][0]),compare);
    // qsort((void *)EigenValIdxs, N, sizeof(EigenValIdxs[0]),compare_struct);

    // // SIGMA_INV is MxN
    // double* SIGMA_INV = (double *) calloc( M * N, sizeof(double))  ;
    // for (int i = 0; i < N; ++i)
    // {
    //  SIGMA_INV[i*N + i] = 1.0/SIGMA[0][i];
    // }
    // // display_matrix(SIGMA_INV, M, N);


    // double* EVectors_T = (double *) malloc(sizeof(double) * M * M)  ;
    // // transpose(EVectors, M, M, &EVectors_T);
    
    // // Fill the eigenvectors in V_T row-wise

    // #pragma omp parallel for collapse(2)
    // for (int i = 0; i < N; ++i)
    // {
    //  // int idx = EigenValIdxs[i].index;
    //  // REMOVED ADDED OMP
    //  // assign(V_T, idx, EVectors_T, M);
    //  for (int j = 0; j < M; ++j)
    //  {
    //      if (i<N){
    //          V_T[0][EigenValIdxs[i].index*M+j] = EVectors[j*M+EigenValIdxs[i].index];
    //      }
    //      else{
    //          V_T[0][i*M + j] = 0;
    //      }
    //  }
    // }
    // // (double ** U , int ColNo, double * A , int n){

    
    


    // // Fill rows with index >= N of V_T with zeros 
    // // for (int i = N; i < M; ++i)
    // // {
    // //   for (int j = 0; j < M; ++j)
    // //   {
    // //       V_T[0][i*M + j] = 0;
    // //   }
    // // }

    // double* V = (double *) malloc(sizeof(double) * M * M) ;
    // transpose(V_T[0], M, M, &V);


    

    // double* MV = (double *) malloc(sizeof(double) * N * M) ;
    // multilply_matrices(D_T, N, M, V, M, M, &MV);
    // multilply_matrices(MV, N, M, SIGMA_INV, M, N, U);

    // // cerr << "This is Matrix U :"<<endl;
    // // display_matrix(U[0], N, N);
    // // cerr << "These are the singularvalues"<<endl;
    // for (int i = 0; i < N; ++i){cerr << SIGMA[0][i]<<" ";}cerr << endl;
    // // cerr << "This is Matrix V :"<<endl;
    // // display_matrix(V, M, M);

    // // check_product(U[0], SIGMA[0], V_T[0], D_T, M, N);

    // return;
}



__host__ void PCA(int retention, int M, int N, double* D, double* U, double* SIGMA, double** D_HAT, int *K){

    
    double ** D_2d = make2D(D,M,N);

    double EgnValSum = 0.0;
    // #pragma omp parallel for reduction(+:EgnValSum)
    for (int i = 0; i < N; ++i)
    {
        EgnValSum += SIGMA[i]*SIGMA[i];
    }
    // cerr << "EigenValSum : "<<EgnValSum<<endl;

    K[0] = 0;
    double componentSum = 0.0;
    
    
    for (int i = 0; i < N; ++i)
    {   
        
        componentSum = (componentSum + SIGMA[K[0]]*SIGMA[K[0]]);
        K[0] += 1;
        // cerr <<"componentSum : "<<i<<" " <<componentSum << endl;
        if(componentSum/EgnValSum * 100 >= retention){
            break;
        }
        
    }

    int k = K[0];
    double * W = (double *) malloc(sizeof(double)* N * k);
    copyMatrix2(&W, U, N, k, N);

    double ** W_2d = make2D(W,N,k);

    // cerr << "This is Matrix W : "<<endl  ;
    // display_matrix(W, N, k);

    D_HAT[0] = (double *) malloc(sizeof(double)* M * k);

    // multilply_matrices(D, M, N, W, N, K[0], D_HAT);
    double **D_HAT_temp  =   mat_mul(D_2d, M, N, W_2d, N, K[0] );


    // printf("M %d k %d\n",M,k);
    copytheOutput(D_HAT_temp,D_HAT,M,k);




    // cerr << "This is Matrix D_HAT : "<<endl ;
    // display_matrix(D_HAT[0], M, k);
    // cout <<"K is : "<< k<<endl;
    return;

}








// /*
//  *****************************************************
//      TODO -- You must implement this function
//  *****************************************************
// */
__host__ void SVD_and_PCA (int M, 
        int N, 
        double* D, 
        double** U, 
        double** SIGMA, 
        double** V_T, 
        int* SIGMAm,
        int* SIGMAn, 
        double** D_HAT, 
        int *K,
        int retention) {
    // write your code here



    //  copytheOutput and make2D check 
//     double * MO = (double *) malloc(10*sizeof(double));
//     for (int i = 0; i < 10; ++i)
//     {
//        MO[i] = i;
//         if(i%5==0){
//             printf("\n");
//         }
//        printf("%lf ",MO[i] );

//     }
//     printf("\n");


//     double ** M_2d = make2D(MO,2,5);


//     for (int i = 0; i < 2; ++i)
//     {
//         for (int j =0;j<5;j++){
//             printf("%lf ",M_2d[i][j] );
//         }
//         printf("\n");
//     }

//     double * MO2 = (double *) malloc(10*sizeof(double));

//     copytheOutput(M_2d,&MO2,2,5);
// for (int i = 0; i < 10; ++i)
//     {
//        // MO[i] = i;
//         if(i%5==0){
//             printf("\n");
//         }
//        printf("%lf ",MO2[i] );

//     }
//     printf("\n");

//     exit(0);

//     double * A = (double *)malloc(sizeof(double)*4*4 );
//     double * B = (double *)malloc(sizeof(double)*4*4 );
//     double * C = (double *)malloc(sizeof(double)*4*4 );

//     for (int i = 0; i < 4; ++i)
//     {
//         for (int j = 0; j < 4; ++j)
//         {
//             A[4*i+j] = i;
//                 B[4*i+j] = i;
//             // if (i==j){
//             //     A[4*i+j] = 1;
//             //     B[4*i+j] = 1;

//             // }
//             // else{
//             //     A[4*i+j] = 0;
//             //     B[4*i+j] = 0;
//             // }
//         }
//     }


//     Matrix Am , Bm, Cm;
//     Am.width = 4;
//     Am.height = 4;
//     Am.stride = 4;
// Am.elements = A;
//     Bm.width = 4;
//     Bm.height = 4;
//     Bm.stride = 4;
// Bm.elements = B;

// Cm.width = 4;
//     Cm.height = 4;
//     Cm.stride = 4;
// Cm.elements = C;


//     MatMulCu(Am,Bm,Cm);

//    printf("\n");
//     for (int i = 0; i < 4; ++i)
//     {
//         for (int j = 0; j < 4; ++j)
//         {
//             printf("%lf ",Am.elements[i*4+j] );
//         }
//         printf("\n");
//     }printf("\n");
//     for (int i = 0; i < 4; ++i)
//     {
//         for (int j = 0; j < 4; ++j)
//         {
//             printf("%lf ",Bm.elements[i*4+j] );
//         }
//         printf("\n");
//     }
//     printf("\n");
//     for (int i = 0; i < 4; ++i)
//     {
//         for (int j = 0; j < 4; ++j)
//         {
//             printf("%lf ",Cm.elements[i*4+j] );
//         }
//         printf("\n");
//     }

// return;
//     printf("ENtered SVD M %d N %d \n",M,N);


//      cudaEvent_t start, stop;

// cudaEventCreate(&start);
//     cudaEventCreate(&stop);

//     float computation_time;

//     double * eigenvalues , **eigenvectors;
    // printf("NEW\n"); 
    SVD(M, N, D,  U, SIGMA,  V_T);

//     printf("Exited SVD\n"); 
//      cudaEventRecord(start);
// cudaEventRecord(stop);
//     cudaEventSynchronize(stop);
//     cudaEventElapsedTime(&computation_time, start, stop);
//     printf("Time SVD: %f\n",computation_time );

//     printf("ENtered PCA\n");

// cudaEventCreate(&start);
//     cudaEventCreate(&stop);
    PCA(retention, M, N, D, U[0], SIGMA[0], D_HAT, K);
//     printf("Exited PCA\n"); 
//     cudaEventRecord(start);
// cudaEventRecord(stop);
//     cudaEventSynchronize(stop);
//     cudaEventElapsedTime(&computation_time, start, stop);
//     printf("Time PCA: %f\n",computation_time );
    
}

















