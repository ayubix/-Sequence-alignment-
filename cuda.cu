#include "cuda.h"
#include <string.h>
#define THREADS_PER_BLOCK 256




// The same function as in the host file , only with proper decalration

__device__ void cudaUpdateSymbolArray(char symbol, int *arr)
{ // Function update the symbol array ,and add 1 to the proper index
    if (symbol == '$')
        arr[0]++;
    else if (symbol == '%')
        arr[1]++;
    else if (symbol == '#')
        arr[2]++;
    else
        arr[3]++;
}




// The same function as in the host file , only with proper decalration

__device__ int cudaAlignmentScore(int *weights, int *symbolCount)
{ // calculates the alignment score of two strings according to given equation
    int score, positiveScore, negetiveScore;
    positiveScore = weights[0] * symbolCount[0];
    negetiveScore = weights[1] * symbolCount[1] + weights[2] * symbolCount[2] + weights[3] * symbolCount[3];
    score = positiveScore - negetiveScore;
    return score; // return the ultimate score of the to strings
}







__global__ void kernelFunc(char *seq1, char *seq2, char *matrix, int *weights, int *scoreArr, int num_offsets, int seq2_len, int num_threads){
    // the kernel function , run on the gpu each thread will calculate aligmnet score and save it in the device score array

    int tid = blockDim.x * blockIdx.x + threadIdx.x; // calculate the thread index in the grid 
    if(tid > num_threads) // make sure we do not exceeded from the grid size , only relevent for the last block 
        return;

    int hyphen,offset; // initialize variables
    char result;

    hyphen = tid /num_offsets;      // calculate the values of the offset and the mutant for current thread
    offset = tid % num_offsets;

    int symbolCount[4] = {0};        // Hold for each symbol the number of appearences in comparision ,for calculate score

    for (int k = 0; k < seq2_len; k++) // go over str2 and compare to seq1 char by char and calculate aligment score(same logic as in the host)
                                        // but only with the proper functions for the device 
    {
        if (k >= hyphen)
        {
            result = matrix[(seq1[k + offset + 1] - 'A')*MATRIXSIZE+ seq2[k] - 'A'];
            cudaUpdateSymbolArray(result, symbolCount);
        }
        else
        {
            result = matrix[(seq1[k + offset] - 'A')*MATRIXSIZE + seq2[k] - 'A'];
            cudaUpdateSymbolArray(result, symbolCount); 
        }
    }
    scoreArr[tid] = cudaAlignmentScore(weights, symbolCount); 
   
}






void findMaxScoreCuda(char *seq1, char *seq2, char matrix[MATRIXSIZE][MATRIXSIZE], int *weights, int *scoreArr, int start_hyp, int end_hyp){
    char* dev_seq1,*dev_seq2, *dev_matrix;
    int * dev_weights,*dev_score,*host_score;
    int seq1_len, seq2_len, offsets;        // initialize variables
    seq1_len = strlen(seq1);
    seq2_len = strlen(seq2);        // calculate the lengths of the two strings

    offsets = seq1_len- seq2_len;   // caculate the difference between the lengths , for the offset size

    int numthreads = (end_hyp-start_hyp)*offsets; // calculate the number of needed threads, for the mission as the size of hyphen * offset
    int threadsPerBlock = 256; // default size for the number of threds in each block
    int numBlocks = numthreads / threadsPerBlock; // calculate how many blocks we need for that task


    if(numthreads % threadsPerBlock != 0) // if we have a remainder,we will use one extra block
        numBlocks++;


    
    cudaMalloc((void**)&dev_seq1,sizeof(char)*(seq1_len+1));
    cudaMalloc((void**)&dev_seq2,sizeof(char)*(seq2_len+1));
    cudaMalloc((void**)&dev_matrix,sizeof(char)*MATRIXSIZE*MATRIXSIZE);        // allocate memory on device space for all this data stractures
    cudaMalloc((void**)&dev_weights,sizeof(int)*4);                             // for the kernel function to calculate the scores
    cudaMalloc((void**)&dev_score,sizeof(int)*(numthreads));
    cudaMemcpy((void*)(dev_seq1),seq1,(seq1_len+1)*sizeof(char),cudaMemcpyHostToDevice);
    cudaMemcpy((void*)(dev_seq2),seq2,(seq2_len+1)*sizeof(char),cudaMemcpyHostToDevice);
    cudaMemcpy((void*)(dev_matrix),matrix,MATRIXSIZE*MATRIXSIZE*sizeof(char),cudaMemcpyHostToDevice);
    cudaMemcpy((void*)(dev_weights),weights,4*sizeof(int),cudaMemcpyHostToDevice);     // copying all the data to the device memory from the host

    host_score = (int*)malloc(sizeof(int)*(numthreads)); // allocate an array on host ,that store the results calcualted on the device

    kernelFunc<<<numBlocks,threadsPerBlock>>>(dev_seq1,dev_seq2,dev_matrix,dev_weights,dev_score,offsets,seq2_len,numthreads); //execute the kernel function
                                                                                                                                // with all the relevant arguments

    cudaDeviceSynchronize(); // makes the cpu wait until all gpu threads will end their tasks
    cudaMemcpy(host_score,dev_score,sizeof(int)*numthreads,cudaMemcpyDeviceToHost);//transfer the data from the device array to the host


    for(int i = start_hyp; i < end_hyp;i++){ // find the maximum score from all the scores that the cuda calculate , and save it in the scoreArr
        for(int j = 0; j< offsets;j++)
        {
            if(host_score[(i-start_hyp)*offsets+j] > scoreArr[0]){
                scoreArr[0] = host_score[(i-start_hyp)*offsets+j];
                scoreArr[1] = j;
                scoreArr[2] = i;
            }
        }
    }

    
    cudaFree(dev_seq1);
    cudaFree(dev_seq2);
    cudaFree(dev_matrix);      //releasing all the allocated space,of device and later of the host
    cudaFree(dev_score);
    cudaFree(dev_weights);
    free(host_score);

}