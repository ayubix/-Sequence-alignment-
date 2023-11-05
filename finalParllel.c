#include <mpi.h>
#include "finalParallel.h"
#include "cuda.h"

char PrimaryGroup[PRIMARY_GROUP_LENGTH][MAX_PRIMARY] = {"NDEQ", "NEQK", "STA", "MILV", "QHRK", "NHQK", "FYW", "HY", "MILF"}; // defined the two word groups
char SecondaryGroup[SECONDARY_GROUP_LENGTH][MAX_SECONDARY] = {"SAG", "ATV", "CSA", "SGND", "STPA", "STNK", "NEQHRK", "NDEQHK", "SNDEQK", "HFY", "FVLIM"};





char **readInput(int *weights, char **seq1, int *numberStrings)
{ // read the given input file and store all the data for the next stages

    char seq1Temp[MAX_SEQ1]; // init varibales
    char **seq2Arr;
    FILE *f = fopen(INPUT_FILE, "r"); // open file for reading
    if (f == NULL)                    // check if open succeed
    {
        perror("Error");
        printf("Could not open file %s", INPUT_FILE);
    }
    fscanf(f, "%d %d %d %d", weights, weights + 1, weights + 2, weights + 3); // get the weights

    fscanf(f, "%s", seq1Temp); // save the seq1 string , and allocate memory for it
    *seq1 = (char *)malloc(sizeof(char) * strlen(seq1Temp) + 1);
    if (*seq1 == NULL)
    {
        printf("error with seq1\n");
    }
    strcpy(*seq1, seq1Temp); // copying the temporary string to permanent
    fscanf(f, "%d", numberStrings); // get the number of strings in the seq2 part

    seq2Arr = (char **)malloc(sizeof(char *) * *numberStrings); // allocation for strings array
    if (seq2Arr == NULL)
    { // check if  the allocation worked
        printf("error with allocation sequences \n");
    }

    for (int i = 0; i < *numberStrings; i++)
    {
        seq2Arr[i] = (char *)malloc(MAX_SEQ2 * sizeof(char)); // allocate space for each string
        if (seq2Arr[i] == NULL)
        { // check if  the allocation worked
            printf("error with allocation string number %d \n", i + 1);
        }
        fscanf(f, "%s", seq2Arr[i]); // get the string
    }
    fclose(f);      // closing the file
    return seq2Arr; // return the strings array to main for future usage
}





void printInputData(int *weights, int numOfSequenes, char *seq1, char **seq2)
{ // print the file input data for debugging purposes

    printf("The weights are %d %d %d %d\n", weights[0], weights[1], weights[2], weights[3]);
    printf("The number of string for comparison with s1 are %d \n", numOfSequenes);
    int len = strlen(seq1);
    printf("Seq1 is ");
    for (int i = 0; i < len; i++)
    {
        printf("%c", seq1[i]);
    }
    printf("\n");

    for (int i = 0; i < numOfSequenes; i++)
    {
        printf("Seq2 number %d is  %s \n", i + 1, seq2[i]);
        
    }
}




int alignmentScore(int *weights, int *symbolCount)
{ // calculates the alignment score of two strings according to given equation
    int score, positiveScore, negetiveScore;
    positiveScore = weights[0] * symbolCount[0];
    negetiveScore = weights[1] * symbolCount[1] + weights[2] * symbolCount[2] + weights[3] * symbolCount[3];
    score = positiveScore - negetiveScore;
    return score; // return the ultimate score of the to strings
}




void fillMatrix(char matrix[MATRIXSIZE][MATRIXSIZE])
{ // Filling the matrix with the relevant symbol for each pair of letters according the instructions
  // and assign it in the matrix[charA][charB] cell. This allow us to calculate the relation between
  // the letters only one time and be able to compare between all the strings in high efficiency

    
    for (char i = 'A'; i <= 'Z'; i++) // running on the matrix for assigning all the symbols
    {

        for (char j = 'A'; j <= 'Z'; j++)
        {
            matrix[i - 'A'][j - 'A'] = sortPairToGroup(i, j); // the result from the comparison function between chars i and j
        }
    }
}




char sortPairToGroup(char a, char b)
{ // returns to the filling matrix the appropriate sign for the pair(a,b)

    if (a == b)
        return '$'; // return $ in case the chars are the same
    else if (belongToPrimary(a, b) == 1)
        return '%'; // return % if the chars belong to primary group

    else if (belongToSecondary(a, b) == 1)
        return '#'; // return # if the chars belong to secodary group
    return ' ';     // return space when the chars do not meet any of the previous conditions
}




int belongToPrimary(char a, char b)
{ // Return 1 if the letters belong to the primary group , else return 0

    for (int i = 0; i < PRIMARY_GROUP_LENGTH; i++) // running on the primary group words and ask whether a and b part of it
    {
        if (strchr(PrimaryGroup[i], a) != NULL && strchr(PrimaryGroup[i], b) != NULL) // if both of the letters exist in the current word return 1
                                                                                      //,else move forward to the next word
            return 1;
    }

    return 0; // In case none of the group words contains both of the letters
}




int belongToSecondary(char a, char b)
{ // Return 1 if the letters belong to the secondary group , else return 0.(Same logic as primary function
    for (int i = 0; i < SECONDARY_GROUP_LENGTH; i++)
    {
        if (strchr(SecondaryGroup[i], a) != NULL && strchr(SecondaryGroup[i], b) != NULL)
            return 1;
    }

    return 0;
}




void compareStrings(char *seq1, char **seq2, int seq2Size, char matrix[MATRIXSIZE][MATRIXSIZE], int *weights, int id, int procNumber)
{  // The function send the input strings for the comparison function and print the maximum score argument for each of them

    int maxArraySize = 3 * seq2Size; // The size of the array that store the local maximum parmaters of processes,local maximum for each sting in seq2
    int *masterMaxArray; // An array that the master will use for collecting all the local maximum from the other slaves,and find the global maximum

    if (id == 0) // Only the master allocate this array for find and print the best score parameters of every string in seq2
    {
        masterMaxArray = (int *)malloc(sizeof(int) * maxArraySize * procNumber);
        if (masterMaxArray == NULL)
        {
            printf("allocation error with maxScore array\n"); // In case the allocation failed
            return;
        }
    }

    int *maxScoreArr; // array that store the local maximum parmaters of processes,local maximum for each sting in seq2

    maxScoreArr = (int *)malloc(sizeof(int) * maxArraySize); // each slave allocate this array
    if (maxScoreArr == NULL)
    {
        printf("allocation error with maxScore array\n"); // In case the allocation failed
        return;
    }

    for (int i = 0; i < seq2Size; i++) // running on the seq2 array , and compare each string with seq1
    {
        maxScoreArr[i * 3] = INT_MIN; // prevent the situation that all the score are negetive and the initial parameter value is 0
        findMaxScoreParameters(seq1, seq2[i], matrix, weights, maxScoreArr + i * 3, id, procNumber); // execute the function that update the max array with the
                                                                                                        // relevnat arguments
    }

    MPI_Gather(maxScoreArr, maxArraySize, MPI_INT, masterMaxArray, maxArraySize, MPI_INT, 0, MPI_COMM_WORLD);// gather all the results from the other slaves 
                                                                                                                // to the master for another check who is the 
                                                                                                                //global maximum of each string
    if (id == 0) //Only master performs this
    { 
        for (int i = 0; i < seq2Size; i++)
        {
            for (int j = 0; j < procNumber; j++)
            {
                if (masterMaxArray[j * maxArraySize + i * 3] > maxScoreArr[i * 3])
                {
                                                                                            //All of these logic , in order to figure out what is the max score of 
                                                                                        // the sequences,running on the array size with nested loop for these task
                                                                                    //using the previous array for storing the global max values and print them

                    maxScoreArr[i * 3] = masterMaxArray[j * maxArraySize + i * 3];
                    maxScoreArr[i * 3 + 1] = masterMaxArray[j * maxArraySize + i * 3 + 1];
                    maxScoreArr[i * 3 + 2] = masterMaxArray[j * maxArraySize + i * 3 + 2];
                }
            }
            printf("The best alignment score of seq2 number %d is %d , offset = %d ,hypen = %d\n", i + 1, maxScoreArr[i * 3], maxScoreArr[i * 3 + 1], maxScoreArr[i * 3 + 2]);
            // printing the best score paramater of these sequence
        }
        free(masterMaxArray);           // free all the allocated memory in the function
    }
    free(maxScoreArr);
}




void findMaxScoreParameters(char *seq1, char *str, char matrix[MATRIXSIZE][MATRIXSIZE], int *weights, int *scoreArr, int id, int procNumber)
{ // The function compare between two strings
    // char by char and find the offset and hyphen location that produce
    // the best possible alignment score,save it in the argument array for future use , each process work on different part of the string
    // and using also opem mp threads and cuda functionallity

    int strLen = strlen(str);
    int lenDiff = strlen(seq1) - strLen;
    char result;

    //The logic here is the same as in the sequential version , but in this parallel version , I changed the outer loop from the way it was in the serial problem  , so each process work on 
    // number of chars in the string. The calculation is length of the string divide to the num of processes(it is the offset loop)
    //and later each of the omp threads will be responsible to part of this chunk or the cuda will calculate it instead





    int workPerProcess = strLen / procNumber;
    int processStart, processEnd;			 // initialize the start and end pointers for each process
    processStart = workPerProcess * id;
    processEnd = processStart + workPerProcess;

    if (id == procNumber - 1)
    {
        processEnd = strLen;    // In case we have a remainder
    }
    if (id == 0)
    {
        processStart = 1;    // The loop should start from 1
    }



#pragma omp parallel private(result) // entering a parllel section,declare on result as private
    {
        int threadScoreArr[3]; // new array decalared inside the parallel section , prevent any race condition section issues
        threadScoreArr[0] = INT_MIN; // prevent the situation that all the score are negetive

        int thread_start, thread_end, thread_work, num_threads, tid; // variables for deviding the work between the program threads
        tid = omp_get_thread_num(); // get the number of the current thread for future use
        num_threads = omp_get_num_threads() - 1; // thread 0 doesn't participate with other threads load

        if (tid == 0) // only thread number 0 do this part,half of the process work will be done by cuda ,and this thread will execute it
        {
            thread_start = processStart; 
            thread_end = processStart + processEnd / 2; // initialize the start and end pointers for this process
            findMaxScoreCuda(seq1, str, matrix, weights, threadScoreArr, thread_start, thread_end); // calling to the cuda function
        }



        else // the other threads will be responsible to the second half of the process work , each of them will work on equal part
			// processwork : 2 : num of theres is the calculation.
        {
            thread_work = workPerProcess / num_threads;
            thread_start = processStart + processEnd / 2 + (tid - 1) * thread_work; // initialize the start and end pointers for this process  
            thread_end = thread_start + thread_work;

            if (tid == num_threads - 1)
                thread_end = processEnd; // In case we have a remainder


            for (int i = thread_start; i <= thread_end; i++) // for all the hyphen options,in this version go over the partial string of the process
            {
                for (int j = 0; j <= lenDiff; j++) // for all the offset options
                {
                    int symbolCount[4] = {0};        // Hold for each symbol the number of appearences in comparision ,for calculate score
                    for (int k = 0; k < strLen; k++) // go over str2 and compare to seq1 char by char and calculate score
                    {
                        if (k >= i)
                        {
                            result = matrix[seq1[k + j + 1] - 'A'][str[k] - 'A']; // get the comparison symbol from the filled matrix
                            updateSymbolArray(result, symbolCount); // update the symbol array with the given sign
                        }
                        else
                        {
                            result = matrix[seq1[k + j] - 'A'][str[k] - 'A']; // get the comparison symbol from the filled matrix
                            updateSymbolArray(result, symbolCount); // update the symbol array with the given sign
                        }
                    }
                    int currentScore = alignmentScore(weights, symbolCount); // getting the current aligment score
                    if (currentScore > threadScoreArr[0]) // check if it is greater than the local maximum
                    {
                        threadScoreArr[0] = currentScore; // store the max alignment score of two strings
                        threadScoreArr[1] = j;            // save the offset value
                        threadScoreArr[2] = i;            // save the hypen value
                    }
                }
            }
        }




#pragma omp critical // Here I passed the maximum values to the argument array in the critical section , to avoid any race condition
                    // only one thread can do it in a time
        {
            if (threadScoreArr[0] > scoreArr[0]) // check whether the new score bigger than the current maximum score
            {
                scoreArr[0] = threadScoreArr[0]; // store the max alignment score of two strings
                scoreArr[1] = threadScoreArr[1]; // save the offset value
                scoreArr[2] = threadScoreArr[2]; // save the hypen value
            }
        }
    }
}





void updateSymbolArray(char symbol, int *arr) 
{ 
    // Function update the symbol array ,and add 1 to the proper index
    if (symbol == '$')
        arr[0]++;
    else if (symbol == '%')
        arr[1]++;
    else if (symbol == '#')
        arr[2]++;
    else
        arr[3]++;
}




int main(int argc, char *argv[]) // Main part
{
    int my_rank, num_procs; // initialize variables and all the MPI required arguments
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    double t = MPI_Wtime(); // for calculating program time

    int weights[4]; // define all the data structures for input
    int numOfSequenes;
    char *seq1;
    char **seq2;
    int seq1Len, tempStrLen; // the second variable will represent the length of seq2 strings
    char matSymbol[MATRIXSIZE][MATRIXSIZE];

    if (my_rank == 0)
    { // only the master read the file

        seq2 = readInput(weights, &seq1, &numOfSequenes);   // execute the file reading function
        MPI_Bcast(weights, 4, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&numOfSequenes, 1, MPI_INT, 0, MPI_COMM_WORLD);
        seq1Len = strlen(seq1);
        MPI_Bcast(&seq1Len, 1, MPI_INT, 0, MPI_COMM_WORLD);          // sharing with the other slaves the file data in the proper order
        MPI_Bcast(seq1, seq1Len, MPI_CHAR, 0, MPI_COMM_WORLD);
        for (int i = 0; i < numOfSequenes; i++)
        {
            tempStrLen = strlen(seq2[i]);
            MPI_Bcast(&tempStrLen, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(seq2[i], tempStrLen, MPI_CHAR, 0, MPI_COMM_WORLD);
        }
    }

    // workers section , get the input data from the master and allocate memory for it
    else
    {
        MPI_Bcast(weights, 4, MPI_INT, 0, MPI_COMM_WORLD); // get the weight array values
        MPI_Bcast(&numOfSequenes, 1, MPI_INT, 0, MPI_COMM_WORLD); // save the number of seq2 strings
        MPI_Bcast(&seq1Len, 1, MPI_INT, 0, MPI_COMM_WORLD);
        seq1 = (char *)malloc(sizeof(char) * seq1Len + 1);
        if (seq1 == NULL)
        {
            printf("error with seq1 allocation\n");
            return 0;
        }
        MPI_Bcast(seq1, seq1Len, MPI_CHAR, 0, MPI_COMM_WORLD); // saving seq1 value after allocation enough space for it
        seq1[seq1Len] = '\0'; // add the zero char to end of string

        seq2 = (char **)malloc(sizeof(char *) * numOfSequenes); // allocation the strings array size
        if (seq2 == NULL)
        {
            printf("seq2 allocation failed\n");
            return 0;
        }
        for (int i = 0; i < numOfSequenes; i++)
        {
            MPI_Bcast(&tempStrLen, 1, MPI_INT, 0, MPI_COMM_WORLD);
            seq2[i] = (char *)malloc(sizeof(char) * tempStrLen + 1);
            if (seq2[i] == NULL)
            {
                printf("seq2 string number %d allocation failed\n", i + 1);
                return 0;
            }

            MPI_Bcast(seq2[i], tempStrLen, MPI_CHAR, 0, MPI_COMM_WORLD); // Getting seq2 value after allocation enough space for it
            seq2[i][tempStrLen] = '\0'; // add the zero char to end of string
            
        }
    }

    fillMatrix(matSymbol); // execute the filling matrix function
    compareStrings(seq1, seq2, numOfSequenes, matSymbol, weights, my_rank, num_procs); // execute the comparison function

    free(seq1);
    for (int i = 0; i < numOfSequenes; i++) {    // free all the allocated memory in the program
        free(seq2[i]);
    }
    free(seq2);


    if (my_rank == 0) 
         printf("Program parallel time is %.4f\n", MPI_Wtime() - t); // print the program time
       
    MPI_Finalize();
    return 0;
}
