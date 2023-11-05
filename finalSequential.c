#include <mpi.h>
#include <stdio.h> // includes
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <ctype.h>     
#include <math.h>




                        // Defines and constants
#define MAX_SEQ1 3000
#define MAX_SEQ2 2000
#define INPUT_FILE "inputNew.txt"
#define PRIMARY_GROUP_LENGTH 9
#define MAX_PRIMARY 5
#define SECONDARY_GROUP_LENGTH 11
#define MAX_SECONDARY 7
#define MATRIXSIZE 26




                    // Globals

char PrimaryGroup[PRIMARY_GROUP_LENGTH][MAX_PRIMARY] = { "NDEQ","NEQK","STA","MILV","QHRK","NHQK","FYW","HY","MILF"}; // defined the two word groups
char SecondaryGroup[SECONDARY_GROUP_LENGTH][MAX_SECONDARY] = {"SAG","ATV","CSA","SGND","STPA","STNK","NEQHRK","NDEQHK","SNDEQK","HFY","FVLIM"};




                    // functions declaration  

char** readInput(int *weights,char **seq1,int *numberStrings) ;
int alignmentScore(int *weights,int *symbolCount);
void printInputData(int * weights,int numOfSequenes,char* seq1,char** seq2);
void fillMatrix(char matrix[MATRIXSIZE][MATRIXSIZE]);
char sortPairToGroup(char a, char b);
int belongToPrimary(char a,char b);
int belongToSecondary(char a,char b);
void compareStrings(char *seq1,char **seq2,int seq2Size,char matrix[MATRIXSIZE][MATRIXSIZE],int *weights);
void findMaxScoreParameters(char *seq1,char *str,char matrix[MATRIXSIZE][MATRIXSIZE],int *weights,int *scoreArr);
void updateSymbolArray(char symbol,int *arr);






char** readInput(int *weights,char **seq1,int *numberStrings) {  // read the given input file and store all the data for the next stages

    char seq1Temp[MAX_SEQ1]; // init varibales
    char **seq2Arr;
    FILE *f = fopen(INPUT_FILE, "r"); // open file for reading
		if (f == NULL) // check if open succeed
		{
			perror("Error");
			printf("Could not open file %s", INPUT_FILE);
			
		}

		fscanf(f, "%d %d %d %d", weights, weights+1, weights+2, weights+3); // get the weights 
        

		fscanf(f, "%s", seq1Temp); // save the seq1 string , and allocate memory for it
        *seq1 = (char*)malloc(sizeof(char) * strlen(seq1Temp) +1);
        if(*seq1 == NULL) {
            printf("error with seq1\n");
            
        }
        strcpy(*seq1,seq1Temp); //copying the temporary string to permanent


		fscanf(f, "%d", numberStrings); //get the number of strings in the seq2 part
		seq2Arr = (char **)malloc(sizeof(char *) * *numberStrings); // allocation for strings array
        if(seq2Arr == NULL) { // check if  the allocation worked
        printf("error with allocation sequences \n");
        }
		for (int i = 0; i < *numberStrings; i++)
		{
			seq2Arr[i] = (char *)malloc(MAX_SEQ2 * sizeof(char)); // allocate space for each string 
            if(seq2Arr[i] == NULL) { // check if  the allocation worked
                printf("error with allocation string number %d \n",i+1);
            }
			fscanf(f, "%s", seq2Arr[i]); // get the string
		}

		fclose(f); //closing the file
        return seq2Arr; //return the strings array to main for future usage
    
}






int alignmentScore(int *weights,int *symbolCount) { // calculates the alignment score of two strings according to given equation
    int score, positiveScore, negetiveScore;
    positiveScore = weights[0] * symbolCount[0];
    negetiveScore = weights[1] * symbolCount[1] + weights[2] * symbolCount[2] + weights[3] * symbolCount[3];
    score = positiveScore - negetiveScore;
    return score; // return the ultimate score of the to strings
}



void printInputData(int * weights,int numOfSequenes,char* seq1,char** seq2) { // print the file input data for debugging purposes

    printf("The weights are %d %d %d %d\n",weights[0],weights[1],weights[2],weights[3]);
    printf("The number of string for comparison with s1 are %d \n" ,numOfSequenes);
    int len = strlen(seq1);
    printf("Seq1 is ");
    for (int i = 0; i <len ; i++)
    {
        printf("%c",seq1[i]);
    }
    printf("\n");
    
     for (int i = 0; i <numOfSequenes ; i++) {
            printf("Seq2 number %d is  %s \n",i+1,seq2[i]);
     }

}




void fillMatrix(char matrix[MATRIXSIZE][MATRIXSIZE]) { // Filling the matrix with the relevant symbol for each pair of letters according the instructions
                                                        // and assign it in the matrix[charA][charB] cell. This allow us to calculate the relation between
                                                        // the letters only one time and be able to compare between all the strings in high efficiency

    
    for (char i = 'A'; i <= 'Z'; i++) // running on the matrix for assigning all the symbols
    {
            for(char j = 'A'; j <='Z'; j++) 
            {
                matrix[i-'A'][j-'A'] = sortPairToGroup(i,j); // the result from the comparison function between chars i and j
            }
            
    }

 }





char sortPairToGroup(char a, char b) { // returns to the filling matrix the appropriate sign for the pair(a,b) 

    if(a == b) 
        return '$'; // return $ in case the chars are the same
    else if(belongToPrimary(a,b) == 1)
        return '%'; // return % if the chars belong to primary group

    else if(belongToSecondary(a,b) == 1)
        return '#'; // return # if the chars belong to secodary group
    return ' '; // return space when the chars do not meet any of the previous conditions


}





int belongToPrimary(char a,char b) { // Return 1 if the letters belong to the primary group , else return 0
   
    for (int i = 0; i < PRIMARY_GROUP_LENGTH; i++) // running on the primary group words and ask whether a and b part of it
    {
        if(strchr(PrimaryGroup[i],a) != NULL && strchr(PrimaryGroup[i],b) != NULL) // if both of the letters exist in the current word return 1
                                                                                    //,else move forward to the next word
            return 1;
            
       }
       
return 0; // In case none of the group words contains both of the letters

}




int belongToSecondary(char a,char b) { // Return 1 if the letters belong to the secondary group , else return 0.(Same logic as primary function
                                        // hence making one function for both of them and sending the group details as arguments considered
                                        // but passing the static global matrix with different sizes turned out as problem so I stayed with these version

     for (int i = 0; i < SECONDARY_GROUP_LENGTH; i++) 
    {
        if(strchr(SecondaryGroup[i],a) != NULL && strchr(SecondaryGroup[i],b) != NULL) 
            return 1;
            
       }
       
return 0; 

}





void compareStrings(char *seq1,char **seq2,int seq2Size,char matrix[MATRIXSIZE][MATRIXSIZE],int *weights) { // The function send the input strings 
// for the comparison function and print the maximum score argument for each of them

     
    for (int i = 0; i < seq2Size; i++) // running on the seq2 array , and working on each string for single iteration
    {
        int maxScoreArr[3] = {0}; //Store the max score , and the offset and hypen values of these score
        findMaxScoreParameters(seq1,seq2[i],matrix,weights,maxScoreArr);
        printf("The best alignment score of seq2 number %d is %d , offset = %d ,hypen = %d\n",i+1,maxScoreArr[0],maxScoreArr[1],maxScoreArr[2]);
        
    }
    
}




void findMaxScoreParameters(char *seq1,char *str,char matrix[MATRIXSIZE][MATRIXSIZE],int *weights,int *scoreArr) { // The function compare between two strings   
                                                             // char by char and find the offset and hyphen location that produce
                                                            // the best possible alignment score,save it in the argument array for future use

   
   
    int strLen = strlen(str);
    int lenDiff = strlen(seq1) - strLen;
    char result;
    

                                    // In order to check which combinations of hyphen and offset
                                    // will yield the best score ,we have to check all the options
                                    // so it will be (H*F) hyphen size can be between 1 to length of str
                                    // and the size of offset between 0 to the difference of str1 to str2
                                    // therefore the next loops assign in this way

    for (int i = 1; i <= strLen ; i++) // for all the hyphen options
    {
        for (int j = 0; j <= lenDiff; j++) // for all the offset options
        {
             int symbolCount[4] = {0}; // Hold for each symbol the number of appearences in comparision ,for calculate score
            for (int k = 0; k < strLen; k++) // go over str2 and compare to seq1 char by char and calculate score
            {
                if(k>=i) {
                     result = matrix[seq1[k+j+1]-'A'][str[k]-'A']; // get the comparison symbol from the filled matrix
                     updateSymbolArray(result,symbolCount); // update the symbol array with the given sign
                     
                }
                else {
                    result = matrix[seq1[k+j]-'A'][str[k]-'A']; // get the comparison symbol from the filled matrix
                    updateSymbolArray(result,symbolCount); // update the symbol array with the given sign
                   
                    
                }
                 
            }
             int currentScore = alignmentScore(weights,symbolCount); // calculate the score of these combination
              
            
             if(currentScore > scoreArr[0]) { // check whether the new score bigger than the current maximum score
                scoreArr[0] = currentScore; // store the max alignment score of two strings
                scoreArr[1] = j; // save the offset value 
                scoreArr[2] = i; //save the hypen value
             }
            

        }

    
    }
    

}








void updateSymbolArray(char symbol,int *arr){ // Function updates the symbol array ,and add 1 to the proper index for future score calculation
    if(symbol == '$')
        arr[0]++;
    else if(symbol == '%')
        arr[1]++;
     else if(symbol == '#')
        arr[2]++;
    else
        arr[3]++;

}








int main(int argc, char *argv[])    // Main part
{
    MPI_Init(&argc, &argv);
    int weights[4];
    int numOfSequenes;
    char* seq1;     // init varibales
    char** seq2;
    char matSymbols[MATRIXSIZE][MATRIXSIZE]; // store in cell mat[i][j] the comparison result between the chars i,j(can be % / $ / # / spcae)

    seq2 = readInput(weights,&seq1,&numOfSequenes); // execute the file reading function
    double t = MPI_Wtime();  // for calculating program time
    fillMatrix(matSymbols); // execute the filling matrix function
    compareStrings(seq1,seq2,numOfSequenes,matSymbols,weights); // comapre all seq2 strings with seq1 and print the best score 
    printf("Program sequential time is %.4f\n",MPI_Wtime()-t); // print the time of the program

    free(seq1);
    for (int i = 0; i < numOfSequenes; i++) {   // free all the allocated memory in the program
        free(seq2[i]);
    }
    free(seq2);
    
    
    MPI_Finalize(); // close mpi protocol
    return 0;
}

