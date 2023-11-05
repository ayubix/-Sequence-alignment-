
#include <stdio.h> // includes
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <omp.h>
#include <limits.h>



#define MAX_SEQ1 3000 // Defines and constants
#define MAX_SEQ2 2000
#define INPUT_FILE "inputNew.txt"
#define PRIMARY_GROUP_LENGTH 9
#define MAX_PRIMARY 5
#define SECONDARY_GROUP_LENGTH 11
#define MAX_SECONDARY 7
#define MATRIXSIZE 26




char **readInput(int *weights, char **seq1, int *numberStrings); // functions declaration  
void printInputData(int *weights, int numOfSequenes, char *seq1, char **seq2);
void fillMatrix(char matrix[MATRIXSIZE][MATRIXSIZE]);
char sortPairToGroup(char a, char b);
char sortPairToGroup(char a, char b);
int belongToPrimary(char a, char b);
int belongToSecondary(char a, char b);
void compareStrings(char *seq1, char **seq2, int seq2Size, char matrix[MATRIXSIZE][MATRIXSIZE], int *weights, int id, int procNumber);
void findMaxScoreParameters(char *seq1, char *str, char matrix[MATRIXSIZE][MATRIXSIZE], int *weights, int *scoreArr, int id, int procNumber);
int alignmentScore(int *weights, int *symbolCount);
void updateSymbolArray(char symbol, int *arr);
