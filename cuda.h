#include "finalParallel.h"  //include
#include <stdlib.h>
#include <string.h>

                        // header for the device function
void findMaxScoreCuda(char *seq1, char *seq2, char matrix[MATRIXSIZE][MATRIXSIZE], int *weights, int *scoreArr, int start, int end);