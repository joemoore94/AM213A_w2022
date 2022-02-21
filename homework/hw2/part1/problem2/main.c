/* File: main.c
 * Author: Joseh Moore
 * Purpose: test some basic functions
 */

#include <stdlib.h>
#include <stdio.h>

#include "../ColMajorMat.h"
#include "../MatFunc.h"

/*
Pretty straight forward space separated answers output to the terminal
*/

int main(int argc, char const *argv[]) {


    // You can swap the comment on the two rows below and exicute: ./run ../Amat.dat
    // ColMajorMat A = CreateMatrix(argv[1]);
    ColMajorMat A = CreateMatrix("../Amat.dat");

    PrintMat(A);

    double trace = MatTrace(A);
    printf("\ntrace(A) = %f\n\n", trace);

    for (int i = 0; i < A->n; i++) {
        printf("||A%d|| = %f\n", i, TwoNormCol(A, i));
    }

    return 0;
}
