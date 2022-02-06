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
    ColMajorMat A = CreateMatrix("../Amat.dat");

    PrintMat(A, A->m, A->n);

    double trace = 0;
    MatTrace(A, A->m, &trace);
    printf("\ntrace(A) = %f\n\n", trace);

    for (int i = 0; i < A->n; i++) {
        double norm = 0;
        TwoNormCol(A, A->m, i, &norm);
        printf("||A%d|| = %f\n", i, norm);
    }

    return 0;
}
