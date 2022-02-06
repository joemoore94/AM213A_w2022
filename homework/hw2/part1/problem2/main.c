/* File: main.c
 * Author: Joseh Moore
 * Purpose: test matrix uploader
 */

#include <stdlib.h>
#include <stdio.h>

#include "../ColMajorMat.h"
#include "../MatFunc.h"

int main(int argc, char const *argv[]) {
    ColMajorMat A = CreateMatrix("../Amat.dat");
    ColMajorMat C = CreateMatrix("../Cmat.dat");

    double trace = 0;
    MatTrace(A, A->m, &trace);
    printf("%f\n", trace);

    double norm = 0;
    TwoNorm(C, C->m, &norm);
    printf("%f\n", norm);

    PrintMat(A, A->m, A->n);

    return 0;
}
