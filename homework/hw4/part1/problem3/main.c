/* File: main.c
 * Author: Joseph Moore
 * Purpose: using Cholesky factorization to fit a curve to some data
 */

#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <math.h>

#include "../ColMajorMat.h"
#include "../MatFunc.h"

/*
Pretty straight forward space separated answers output to the terminal
*/

int main(int argc, char const *argv[]) {


    ColMajorMat A = CreateMatrix("../Amat2.dat");
    A->data[MatIdx(A,1,1)] = -3;
    A->data[MatIdx(A,2,3)] = -2;
    A->data[MatIdx(A,3,2)] = -2;
    A->data[MatIdx(A,3,3)] = -1;
    printf("Original matrix\n");
    PrintMat(A);

    double lambda1 = -8.0286;
    // double lambda2 = 7.9329;
    // double lambda3 = 5.6689; // this one must be stoped manually 
    // double lambda4 = -1.5732;

    ColMajorMat v = InverseIt(A, lambda1);
    printf("\nEigenvector for %f\n", lambda1);
    PrintMat(v);

    return 0;
}
