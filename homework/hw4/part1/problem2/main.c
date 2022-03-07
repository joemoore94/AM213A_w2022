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


    ColMajorMat A = CreateMatrix("../Amat1.dat");

    printf("original matrix\n");
    PrintMat(A);

    ColMajorMat Q = malloc(sizeof(*Q));
    Q->m = A->m; Q->n = A->m;
    Q->data = calloc(Q->m*Q->n,sizeof(*Q->data));


    // QR without shift
    double error = 1;
    int i = 0;
    while (error > 1e-8) {
        Householder(A, Q);
        A = MatMatMult(A, Q);

        error = fabs(A->data[MatIdx(A,0,0)] - TwoNormCol(A, 0));
        i++;
    }

    printf("\n%d", i);

    printf("\nAfter QR\n");
    PrintMat(A);


    // QR with shift
    A = CreateMatrix("../Amat1.dat");

    double mu = A->data[MatIdx(A,A->m-1,A->n-1)];
    error = 1;
    i = 0;
    while (error > 1e-8) {
        MatMatSub(A, I(A->m, mu));
        Householder(A, Q);
        A = MatMatMult(A, Q);
        MatMatSub(A, I(A->m, -mu));

        error = fabs(A->data[MatIdx(A,0,0)] - TwoNormCol(A, 0));
        mu = A->data[MatIdx(A,A->m-1,A->n-1)];
        i++;
    }

    printf("\n%d", i);

    printf("\nAfter QR with shift\n");
    PrintMat(A);

    return 0;
}
