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

    system("python least_squares_prep.py");

    ColMajorMat A = CreateMatrix("../Amat.dat");
    ColMajorMat B = CreateMatrix("../Bmat.dat");

    printf("Original matrices\n");
    printf("A: ");
    PrintMat(A);
    printf("B: ");
    PrintMat(B);

    bool singular;
    bool SPD;
    CF(A, &singular, &SPD);
    printf("\nAfter Cholesky factorization\n");
    printf("A: ");
    PrintMat(A);
    printf("B: ");
    PrintMat(B);

    CFBacksub(A, B);
    printf("\nAfter backsubstitution\n");
    printf("A: ");
    PrintMat(A);
    printf("B: ");
    PrintMat(B);

    FILE *fptr;
    fptr = fopen("../Cmat.dat","w");
    fprintf(fptr,"%d 1\n", B->m);
    for (int i = 0; i < B->m; i++) {
        fprintf(fptr,"%f\n",B->data[MatIdx(B,i,0)]);
    }
    fclose(fptr);

    B = CreateMatrix("../Bmat.dat");
    A = CreateMatrix("../Amat.dat");
    MatMatSub(B, MatMatMult(A, CreateMatrix("../Cmat.dat")));
    printf("\nError vector (b - Ax)\n");
    PrintMatSciNot(B);
    printf("\n2-norm error: %e\n", TwoNorm(B));

    system("python results_plot.py");

    return 0;
}
