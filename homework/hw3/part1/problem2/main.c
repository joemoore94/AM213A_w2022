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


    ColMajorMat Q = malloc(sizeof(*Q));
    Q->m = A->m; Q->n = A->m;
    Q->data = calloc(Q->m*Q->n,sizeof(*Q->data));
    Householder(A, Q);
    ColMajorMat D = CreateMatrix("../Amat.dat");
    MatMatSub(D, MatMatMult(Q, A));
    printf("A - QR\n");
    PrintMatSciNot(D);


    float sum = 0;
    for (int i = 0; i < D->n; i++) {
        sum += pow(TwoNormCol(D, i), 2);
    }
    sum = sqrt(sum);
    printf("\n||A - QR|| = %e\n", sum);


    ColMajorMat QT = malloc(sizeof(*QT));
    QT->m = A->m; QT->n = A->m;
    QT->data = calloc(QT->m*QT->n,sizeof(*QT->data));
    ColMajorMat A2 = CreateMatrix("../Amat.dat");
    Householder(A2, QT);
    ColMajorMat QTQ = MatMatMult(QT, Q);
    ColMajorMat I = malloc(sizeof(*I));
    I->m = A->m; I->n = A->m;
    I->data = calloc(I->m*I->n,sizeof(*I->data));
    for (int col = 0; col < I->m; col++) {
        for (int row = 0; row < I->m; row++) {
            if (col == row) I->data[MatIdx(I,row,col)] = 1;
            else I->data[MatIdx(I,row,col)] = 0;
        }
    }
    MatMatSub(QTQ, I);
    printf("\nQ^T*Q - I\n");
    printf("Uncomment line 61 to see this matrix (its big!)\n");
    // PrintMatSciNot(QTQ);


    sum = 0;
    for (int i = 0; i < QTQ->n; i++) {
        sum += pow(TwoNormCol(QTQ, i), 2);
    }
    sum = sqrt(sum);
    printf("\n||Q^T*Q - I|| = %e\n", sum);

    HHBacksub(A, Q, B);
    ColMajorMat C = CreateMatrix("../Cmat.dat");
    printf("\nx\n");
    PrintMat(C);

    ColMajorMat Ax = MatMatMult(CreateMatrix("../Amat.dat"), C);
    MatMatSub(Ax, CreateMatrix("../Bmat.dat"));
    printf("\n2-norm: %f\n", TwoNorm(Ax));

    system("python results_plot.py");

    return 0;
}
