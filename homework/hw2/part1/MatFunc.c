/* File: MatFunc.h
 * Author: Joseph Moore
 * Purpose: more functions for matrices
 */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "MatFunc.h"
#include "ColMajorMat.h"

void MatTrace(ColMajorMat A, double size, double *trace) {
    *trace = 0;
    for(int i = 0; i < size; i++) {
        *trace = *trace + A->data[MatIdx(A,i,i)];
    }
}

void TwoNorm(ColMajorMat A, double size, double* norm) {
    *norm = 0;
    for (int i = 0; i < size; i++) {
        *norm = *norm + A->data[MatIdx(A,i,0)]*A->data[MatIdx(A,i,0)];
    }
    *norm = sqrt(*norm);
}

void PrintMat(ColMajorMat A, int m, int n) {
    printf("%dx%d Matrix\n", m, n);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            printf("%f ", A->data[MatIdx(A,j,i)]);
        }
        printf("\n");
    }
}
