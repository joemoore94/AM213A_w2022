/* File: main.c
 * Author: Joseh Moore
 * Purpose: ref some matrices
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

#include "../ColMajorMat.h"
#include "../MatFunc.h"

int main(int argc, char const *argv[]) {

    ColMajorMat A = CreateMatrix("../Amat.dat");
    ColMajorMat B = CreateMatrix("../Bmat.dat");

    printf("Original matrices\n");
    printf("A: ");
    PrintMat(A);
    printf("B: ");
    PrintMat(B);

    bool singular;
    GEPP(A, B, &singular);

    printf("\nMatrices after gaussian elimination\n");
    printf("A: ");
    PrintMat(A);
    printf("B: ");
    PrintMat(B);

    BackSub(A, B);

    printf("\nAfter back-substitution\n");
    printf("A: ");
    PrintMat(A);
    printf("B: ");
    PrintMat(B);

    printf("\nError matix\n");
    ColMajorMat E = MatMatMult(CreateMatrix("../Amat.dat"), B);
    MatMatSub(E, CreateMatrix("../Bmat.dat"));
    PrintMatSciNot(E);

    printf("\nError matix column norms\n");
    for (int i = 0; i < E->n; i++) {
        printf("||E%d|| = %e\n", i, TwoNormCol(E, i));
    }

    return 0;
}
