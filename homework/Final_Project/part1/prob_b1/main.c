/* File: main.c
 * Author: Joseph Moore
 * Purpose: SVD
 */

#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <math.h>

#include "../ColMajorMat.h"
#include "../MatFunc.h"

int main(int argc, char const *argv[]) {

    ColMajorMat A = CreateMatrix("../Amat.dat");
    ColMajorMat b = CreateMatrix("../Bmat.dat");

    int choice;
    printf("Which algorithm would you like to use?\n(1) Gauss-Jacobi\n(2) Gauss-Seidel\n");
    printf("Enter choice: ");
    scanf("%d", &choice);

    // int number;
    // printf("Enter integer value of diagonal elements: ");
    // scanf("%d", &number);
    //
    // for (int i = 0; i < A->m; i++) {
    //     A->data[MatIdx(A,i,i)] = number;
    // }

    for (int i = 0; i < A->m; i++) {
        A->data[MatIdx(A,i,i)] = i+1;
    }

    ColMajorMat x;

    switch(choice) {
        case 1:
            x = GaussJacobi(A, b, 1e-5);
            printf("\nx\n");
            PrintMat(x);
            break;

        case 2:
            x = GaussSeidel(A, b, 1e-5);
            printf("\nx\n");
            PrintMat(x);
            break;

        default:
            printf("Invalid choice\n" );
    }

    return 0;
}
