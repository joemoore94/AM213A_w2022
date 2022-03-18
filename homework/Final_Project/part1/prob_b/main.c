/* File: main.c
 * Author: Joseph Moore
 * Purpose: iterative methods
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
    printf("(3) Conjugate Gradientn\n");
    printf("Enter choice: ");
    scanf("%d", &choice);

    int number;
    printf("Enter integer value of diagonal elements: ");
    scanf("%d", &number);

    for (int i = 0; i < A->m; i++) {
        A->data[MatIdx(A,i,i)] = number;
    }

    // // 10 by 10 matrix
    // for (int i = 0; i < A->m; i++) {
    //     A->data[MatIdx(A,i,i)] = i+1;
    // }

    // // 100 by 100 matrix
    // A->m = 100; A->n = 100;
    // free(A->data);
    // A->data = calloc(A->m*A->n,sizeof(*A->data));
    // for (int i = 0; i < A->m; i++) {
    //     for (int j = 0; j < A->m; j++) {
    //         A->data[MatIdx(A,i,j)] = 1;
    //     }
    // }
    // for (int i = 0; i < A->m; i++) {
    //     A->data[MatIdx(A,i,i)] = i+1;
    // }

    ColMajorMat x = malloc(sizeof(*x));

    switch(choice) {
        case 1:
            x = GaussJacobi(A, b, 1e-5);
            break;

        case 2:
            x = GaussSeidel(A, b, 1e-5);
            break;

        case 3:
            x = CG(A, b, 1e-5);
            break;

        default:
            printf("Invalid choice\n" );
    }

    printf("\nx\n");
    PrintMat(x);

    return 0;
}
