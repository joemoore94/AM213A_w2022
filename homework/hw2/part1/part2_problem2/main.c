/* File: main.c
 * Author: Joseh Moore
 * Purpose: testing out numerical limits
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

#include "../ColMajorMat.h"
#include "../MatFunc.h"

int main(int argc, char const *argv[]) {
    ColMajorMat A = malloc(sizeof(*A));
    A->m = 2; A->n = 2;
    A->data = calloc(A->m*A->n,sizeof(*A->data));
    A->data[MatIdx(A,0,0)] = 1;
    A->data[MatIdx(A,0,1)] = 1;
    A->data[MatIdx(A,1,0)] = 1e-17;
    A->data[MatIdx(A,1,1)] = 1;

    ColMajorMat b = malloc(sizeof(*b));
    b->m = 2; b->n = 1;
    b->data = calloc(b->m*b->n,sizeof(*b->data));
    b->data[MatIdx(b,0,0)] = 2;
    b->data[MatIdx(b,1,0)] = 1;

    PrintMatSciNot(A);

    // A->data[MatIdx(A,1,0)] *= 1e18;
    // A->data[MatIdx(A,1,1)] *= 1e18;
    //
    // b->data[MatIdx(b,1,0)] *= 1e18;

    PrintMatSciNot(A);

    bool singular;
    GEPP(A, b, &singular);
    // BackSub(A, b);

    PrintMatSciNot(A);
    PrintMatSciNot(b);

    return 0;
}
