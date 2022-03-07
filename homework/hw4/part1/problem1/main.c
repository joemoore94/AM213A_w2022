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


    ColMajorMat A = CreateMatrix("../Amat.dat");

    PrintMat(A);

	Hessenberg(A);
	Hessenberg(A);
    Hessenberg(A);
    Hessenberg(A);

    PrintMat(A);

    return 0;
}
