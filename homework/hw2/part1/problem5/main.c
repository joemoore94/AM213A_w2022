/* File: main.c
 * Author: Joseh Moore
 * Purpose:
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

#include "../ColMajorMat.h"
#include "../MatFunc.h"

int main(int argc, char const *argv[]) {

    // ColMajorMat A = CreateMatrix("../Amat.dat");
    ColMajorMat C = CreateMatrix("../Cmat.dat");
    ColMajorMat D = CreateMatrix("../Dmat.dat");
    C->data[MatIdx(C,2,0)] = M_PI;
    C->data[MatIdx(C,2,1)] = M_E;
    C->data[MatIdx(C,2,2)] = -M_SQRT2;

    bool singular;
    GEPP(C, D, &singular);
    BackSub(C, D);

    PrintMat(D);

    system("python planeplot.py");

    return 0;
}
