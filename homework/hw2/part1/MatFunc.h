/* File: MatFunc.h
 * Author: Joseph Moore
 * Purpose: more functions for matrices
 */

#ifndef MATFUNC_H
#define MATFUNC_H

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <stdbool.h>

#include "ColMajorMat.h"

/* We will only use a pointer to the struct, so make that the typedef */
typedef struct _p_ColMajorMat *ColMajorMat;

double MatTrace(ColMajorMat);
double TwoNorm(ColMajorMat);
double TwoNormCol(ColMajorMat,double);
void PrintMat(ColMajorMat);
void PrintMatSciNot(ColMajorMat);
void GEPP(ColMajorMat,ColMajorMat,bool*);
void RowSwap(ColMajorMat,int,int);
void BackSub(ColMajorMat,ColMajorMat);
ColMajorMat MatMatMult(ColMajorMat,ColMajorMat);
void MatMatSub(ColMajorMat,ColMajorMat);

#endif
