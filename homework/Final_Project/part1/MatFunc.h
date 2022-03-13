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
double TwoNormCol(ColMajorMat,int);
void PrintMat(ColMajorMat);
void PrintMatSciNot(ColMajorMat);
void GEPP(ColMajorMat,ColMajorMat,bool*);
void RowSwap(ColMajorMat,int,int);
void BackSub(ColMajorMat,ColMajorMat);
void BackSubLU(ColMajorMat,ColMajorMat,ColMajorMat);
ColMajorMat MatMatMult(ColMajorMat,ColMajorMat);
void MatMatSub(ColMajorMat,ColMajorMat);
ColMajorMat LUDecomp(ColMajorMat,bool*);

void CF(ColMajorMat,bool*,bool*);
void CFBacksub(ColMajorMat, ColMajorMat);
void Householder(ColMajorMat, ColMajorMat);
void HHBacksub(ColMajorMat, ColMajorMat, ColMajorMat);
void Transpose(ColMajorMat);
ColMajorMat T(ColMajorMat);

void Hessenberg(ColMajorMat);
ColMajorMat I(double, double);
ColMajorMat InverseIt(ColMajorMat, double);
void CopyMat(ColMajorMat, ColMajorMat);

#endif
