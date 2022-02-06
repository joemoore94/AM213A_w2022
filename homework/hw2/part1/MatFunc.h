/* File: MatFunc.h
 * Author: Joseph Moore
 * Purpose: more functions for matrices
 */

#ifndef MATFUNC_H
#define MATFUNC_H

/* We will only use a pointer to the struct, so make that the typedef */
typedef struct _p_ColMajorMat *ColMajorMat;

void MatTrace(ColMajorMat,double,double*);
void TwoNorm(ColMajorMat,double,double*);
void PrintMat(ColMajorMat,int,int);

#endif
