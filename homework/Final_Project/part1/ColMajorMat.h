/* File: ColMajorMat.h
 * Author: Ian May
 * Purpose: Another iteration on the previous example, now modeling the OOC
 *          approach
 */

#ifndef COLMAJORMAT_H
#define COLMAJORMAT_H

/* We will only use a pointer to the struct, so make that the typedef */
typedef struct _p_ColMajorMat *ColMajorMat;

/* Define a structure that holds the matrix dimensions and a flattened array */
struct _p_ColMajorMat{
  int m, n;
  double *data;
};

/* Functions to create and destroy the matrix struct */
ColMajorMat CreateMatrix(const char*);
ColMajorMat CreateMatrixImage(const char*,int,int);
void SaveMatrix(ColMajorMat,const char*);
void DestroyMatrix(ColMajorMat*);

/* Functions to interact with this matrix */
void FillMatrix(ColMajorMat);
void MatVecMult(ColMajorMat,double[*],double[*]);

/* Inlined element access */
static inline int MatIdx(ColMajorMat A,int i,int j)
{
  return A->m*j + i;
}


#endif
