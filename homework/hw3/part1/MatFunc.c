/* File: MatFunc.h
 * Author: Joseph Moore
 * Purpose: more functions for matrices
 */

#include "MatFunc.h"

/*
I avoided matrix dimensions as inputs since they are included in the data structure of the matrix.
The dimensions are being passed, just through the data structure.
*/

/*calculates the trace of a vector*/
float MatTrace(ColMajorMat A) {
    float trace = 0;
    for(int i = 0; i < A->m; i++) {
        trace = trace + A->data[MatIdx(A,i,i)];
    }
    return trace;
}

/*calculates the 2-norm of a vector*/
float TwoNorm(ColMajorMat A) {
    float norm = 0;
    for (int i = 0; i < A->m; i++) {
        norm = norm + A->data[MatIdx(A,i,0)]*A->data[MatIdx(A,i,0)];
    }
    return sqrt(norm);
}

/*calculates the 2-norm of a specified column vector of a matrix*/
float TwoNormCol(ColMajorMat A, float col) {
    float norm = 0;
    for (int i = 0; i < A->m; i++) {
        norm = norm + A->data[MatIdx(A,i,col)]*A->data[MatIdx(A,i,col)];
    }
    return sqrt(norm);
}

/*prints all the elements of a matrix to the screen*/
void PrintMat(ColMajorMat A) {
    printf("%dx%d Matrix\n", A->m, A->n);
    for (int i = 0; i < A->m; i++) {
        for (int j = 0; j < A->n; j++) {
            printf("%f ", A->data[MatIdx(A,i,j)]);
        }
        printf("\n");
    }
}

/*prints all the elements of a matrix to the screen in scientific notation*/
void PrintMatSciNot(ColMajorMat A) {
    printf("%dx%d Matrix\n", A->m, A->n);
    for (int i = 0; i < A->m; i++) {
        for (int j = 0; j < A->n; j++) {
            printf("%e ", A->data[MatIdx(A,i,j)]);
        }
        printf("\n");
    }
}

/*Gaussian elimination with Partial Pivoting:*/
void GEPP(ColMajorMat A, ColMajorMat B, bool* singular) {
    *singular = false;
    float eps = 0.001;
    float scale = 0;
    for (int col = 0; col < A->m; col++) {
        int maxIdx = col;
        for (int row = col; row < A->m; row++) {
            if (fabs(A->data[MatIdx(A,row,col)]) > fabs(A->data[MatIdx(A,maxIdx,col)])) {
                maxIdx = row;
            }
        }
        if (A->data[MatIdx(A,maxIdx,col)] < eps && A->data[MatIdx(A,maxIdx,col)] > -eps) {
            *singular = true;
            break;
        }
        RowSwap(A, col, maxIdx);
        RowSwap(B, col, maxIdx);
        for (int row = col+1; row < A->m; row++) {
            scale = A->data[MatIdx(A,row,col)]/A->data[MatIdx(A,col,col)];
            for (int i = 0; i < A->n; i++) {
                A->data[MatIdx(A,row,i)] -= scale*A->data[MatIdx(A,col,i)];
            }
            for (int i = 0; i < B->n; i++) {
                B->data[MatIdx(B,row,i)] -= scale*B->data[MatIdx(B,col,i)];
            }
        }
    }
}

/*zero indexed row swapping*/
void RowSwap(ColMajorMat A, int row1, int row2) {
    float temp;
    for (int i = 0; i < A->n; i++) {
        temp = A->data[MatIdx(A,row1,i)];
        A->data[MatIdx(A,row1,i)] = A->data[MatIdx(A,row2,i)];
        A->data[MatIdx(A,row2,i)] = temp;
    }
}

/*gives the solution of Ux = y, where U is an upper triagular matrix*/
void BackSub(ColMajorMat A, ColMajorMat B) {
    float sum;
    for (int Bcol = 0; Bcol < B->n; Bcol++) {
        B->data[MatIdx(B,B->m-1,Bcol)] /= A->data[MatIdx(A,A->m-1,A->n-1)];

        for (int Arow = A->m-2; Arow > -1; Arow--) {
            sum = 0.0;
            for (int i = Arow+1; i < A->n; i++) {
                sum += A->data[MatIdx(A,Arow,i)]*B->data[MatIdx(B,i,Bcol)];
            }
            B->data[MatIdx(B,Arow,Bcol)] -= sum;
            B->data[MatIdx(B,Arow,Bcol)] /= A->data[MatIdx(A,Arow,Arow)];
        }
    }
}

/*gives the solution of Ux = y, where U is an upper triagular matrix*/
void BackSubLU(ColMajorMat A, ColMajorMat B, ColMajorMat s) {
    ColMajorMat mat = malloc(sizeof(*mat));
    mat->m = B->m; mat->n = B->n;
    mat->data = calloc(mat->m*mat->n,sizeof(*mat->data));
    for (int row = 0; row < B->m; row++) {
        for (int col = 0; col < B->n; col++) {
            mat->data[MatIdx(mat,row,col)] = B->data[MatIdx(B,s->data[MatIdx(s,row,0)],col)];
        }
    }
    B->data = mat->data;
    for (int Bcol = 0; Bcol < B->n; Bcol++) {
        for (int j = 0; j < B->m-1; j++) {
            for (int i = j+1; i < B->m; i++) {
                B->data[MatIdx(B,i,Bcol)] -= B->data[MatIdx(B,j,Bcol)]*A->data[MatIdx(A,i,j)];
            }
        }
    }
    float sum;
    for (int Bcol = 0; Bcol < B->n; Bcol++) {
        B->data[MatIdx(B,B->m-1,Bcol)] /= A->data[MatIdx(A,A->m-1,A->n-1)];

        for (int Arow = A->m-2; Arow > -1; Arow--) {
            sum = 0.0;
            for (int i = Arow+1; i < A->n; i++) {
                sum += A->data[MatIdx(A,Arow,i)]*B->data[MatIdx(B,i,Bcol)];
            }
            B->data[MatIdx(B,Arow,Bcol)] -= sum;
            B->data[MatIdx(B,Arow,Bcol)] /= A->data[MatIdx(A,Arow,Arow)];
        }
    }
}

/*multiplies two matrices and returns a new matrix*/
ColMajorMat MatMatMult(ColMajorMat A, ColMajorMat B) {
    ColMajorMat mat = malloc(sizeof(*mat));
    mat->m = A->m; mat->n = B->n;
    mat->data = calloc(mat->m*mat->n,sizeof(*mat->data));
    float sum;
    for (int row = 0; row < mat->m; row++) {
        for (int col = 0; col < mat->n; col++) {
            sum = 0.0;
            for (int i = 0; i < A->n; i++) {
                sum += A->data[MatIdx(A,row,i)]*B->data[MatIdx(A,i,col)];
            }
            mat->data[MatIdx(mat,row,col)] = sum;
        }
    }
    return mat;
}

/*performs a matrix subtraction and replaces the first matrix with the results*/
void MatMatSub(ColMajorMat A,ColMajorMat B) {
    for (int row = 0; row < A->m; row++) {
        for (int col = 0; col < A->n; col++) {
            A->data[MatIdx(A,row,col)] -= B->data[MatIdx(B,row,col)];
        }
    }
}

/*decomposes a matrix into a lower and upper triangular matrix*/
ColMajorMat LUDecomp(ColMajorMat A, bool* singular) {
    *singular = false;
    float eps = 0.001;
    float scale = 0;

    ColMajorMat mat = malloc(sizeof(*mat));
    mat->m = A->m; mat->n = 1;
    mat->data = calloc(mat->m*mat->n,sizeof(*mat->data));
    for (int i = 0; i < mat->m; i++) {mat->data[MatIdx(mat,i,0)] = i;}

    for (int col = 0; col < A->m; col++) {
        // testing if A is singular
        int maxIdx = col;
        for (int row = col; row < A->m; row++) {
            if (fabs(A->data[MatIdx(A,row,col)]) > fabs(A->data[MatIdx(A,maxIdx,col)])) {
                maxIdx = row;
            }
        }
        if (A->data[MatIdx(A,maxIdx,col)] < eps && A->data[MatIdx(A,maxIdx,col)] > -eps) {
            *singular = true;
            break;
        }

        RowSwap(A, col, maxIdx);
        RowSwap(mat, col, maxIdx);
        for (int row = col+1; row < A->m; row++) {
            scale = A->data[MatIdx(A,row,col)]/A->data[MatIdx(A,col,col)];
            for (int i = col; i < A->n; i++) {
                A->data[MatIdx(A,row,i)] -= scale*A->data[MatIdx(A,col,i)];
            }
            A->data[MatIdx(A,row,col)] = scale;
        }
    }
    return mat;
}

/*decomposes a matrix into a lower triangular matrix using Cholesky factorization*/
void CF(ColMajorMat A, bool *singular, bool *SPD) {
    float eps = 0.001;
    *singular = false;
    *SPD = false;

    for (int j = 0; j < A->m; j++) {
        // testing if A is singular
        int maxIdx = j;
        for (int row = j; row < A->m; row++) {
            if (fabs(A->data[MatIdx(A,row,j)]) > fabs(A->data[MatIdx(A,maxIdx,j)])) {
                maxIdx = row;
            }
        }
        if (A->data[MatIdx(A,maxIdx,j)] < eps && A->data[MatIdx(A,maxIdx,j)] > -eps) {
            *singular = true;
            break;
        }

        for (int k = 0; k < j; k++) {
            A->data[MatIdx(A,j,j)] -= A->data[MatIdx(A,j,k)]*A->data[MatIdx(A,j,k)];
        }

        // testing if A is SPD
        if (A->data[MatIdx(A,j,j)] < 0) {
            *SPD = true;
            break;
        }

        A->data[MatIdx(A,j,j)] = sqrt(A->data[MatIdx(A,j,j)]);
        for (int i = j+1; i < A->m; i++) {
            for (int k = 0; k < j; k++) {
                A->data[MatIdx(A,i,j)] -= A->data[MatIdx(A,i,k)]*A->data[MatIdx(A,j,k)];
            }
            A->data[MatIdx(A,i,j)] /= A->data[MatIdx(A,j,j)];
        }
    }
}

void CFBacksub(ColMajorMat A, ColMajorMat B) {
    float sum;
    float eps = 0.001;

    for (int Bcol = 0; Bcol < B->m; Bcol++) {

        // Forward substitution, solving Ly = b
        for (int i = 0; i < A->m; i++) {
            sum = B->data[MatIdx(B,i,Bcol)];
            for (int j = 0; j < i; j++) {
                sum -= B->data[MatIdx(B,j,Bcol)]*A->data[MatIdx(A,i,j)];
            }
            B->data[MatIdx(B,i,Bcol)] = sum/A->data[MatIdx(A,i,i)];
        }

        // Backward substitution, solving L⇤x = y
        for (int i = A->m-1; i > -1; i--) {
            if (A->data[MatIdx(A,i,i)] < eps && A->data[MatIdx(A,i,i)] > -eps) {
                break;
            }
            for (int k = i+1; k < A->m; k++) {
                B->data[MatIdx(B,i,Bcol)] -= A->data[MatIdx(A,k,i)]*B->data[MatIdx(B,k,Bcol)];
            }
            B->data[MatIdx(B,i,Bcol)] /= A->data[MatIdx(A,i,i)];
        }
    }
}
