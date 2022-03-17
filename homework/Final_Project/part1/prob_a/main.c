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

extern void dgesvd_(char*,char*,int*,int*,double*,int*,double*,double*,int*,double*,int*,double*,int*,int*);

int main(int argc, char const *argv[]) {

    ColMajorMat A = CreateMatrixImage("../dog_bw_data.dat", 3355, 5295);
    // ColMajorMat A = CreateMatrixImage("../dog_bw_data.dat_811x1280", 3355, 5295);

    /* assign input values */
    char jobu='A', jobvt='A';
    int n = A->n, m = A->m;

    /* Space to store the singular values */
    double Svals[n];

    /* create U and VT matrices */
    int ldu = m, ldvt = n;
    double *u = malloc(m*m*sizeof(*u));
    double *vt = malloc(n*n*sizeof(*vt));

    int lwork = -1;
    double *work = malloc(sizeof(*work));
    int info = 0;

    /* calculate work space */
    dgesvd_(&jobu,&jobvt,&m,&n,A->data,&n,Svals,u,&ldu,vt,&ldvt,work,&lwork,&info);
    lwork = (int)work[0];
    free(work);
    work = malloc(10*lwork*sizeof(*work));

    /* Call dgesvd */
    dgesvd_(&jobu,&jobvt,&m,&n,A->data,&m,Svals,u,&ldu,vt,&ldvt,work,&lwork,&info);

    // for (int i = 0; i < 10; i++) {
    //     printf("k = %d: %f\n", i+1, Svals[i]);
    // }
    // printf("k = %d: %f\n", 20, Svals[19]);
    // printf("k = %d: %f\n", 39, Svals[39]);
    // printf("k = %d: %f\n", 79, Svals[79]);
    // printf("k = %d: %f\n", 159, Svals[159]);
    // printf("k = %d: %f\n", 319, Svals[319]);
    // printf("k = %d: %f\n", 639, Svals[639]);
    // printf("k = %d: %f\n", 1279, Svals[1279]);
    // printf("k = %d: %f\n", 2559, Svals[2559]);
    // printf("k = %d: %f\n", 3354, Svals[3354]);

    /* trun u into a matrix*/
    ColMajorMat U = malloc(sizeof(*U));
    U->m = m; U->n = m;
    U->data = calloc(U->m*U->n,sizeof(*U->data));
    U->data = u;

    /* turn vt into a matrix*/
    ColMajorMat VT = malloc(sizeof(*VT));
    VT->m = n; VT->n = n;
    VT->data = calloc(VT->m*VT->n,sizeof(*VT->data));
    VT->data = vt;

    ColMajorMat S = malloc(sizeof(*S));
    S->m = m; S->n = n;
    S->data = calloc(S->m*S->n,sizeof(*S->data));

    int k = 20;
    // int k = 3355;
    for (int col = 0; col < n; col++) {
        for (int row = 0; row < m; row++) {
            if (col == row && col < k) S->data[MatIdx(S,row,col)] = Svals[col];
            else S->data[MatIdx(S,row,col)] = 0;
        }
    }

    ColMajorMat SVT;
    ColMajorMat USVT;
    char filename[100];

    for (int i = 0; i < 8; i++) {
        for (int col = 0; col < n; col++) {
            for (int row = 0; row < m; row++) {
                if (col == row && col < k) S->data[MatIdx(S,row,col)] = Svals[col];
                else S->data[MatIdx(S,row,col)] = 0;
            }
        }

        SVT = MatMatMult(S, VT);
        USVT = MatMatMult(U, SVT);

        sprintf(filename, "Image_%d.dat", k);

        SaveMatrix(USVT, filename);

        k *= 2;
    }

    /* calculating error*/
    // ColMajorMat AS = CreateMatrixImage("Image_2560.dat", 3355, 5295);
    // MatMatSub(A,AS);
    // double sum = 0;
    // for (int i = 0; i < A->m; i++) {
    //     for (int j = 0; j < A->n; j++) {
    //         sum += pow(A->data[MatIdx(A,i,j)], 2);
    //     }
    // }
    // sum = sqrt(sum) / (A->m*A->n);
    // printf("%f\n", sum);

    return 0;
}
