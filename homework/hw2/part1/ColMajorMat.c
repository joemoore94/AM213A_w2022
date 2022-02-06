/* File: ColMajorMat.c
 * Author: Ian May
 * Purpose: Another iteration on the previous example, now modeling the OOC
 *          approach
 */

#include <stdlib.h>
#include <stdio.h>

#include "ColMajorMat.h"

/* Matrix creation handles allocation */
ColMajorMat CreateMatrix(const char *fname)
{
  /* Allocate struct */
  ColMajorMat mat = malloc(sizeof(*mat));
  /* Open file and get matrix size */
  FILE *fp = fopen(fname,"r");
  if (fp) {
    fscanf(fp,"%d %d\n",&mat->m,&mat->n);
    /* Allocate space for data */
    mat->data = calloc(mat->m*mat->n,sizeof(*mat->data));
    /* Read in all values */
    char line[256];
    for (int i=0; i<mat->m; i++) {
      if(fgets(line,sizeof line,fp)) {
        int off = 0;
        char *scan = &line[0];
        for (int j=0; j<mat->n; j++) {
          if(sscanf(scan,"%lf%n",&mat->data[MatIdx(mat,i,j)],&off) > 0) {
            scan += off;
          } else {
            printf("Warning: row %d had too few entries\n",i);
          }
        }
      } else {
        break;
      }
    }
    /* Close file */
    fclose(fp);
    /* Pass this struct back out */
    return mat;
  } else {
    printf("File %s does not exist. Returning NULL\n",fname);
    free(mat);
    return NULL;
  }
}

/* Destruction frees everything creation allocated */
void DestroyMatrix(ColMajorMat *A)
{
  if (*A) {
    /* Free the data array inside the struct */
    free((*A)->data);
    /* Free the struct */
    free(*A);
  }
  *A = NULL;
}

/* function MatVecMult
   inputs:
     A : ColMajorMat, possibly non-square
     x : Vector compatible with shape of A
   outputs:
     b : Vector matching rows of A, does not need to be filled
   notes:
     This accesses A in the slow order, but is easier to read
*/
void MatVecMult(ColMajorMat A,double x[A->m],double b[A->m])
{
  for(int i=0; i<A->m; i++) {
    /* Zero out b before accumulating */
    b[i] = 0;
    /* dot x with each row */
    for(int j=0; j<A->n; j++) {
      b[i] += A->data[MatIdx(A,i,j)]*x[j];
    }
  }
}
