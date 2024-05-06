#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <lapacke.h>
#include "include/my_SVD.h"

int main(void)
{
    double *E;
    E = (double*)malloc(2*sizeof(double));
    double complex *A;
    A = (double complex*)malloc(4*sizeof(double complex));
    A[0]=0;A[1]=1*I;A[2]=-1*I;A[3]=0;
    m_cprint(A,2,2);
    int info = LAPACKE_zheev(LAPACK_ROW_MAJOR, 'V', 'U', 2, A, 2, E);
    if (info != 0) {
        fprintf(stderr, "LAPACKE_zheev returned info=%d\n", info);
        free(E);
        printf("error!");
    }
    m_print(E,1,2);
    m_cprint(A,2,2);

    A[0]=0;A[1]=1*I;A[2]=-1*I;A[3]=0;
    m_cprint(A,2,2);
    info = LAPACKE_zheev(LAPACK_COL_MAJOR, 'V', 'U', 2, A, 2, E);
    if (info != 0) {
        fprintf(stderr, "LAPACKE_zheev returned info=%d\n", info);
        free(E);
        printf("error!");
    }
    m_print(E,1,2);
    m_cprint(A,2,2);
    free(E);
    free(A);
    return 0;
}


