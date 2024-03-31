#include <stdlib.h>
#include <complex.h>
#include <stdio.h>
#include <math.h>
//print tensor
void m_cprint(double complex *A, int row, int column){
    for(int i = 0; i < row; i++){
        for(int j = 0; j < column; j++){
            printf("%3.4lf+i%3.4lf  ", creal(A[i*column+j]), cimag(A[i*column+j]));
        }printf("\n");
    }printf("\n");
}

//inner product :c = <a, b>
double complex m_ipro(double complex *a, double complex *b, int n){
    double complex c;
    c = 0 + 0*I;

    for(int i = 0; i < n; i++){
        c += conj(a[i])*b[i];
    }

    return c;
}

//is unitary?
void m_isUnitary(double complex *A, int n){
    double complex *U;
    U = (double complex *)malloc(n*n*sizeof(double complex));
    for(int i = 0; i < n*n; i++){
        U[i] = 0+0*I;
    }
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            for(int k = 0; k < n; k++){
                U[i*n+j] += conj(A[k*n+i])*A[k*n+j];
            }
        }
    }
    m_cprint(U, n, n);
    free(U);
}

//multiply
double complex * m_mul(double complex *A, double complex *B, int Arow, int Acolumn, int Bcolumn){
    double complex *C;
    C = (double complex *)malloc(Arow*Bcolumn*sizeof(double complex));
    for(int i = 0; i < Arow*Bcolumn; i++){
        C[i] = 0+0*I;
    }
    for(int i = 0; i < Arow; i++){
        for(int j = 0; j < Bcolumn; j++){
            for(int k = 0; k < Acolumn; k++){
                C[i*Bcolumn+j] += A[i*Acolumn+k]*B[k*Bcolumn+j];
            }
        }
    }
    return C;
    free(C);
}