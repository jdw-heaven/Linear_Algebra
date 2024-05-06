#ifndef _MY_SVD_H
#define _MY_SVD_H

#include <stdio.h>  // 根据需要包含其他的头文件
#include <complex.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>

//print tensor
void m_cprint(double complex *A, int row, int column){
    for(int i = 0; i < row; i++){
        for(int j = 0; j < column; j++){
            printf("%3.4lf+i%3.4lf  ", creal(A[i*column+j]), cimag(A[i*column+j]));
        }printf("\n");
    }printf("\n");
}
void m_print(double *A, int row, int column){
    for(int i = 0; i < row; i++){
        for(int j = 0; j < column; j++){
            printf("%3.4lf  ", A[i*column+j]);
        }printf("\n");
    }printf("\n");
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
// expand A in n*n complex space with r orthogonal normalisation volumn to n volumn 
void m_expand(int D, int r, double complex *A){
    for(int i=r+1; i<D; i++){
        A[i*D+i] = 1;
    }
    //orthognal
    double complex *med;
    med = (double complex*)malloc(D*sizeof(double complex));
    double complex up, down;
    for(int i=r+1; i<D; i++){
        for(int j=0; j<D; j++){
            med[j] = 0+0*I;
        }
        for(int j=0; j<i; j++){
            up = 0+0*I; down = 0+0*I;
            for(int k=0; k<D; k++){
                up += conj(A[k*D+j])*A[k*D+i];
                down += conj(A[k*D+j])*A[k*D+j];
            }
            for(int k=0; k<D; k++){
                med[k] += up*A[k*D+j]/down;
            }
        }
        for(int j=0; j<D; j++){
            A[j*D+i] = A[j*D+i]-med[j];
        }
    }
    //normalize
    double mod;
    for(int i=r+1; i<D; i++){
        mod = 0;
        for(int j=0; j<D; j++){
            mod += conj(A[j*D+i])*A[j*D+i];
        }
        mod = sqrt(mod);
        for(int j=0; j<D; j++){
            A[j*D+i] = A[j*D+i]/mod;
        }
    }
    free(med);
}

// singular value decomposition, A = U*S*Vdag
void m_SVD(int ROW, int COL, double complex *A, double complex *U, double complex *Vdag, double complex *S){
    // calculate AdagA
    for(int i=0; i<COL; i++){
        for(int j=0; j<ROW; j++){
            for(int k=0; k<COL; k++){
                Vdag[i*COL+k] += conj(A[j*COL+i])*A[j*COL+k];
            }
        }
    }
    double *E;
    E = (double*)malloc(COL*sizeof(double));
    //m_cprint(Vdag,COL,COL);

    int info = LAPACKE_zheev(LAPACK_COL_MAJOR, 'V', 'U', COL, Vdag, COL, E);
    if (info != 0) {
        fprintf(stderr, "LAPACKE_zheev returned info=%d\n", info);
        free(E);
        printf("error!");
    }
    m_cprint(Vdag,COL,COL);
    printf("The singular value's square(the eigenvalues of AdagA):\n");
    m_print(E,1,COL);
    // the number of singular value
    int r = 0;
    for(int i=0; i<COL; i++){
        if(E[COL-1-i] < 1e-10){
            r = i-1;
            break;
        }
        S[i*COL+i] = csqrt(E[COL-1-i]);
        r = i;
    }
    printf("%d\n",r);
    printf("Matrix Sigma is:\n");
    m_cprint(S,ROW,COL);
    //get Vdag(exchange Vdag's row, as the result is range as row) 
    double complex med;
    for(int i=0; i<(int)(COL/2); i++){
        for(int j=0; j<COL; j++){
            med = Vdag[i*COL+j];
            Vdag[i*COL+j] = Vdag[(COL-1-i)*COL+j];
            Vdag[(COL-1-i)*COL+j] = med;
        }
    }
    // complex conjugation
    /*
    for(int i=0; i<COL; i++){
        for(int j=0; j<COL; j++){
            Vdag[i*COL+j] = conj(Vdag[i*COL+j]);
        }
    }*/ //why not need?
    printf("Matrix Vdag is:\n");
    m_cprint(Vdag,COL,COL);

    //get U
    for(int i=0; i<ROW; i++){
        for(int j=0; j<r+1; j++){
            for(int k=0; k<COL; k++){
                U[i*ROW+j] += A[i*COL+k]*conj(Vdag[j*COL+k])/S[j*COL+j];
            }
        }
    }
    m_expand(ROW,r,U);
    printf("Matrix U is: \n");
    m_cprint(U,ROW,ROW);
    free(E);
}


#endif 