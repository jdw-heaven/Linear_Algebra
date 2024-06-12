#ifndef _HEADER_FILE_H_
#define _HEADER_FILE_H_

//the fundamental header file need
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <complex.h>
#include <lapacke.h>
#include "/home/heaven/Desktop/doc/Linear_Algebra/m_header_file/random_num/mt19937ar-master/mt19937ar.c"

const int L = 12;   //单周期内的粒子数
const double error = 1e-10; //数值零

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

//inner product :c = <a, b>
double complex m_ipro(double complex *a, double complex *b, int n){
    double complex c;
    c = 0 + 0*I;

    for(int i = 0; i < n; i++){
        c += conj(a[i])*b[i];
    }

    return c;
}

// 将 a 的第 b 位取反 ，最低位编号为 0
int flapBit(int a, int b) { return a ^ (1 << b); }
// 获取 a 的第 b 位，最低位编号为 0
int getBit(int a, int b) { return (a >> b) & 1; }
// H作用在f上
void m_Hf(double complex *v, double complex *f, double Delta){
    // v一定要是空向量
    int num = 0;
    for(int j=0; j<(1<<L); j++){
        for(int i=0; i<L; i++){
            if(getBit(j,i)==getBit(j,(i+1)%L)){
                v[j] += Delta*f[j]/4;
            }else{
                v[j] -= Delta*f[j]/4;
                num = flapBit(flapBit(j,(i+1)%L),i);
                v[num] += f[j]/2;
            }    
        }
    }
}


//lanczos iteration
double m_lanczos(double Delta){

    double* alpha;
    alpha = (double*)calloc((1<<L),sizeof(double));
    double* beta;
    beta = (double*)calloc((1<<L),sizeof(double));

    double complex* f;
    f = (double complex*)calloc((1<<L),sizeof(double complex));
    
    double complex* v;
    v = (double complex*)calloc((1<<L),sizeof(complex double));
    double complex* fp;
    fp = (double complex*)calloc((1<<L),sizeof(double complex));
    double E = 0;
    double Ep = 1;
    
    for(int i=0; i< 1<<L; i++){
        f[i] = genrand_real1()+genrand_real1()*I;
    }
    double norm = sqrt(creal(m_ipro(f,f,1<<L)));
    //printf("%lf\n", norm);
    for(int i=0; i< 1<<L; i++){
        f[i] = f[i]/norm;
    }

    int n = 0;
    do{
        n += 1;
        m_Hf(v,f,Delta);
        for(int i=0; i< 1<<L; i++){
            v[i] = v[i] - beta[n-1]*fp[i];
        }//更好的数值稳定性
        alpha[n] = creal(m_ipro(f,v,1<<L));
        //printf("%lf\n", alpha[n]);
        for(int i=0; i< 1<<L; i++){
            v[i] = v[i] - alpha[n]*f[i];
        }
        beta[n] = sqrt(creal(m_ipro(v,v,1<<L)));
        //printf("%lf\n", beta[n]);
        for(int i=0; i< 1<<L; i++){
            fp[i] = f[i];
            f[i] = v[i]/beta[n];
            v[i] = 0+0*I;
        }
        if(n%17==0){
            Ep = E;
            //求解本征值
            double* A;
            A = (double*)calloc(n*n,sizeof(double));
            for(int i=0; i<n; i++){
                A[i*n+i] = alpha[i+1];
            }
            for(int i=0; i<n-1; i++){
                A[i*n+i+1] = beta[i+1];
                A[(i+1)*n+i] = beta[i+1];
            }
            double *a;
            a = (double *)calloc(n,sizeof(double));
            LAPACKE_dsyev(LAPACK_COL_MAJOR,'N','U',n,A,n,a);
            E = a[0];
            free(a);
            free(A);
        }
    }while( abs(E-Ep) > error );




    free(f);
    free(fp);
    free(v);
    free(alpha);
    free(beta);

    return E;
}


#endif