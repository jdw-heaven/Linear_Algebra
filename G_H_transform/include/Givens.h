#include <stdlib.h>
#include <complex.h>
#include <stdio.h>
#include <math.h>
#include "tensor_cal.h"


// xi to eta
double complex * m_Givens(double complex *xi, double complex *eta, int D){
    double complex *U;
    double complex *Up;


    //坐标变换矩阵T的生成
    double complex mod_eta = csqrt(m_ipro(eta, eta, D));
    double complex *T;
    T = (double complex *)malloc(D*D*sizeof(double complex));
    for(int i = 0; i < D; i++){
        T[i*D] = eta[i]/mod_eta;
    }
    for(int i = 1; i < D; i++){
        T[(i-1)*D+i] = 1;
    }
    printf("the trans. matrix(non-Unitary)\n");
    m_cprint(T, D, D);
    //坐标变换
    m_orth(T, D);
    printf("the trans. matrix(Unitary).\n");
    m_cprint(T, D, D);
    //check is Unitary
    printf("check the T(trans. matrix) is Unitary.\n");
    m_isUnitary(T, D);

    //trans. 
    double complex *Tdag;
    Tdag = (double complex *)malloc(D*D*sizeof(double complex));
    for(int i = 0; i < D; i++){
        for(int j = 0; j < D; j++){
            Tdag[i*D+j] = conj(T[j*D+i]);
        }
    }
    double complex *xip;
    double complex *etap;
    xip = (double complex *)malloc(D*sizeof(double complex));
    etap = (double complex *)malloc(D*sizeof(double complex));
    xip = m_mul(Tdag, xi, D, D, 1);
    etap = m_mul(Tdag, eta, D, D, 1);
    printf("xip and etap.\n");
    m_cprint(xip, D, 1);
    m_cprint(etap, D, 1);
    //xip mod
    printf("xip mod is %lf\n", csqrt(m_ipro(xip, xip, D)));



    free(U);
    free(Up);
    free(xip);
    free(etap);
    free(T);
    free(Tdag);
    return U;
}