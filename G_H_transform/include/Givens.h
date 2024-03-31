#ifndef _GIVENS_H_
#define _GIVENS_H_


#include <stdlib.h>
#include <complex.h>
#include <stdio.h>
#include <math.h>
#include "tensor_cal.h"


// xi to eta
double complex * m_Givens(double complex *xi, double complex *eta, int D){
    double complex *U;
    double complex *Up;
    U = (double complex *)malloc(D*D*sizeof(double complex));
    Up = (double complex *)malloc(D*D*sizeof(double complex));
    for(int i = 0; i < D*D; i++){
        U[i] = 0+0*I;
        Up[i] = 0+0*I;
    }

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


    //Givens 旋转
    double complex c, s, z1, z2;
    double theta1, theta2;
    double complex *Q;
    Q = (double complex *)malloc(D*D*sizeof(double complex));
    for(int i = 0; i < D; i++){
        Up[i*D+i] = 1+0*I;
    }
    for(int i = D-1; i > 0; i--){
        c = cabs(xip[i-1])/sqrt(cabs(xip[i-1])*cabs(xip[i-1])+cabs(xip[i])*cabs(xip[i]));
        s = cabs(xip[i])/sqrt(cabs(xip[i-1])*cabs(xip[i-1])+cabs(xip[i])*cabs(xip[i]));
        theta1 = -carg(xip[i-1]);
        theta2 = -carg(xip[i]);
        z1 = cos(theta1)+sin(theta1)*I;
        z2 = cos(theta2)+sin(theta2)*I;
        for(int i = 0; i < D*D; i++){
            Q[i] = 0;
        }
        for(int i = 0; i < D; i++){
            Q[i*D+i] = 1+0*I;
        }
        Q[(i-1)*D+i-1] = c*z1;
        Q[(i-1)*D+i] = s*z2;
        Q[(i)*D+i-1] = -s*conj(z2);
        Q[(i)*D+i] = c*conj(z1);
        Up = m_mul(Q, Up, D, D, D);

        xip[i-1] = sqrt(cabs(xip[i-1])*cabs(xip[i-1])+cabs(xip[i])*cabs(xip[i]));
        xip[i] = 0;
    }
    U = m_mul(T, Up, D, D, D);
    U = m_mul(U, Tdag, D, D, D);

    //check
    printf("U(Givens) is :\n");
    m_cprint(U, D, D);

    printf("Uxi is :\n");
    m_cprint(m_mul(U, xi, D, D, 1), D, 1);
    printf("well, a great success. its equ. to eta.\n\n");

    free(Q);
    free(U);
    free(Up);
    free(xip);
    free(etap);
    free(T);
    free(Tdag);
    return U;
}

#endif