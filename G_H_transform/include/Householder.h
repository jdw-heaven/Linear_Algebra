#ifndef _HOUSEHOLDER_H_
#define _HOUSEHOLDER_H_

#include <stdlib.h>
#include <complex.h>
#include <stdio.h>
#include <math.h>
#include "tensor_cal.h"

double complex * m_Householder(double complex *xi, double complex *eta, int D){
    double complex *U;
    U = (double complex *)malloc(D*D*sizeof(double complex));
    for(int i = 0; i < D*D; i++){
        U[i] = 0+0*I;
    }

    double complex z, mod;
    double complex *omega;
    omega = (double complex *)malloc(D*sizeof(double complex));
    z = m_ipro(xi, eta, D);
    z = z/cabs(z);
    for(int i = 0; i < D; i++){
        omega[i] = z*xi[i]-eta[i];
    }
    mod = csqrt(m_ipro(omega, omega, D));
    for(int i = 0; i < D; i++){
        omega[i] = omega[i]/mod;
    }
    for(int i = 0; i < D; i++){
        for(int j = 0; j < D; j++){
            U[i*D+j] = z*((i==j)-2*omega[i]*conj(omega[j]));
        }
    }

    //check
    printf("U(Householder) is :\n");
    m_cprint(U, D, D);
    printf("isUnitary? Udag*U is :\n");
    m_isUnitary(U, D);
    printf("Uxi is :\n");
    m_cprint(m_mul(U, xi, D, D, 1), D, 1);
    printf("well, a great success. its equ. to eta.\n\n");

    free(omega);
    free(U);
    return U;
}









#endif