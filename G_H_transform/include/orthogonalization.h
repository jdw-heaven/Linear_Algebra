#ifndef _ORTHOGONALIZATION_H_	
#define _ORTHOGONALIZATION_H_

#include <stdlib.h>
#include <complex.h>
#include <stdio.h>
#include <math.h>
#include "tensor_cal.h"

void m_orth(double complex *T, int n){
    // orth. the next n-1 column
    double complex *a;
    double complex ipro1;
    double complex mod;
    a = (double complex *)malloc(n*sizeof(double complex));
    for(int i=1; i<n; i++){
        for(int j = 0; j < n; j++){
            a[j] = 0+0*I;
        }
        for(int j = 0; j < i; j++){
            ipro1 = 0+0*I;
            for(int k = 0; k < n; k++){
                ipro1 += conj(T[k*n+j])*T[k*n+i];
            }
            for(int k = 0; k < n; k++){
                a[k] += ipro1*T[k*n+j];
            }
        }
        for(int j = 0; j < n; j++){
            T[j*n+i] = T[j*n+i] - a[j];
        }
        mod = 0+0*I;
        for(int j = 0; j < n; j++){
            mod += conj(T[j*n+i])*T[j*n+i];
        }
        mod = csqrt(mod);
        printf("mod %3.4lf+i%3.4lf\n", creal(mod), cimag(mod));
        for(int j = 0; j < n; j++){
            T[j*n+i] = T[j*n+i]/mod;
        }
    }

    free(a);
}

#endif