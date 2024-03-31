#include <stdlib.h>
#include <complex.h>
#include <stdio.h>

void cprint(double complex *Z){
    int n = 3;
    for(int i = 0; i < n; i++){
        printf("%lf + i%lf\n", creal(Z[i]), cimag(Z[i]));
    }
}
void m_orth(double complex *T, int n){
    double length = 0;
    for(int i = 0; i < n; i++){

    }
}
int main(void)
{
    double complex z1 = 1+4*I;
    double complex z2 = 2+1*I;
    double complex z3;
    z3 = z1/2;
    double complex *Z;
    Z = (double complex *)malloc(3*sizeof(double complex));
    Z[0] = z1;Z[1] = z2;Z[2] = z3;
    cprint(Z);
    printf("%lf + i%lf\n", creal(z3), cimag(z3));
    return 0;
}