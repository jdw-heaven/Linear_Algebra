#include "include/orthogonalization.h"
#include "include/Givens.h"
#include "include/Householder.h"
#include "include/tensor_cal.h"

#include "../m_header_file/random_num/mt19937ar-master/mt19937ar.c"

#define D 8

int main(void)
{
    init_genrand((unsigned long)time(NULL));
    double mod_define = 10*genrand_real3();
    printf("mod_define is %3.4lf\n", mod_define);
    double complex *xi;
    double complex *eta;
    xi = (double complex *)malloc(D*sizeof(double complex));
    eta = (double complex *)malloc(D*sizeof(double complex));
    // init 
    for(int i = 0; i < D; i++){
        eta[i] = 10*genrand_real3() + 10*genrand_real3()*I;
        xi[i] = 10*genrand_real3() + 10*genrand_real3()*I;
    }
    // equ. mod
    double complex mod_xi = csqrt(m_ipro(xi, xi, D));
    double complex mod_eta = csqrt(m_ipro(eta, eta, D));
    //check
    printf("initial mod_xi is %3.4lf+i%3.4lf, mod_eta is %3.4lf+i%3.4lf\n", creal(mod_xi), cimag(mod_xi), creal(mod_eta), cimag(mod_eta));

    for(int i = 0; i < D; i++){
        eta[i] = eta[i]*mod_define/mod_eta;
        xi[i] = xi[i]*mod_define/mod_xi;
    }
    mod_xi = csqrt(m_ipro(xi, xi, D));
    mod_eta = csqrt(m_ipro(eta, eta, D));
    printf("eta:\n");
    m_cprint(eta, D, 1);
    printf("final mod_xi is %3.4lf+i%3.4lf, mod_eta is %3.4lf+i%3.4lf\n", creal(mod_xi), cimag(mod_xi), creal(mod_eta), cimag(mod_eta));

    //Givens 变换
    printf("Gievens trans.\n");
    m_Givens(xi, eta, D);

    //Householder 变换
    printf("Householder trans.\n");
    m_Householder(xi, eta, D);


    free(xi);
    free(eta);
    return 0;
}