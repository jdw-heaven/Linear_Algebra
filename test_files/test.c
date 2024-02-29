#include"../m_header_file/m_complex.h"

void main(void)
{
    m_complex a = {1, 3};
    m_complex b = {2, 2};
    m_complex c = {0, 0};
    c = c_add(a, b);
    printf("%4.3lf+%4.3lfi\n", c.real, c.image);
    c = c_sub(a, b);
    printf("%4.3lf+%4.3lfi\n", c.real, c.image);
    c = c_mul(a,b);
    printf("%4.3lf+%4.3lfi\n", c.real, c.image);
    c = c_divi(a,b);
    printf("%4.3lf+%4.3lfi\n", c.real, c.image);
    printf("%4.3lf  %4.3lf\n", c_mod(a), c_mod(b));
    m_complex* A;
    A = (m_complex *)malloc(4*sizeof(m_complex));
    A[0].real =1;A[0].image =1;
    A[1].real =1;A[1].image =2;    
    A[2].real =2;A[2].image =1;
    A[3].real =2;A[3].image =3;
    c_mprintf(A, 2);
    c = c_det(A, 2);
    printf("%4.3lf+%4.3lfi\n", c.real, c.image);
    m_complex Ap[2] = {{1,2}, {4,6}};
    m_complex B[2] = {{2,1}, {3,1}};
    printf("%4.3lf+%4.3lfi\n", c_dot(Ap, B, 2).real, c_dot(Ap, B, 2).image);
    printf("%4.3lf+%4.3lfi\n", c_dot(Ap, Ap, 2).real, c_dot(Ap, Ap, 2).image);


}