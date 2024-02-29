#ifndef _M_COMPLEX_H_	
#define _M_COMPLEX_H_

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
//复数数据类型定义：
typedef struct {
    double real;
    double image;
} m_complex;
/*
复数的基本运算
*/
//复数的加法：
m_complex c_add(m_complex a, m_complex b){
    m_complex c = {0, 0};
    c.real = a.real + b.real;
    c.image = a.image + b.image;

    return c;
}

//复数的减法(a-b)：
m_complex c_sub(m_complex a, m_complex b){
    m_complex c = {0, 0};
    c.real = a.real - b.real;
    c.image = a.image - b.image;

    return c;
}

//复数的乘法：
m_complex c_mul(m_complex a, m_complex b){
    m_complex c = {0, 0};
    c.real = a.real*b.real - a.image*b.image;
    c.image = a.real*b.image + a.image*b.real;

    return c;
}

//复数的除法(c = a/b)：
m_complex c_divi(m_complex a, m_complex b){
    m_complex c = {0, 0};
    c.real = (a.real*b.real + a.image*b.image)/(b.real*b.real + b.image*b.image);
    c.image = (-a.real*b.image + a.image*b.real)/(b.real*b.real + b.image*b.image);

    return c;
}
//复数的共轭
m_complex c_conj(m_complex a){
    m_complex a_conj = {0, 0};
    a_conj.real = a.real;
    a_conj.image = -a.image;

    return a_conj;
}
//复数的模
double c_mod(m_complex a){
    return sqrt(a.real*a.real+a.image*a.image); 
}
/*
复空间内的向量代数
*/
//复矩阵的打印(n*n)：
void c_mprintf(m_complex* A, int n){
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            printf("%4.3lf+%4.3lfi   ", A[i*n+j].real, A[i*n+j].image);
        }printf("\n");
    }printf("\n");
}
//内积(n是向量的维数)<a,b>：
m_complex c_dot(m_complex* a, m_complex* b, int n){
    m_complex c = {0, 0};
    for(int i=0; i<n; i++){
        c = c_add(c, c_mul(c_conj(a[i]), b[i]));
    }
    return c;
}
//行列式（n是向量的维数,此程序不改变原矩阵）：
m_complex c_det(m_complex* A, int n){
    m_complex* Ap;
    Ap = (m_complex *)malloc(n*n*sizeof(m_complex));
    for(int i=0; i<n*n; i++){
        Ap[i] = A[i];
    }
    m_complex det = {0, 0};
    int raw = 0;
    double mod = 0;
    m_complex coi = {1, 0};
    m_complex med = {0, 0};
    m_complex medp = {0, 0};
    m_complex medpp = {0, 0};
    for(int i=0; i<n-1; i++){
        //将模最大的作为主元
        raw = i;
        for(int j=i+1; j<n; j++){
            mod = c_mod(Ap[i*n+i]);
            if( mod < c_mod(Ap[j*n+i])){
                mod = c_mod(Ap[j*n+i]);
                raw = j;
            }
        }
        if( raw!=i ){
            for(int j=i; j<n; j++){
                med = Ap[i*n+j];
                Ap[i*n+j] = Ap[raw*n+j];
                Ap[raw*n+j] = med;
            }
            coi.real = (-1)*coi.real;
        }
        //化为三角矩阵
        for(int j=i+1; j<n; j++){
            med = Ap[j*n+i];
            for(int k=i; k<n; k++){
                medp = c_divi(med, Ap[i*n+i]);
                medpp = c_mul(medp, Ap[i*n+k]);
                Ap[j*n+k] = c_sub(Ap[j*n+k], medpp);
            }
        }
    }
    //c_mprintf(Ap, n);
    //系数+对角元之积
    det = c_mul(coi, Ap[0]);
    for(int i=1; i<n; i++){
        det = c_mul(det, Ap[i*n+i]);
    }

    return det;
}

#endif