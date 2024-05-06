#include "include/my_SVD.h"
#include "/home/heaven/Desktop/doc/Linear_Algebra/m_header_file/random_num/mt19937ar-master/mt19937ar.c"

const int ROW = 5;
const int COL = 2;

int main(void){
    double complex *A;
    A = (double complex*)malloc(ROW*COL*sizeof(double complex));
    double complex *U;
    U = (double complex*)malloc(ROW*ROW*sizeof(double complex));
    double complex *S;
    S = (double complex*)malloc(ROW*COL*sizeof(double complex));
    double complex *Vdag;
    Vdag = (double complex*)malloc(COL*COL*sizeof(double complex));

    init_genrand(time(NULL));
    // make A a rand matrix which numbers are within a circle with radius 10 in complex space
    for(int i=0; i<ROW*COL; i++){
        A[i] = genrand_real3()*10+genrand_real3()*10*I;
    }
    //test
    //A[0]=4;A[1]=11;A[2]=14;A[3]=8;A[4]=7;A[5]=-2; 
    printf("The matrix A is: \n");
    m_cprint(A,ROW,COL);

    // Singular Value Decomposition
    m_SVD(ROW, COL, A, U, Vdag, S);

    //check
    printf("is U and Vdag Unitary?\n");
    printf("U:\n");
    m_isUnitary(U,ROW);
    printf("Vdag:\n");
    m_isUnitary(Vdag,COL);
    printf("USVdag is (theorically A):\n");
    double complex a = 0+0*I;
    for(int i=0; i<ROW; i++){
        for(int j=0; j<COL; j++){
            a = 0+0*I;
            for(int k=0; k<ROW; k++){
                for(int l=0; l<COL; l++){
                    a += U[i*ROW+k]*S[k*COL+l]*Vdag[l*COL+j];
                }
            }
            printf("%3.4lf+i%3.4lf  ", creal(a), cimag(a));
        }printf("\n");
    }printf("\n");


    free(A);
    free(U);
    free(S);
    free(Vdag);
    return 0;
}