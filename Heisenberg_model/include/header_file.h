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

//is unitary?
void m_isUnitary(double complex *A, int n){
    double complex *U;
    U = (double complex *)malloc(n*n*sizeof(double complex));
    for(int i = 0; i < n*n; i++){
        U[i] = 0+0*I;
    }
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            for(int k = 0; k < n; k++){
                U[i*n+j] += conj(A[k*n+i])*A[k*n+j];
            }
        }
    }
    m_cprint(U, n, n);
    free(U);
}

//is ermie?
void m_isermie(double complex *A, int n){
    int isermie = 1;
    for(int i=0; i<n; i++){
        if( fabs(cimag(A[i*n+i])) > 1e-10 ){
            isermie = 0;
        }else{
            ;
        }
    }
    for(int i=0; i<n-1; i++){
        for(int j=1; j<n; j++){
            if( fabs(cimag(A[i*n+j])+cimag(A[j*n+i])) > 1e-10 ){
                isermie = 0;
            }else{
                ;
            }
        }
    }
    if( isermie == 0 ){
        printf("Not ermie.\n");
    }else{
        printf("Ermie.\n");
    }
}

//multiply
double complex * m_mul(double complex *A, double complex *B, int Arow, int Acolumn, int Bcolumn){
    double complex *C;
    C = (double complex *)malloc(Arow*Bcolumn*sizeof(double complex));
    for(int i = 0; i < Arow*Bcolumn; i++){
        C[i] = 0+0*I;
    }
    for(int i = 0; i < Arow; i++){
        for(int j = 0; j < Bcolumn; j++){
            for(int k = 0; k < Acolumn; k++){
                C[i*Bcolumn+j] += A[i*Acolumn+k]*B[k*Bcolumn+j];
            }
        }
    }
    return C;
    free(C);
}

//the Kronecker Product
void m_KP(double complex *A, double complex *B, int Arow, int Acolumn, int Brow, int Bcolumn, double complex *C){
    for(int i=0; i<Arow; i++){
        for(int j=0; j<Acolumn; j++){
            for(int k=0; k<Brow; k++){
                for(int l=0; l<Bcolumn; l++){
                    C[(i*Brow+k)*Acolumn*Bcolumn+(j*Bcolumn+l)] = A[i*Acolumn+j]*B[k*Bcolumn+l];
                }
            }
        }
    }
}

// find the minimum eigenvalue
int m_fmin(double complex *eigenvalues, int n){
    int num = 0;
    for(int i=1; i<n; i++){
        if( creal(eigenvalues[i-1])<creal(eigenvalues[i]) ){
            num = i-1;
        }
    }
    return num;
}

//function to trans state and number
int state_to_num(int *state, int L){
    int num = 0;
    for(int i=0; i < L*L; i++){
        num += state[i]*(int)pow(2,L*L-i-1);
    }

    return num;
}
void num_to_state(int *state, int num, int L){
    int a = 0;
    for(int i=0; i<L*L; i++){
        a = num%(int)pow(2,L*L-i);
        state[i] = (int)(a/(i==L*L-1? 1:pow(2,L*L-1-i)));
    }
}

//the Hamiltonian matrix
void m_Ham(const int L, double complex *the_mean_spin, double complex *the_hamiltonian){
    for(int i=0; i<(int)pow(2,2*L*L); i++){ the_hamiltonian[i] = 0+0*I; }
    //a matrix with 2^the_sublattice_size to stand for the state
    int *state;
    state = (int*)malloc(L*L*sizeof(int));
    // 检查内存分配是否成功
    if (state == NULL) {
        // 内存分配失败，打印错误消息并退出
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    int numj = 0, coe = 1;
    //H_ij
    for(int i=0; i<(int)pow(2,L*L); i++){
        for(int j=0; j<(int)pow(2,L*L); j++){
            //the row nearest particles' coupling
            for(int k=0; k<L; k++){
                for(int l=0; l<L-1; l++){ //each row has the_sublattice_size-1 bonds
                    //S_k*size+l and S_k*size+l+1
                    //trans j to state
                    num_to_state(state, j, L);
                    //Sx and Sy
                    coe = 1;
                    if(state[k*L+l+1]==1){
                        coe = coe*(-1);
                    }else{
                        ;
                    }
                    state[k*L+l+1] = state[k*L+l+1]-1==0? 0:1;
                    if(state[k*L+l]==1){
                        coe = coe*(-1);
                    }else{
                        ;
                    }
                    state[k*L+l] = state[k*L+l]-1==0? 0:1;
    
                    //trans to num
                    numj = state_to_num(state, L);
                    if(i==numj){
                        the_hamiltonian[i*(int)pow(2,L*L)+j] += 0.25-coe*0.25;
                    }else{
                        ;
                    }
                    //Sz
                    if(i==j){
                        the_hamiltonian[i*(int)pow(2,L*L)+j] += coe*0.25;
                    }
                }
            }

            //the column nearest particles' coupling
            for(int k=0; k<L; k++){
                for(int l=0; l<L-1; l++){
                    //S_l*size+k and S_(l+1)*size+k
                    //trans j to state
                    num_to_state(state, j, L);
                    //Sx and Sy
                    coe = 1;
                    if(state[l*L+k]==1){
                        coe = coe*(-1);
                    }else{
                        ;
                    }
                    if(state[(l+1)*L+k]==1){
                        coe = coe*(-1);
                    }else{
                        ;
                    }
                    state[l*L+k] = state[l*L+k]-1==0? 0:1;
                    state[(l+1)*L+k] = state[(l+1)*L+k]-1==0? 0:1;
                    //trans to num
                    numj = state_to_num(state, L);
                    if(i==numj){
                        the_hamiltonian[i*(int)pow(2,L*L)+j] += 0.25-coe*0.25;
                    }else{
                        ;
                    }
                    //Sz
                    if(i==j){
                        the_hamiltonian[i*(int)pow(2,L*L)+j] += coe*0.25;
                    }
                }
            }

            // the up-outside of the sublattice, particle 0 to L-1
            for(int k=0; k<L; k++){
                //the particle upon them is equal to (L-1)*L ot L*L-1
                //trans j to state
                num_to_state(state, j, L);
                //Sx and Sy
                coe = 1;
                if(state[0*L+k]==1){
                    coe = coe*(-1);
                }else{
                    ;
                }
                state[0*L+k] = state[0*L+k]-1==0? 0:1;
                //trans to num
                numj = state_to_num(state, L);
                if(i==numj){
                    the_hamiltonian[i*(int)pow(2,L*L)+j] += 0.5*the_mean_spin[0*L*L+(L-1)*L+k]+coe*0.5*the_mean_spin[1*L*L+(L-1)*L+k]*I;
                }else{
                    ;
                }
                //Sz
                if(i==j){
                    the_hamiltonian[i*(int)pow(2,L*L)+j] += coe*0.5*the_mean_spin[2*L*L+(L-1)*L+k];
                }
            }

            //the down-outside of the sublattice, particle (L-1)*L ot L*L-1
            for(int k=0; k<L; k++){
                //the particle upon them is equal to 0 to L-1
                //trans j to state
                num_to_state(state, j, L);
                //Sx and Sy
                coe = 1;
                if(state[(L-1)*L+k]==1){
                    coe = coe*(-1);
                }else{
                    ;
                }
                state[(L-1)*L+k] = state[(L-1)*L+k]-1==0? 0:1;
                //trans to num
                numj = state_to_num(state, L);
                if(i==numj){
                    the_hamiltonian[i*(int)pow(2,L*L)+j] += 0.5*the_mean_spin[0*L*L+k]+coe*0.5*the_mean_spin[1*L*L+k]*I;
                }else{
                    ;
                }
                //Sz
                if(i==j){
                    the_hamiltonian[i*(int)pow(2,L*L)+j] += coe*0.5*the_mean_spin[2*L*L+k];
                }
            }

            //the left-outside of the sublattice, paticle 0*L to (L-1)*L
            for(int k=0; k<L; k++){
                //the particle upon them is equal to 0*L+L-1 to (L-1)*L+L-1
                //trans j to state
                num_to_state(state, j, L);
                //Sx and Sy
                coe = 1;
                if(state[k*L]==1){
                    coe = coe*(-1);
                }else{
                    ;
                }
                state[k*L] = state[k*L]-1==0? 0:1;
                //trans to num
                numj = state_to_num(state, L);
                if(i==numj){
                    the_hamiltonian[i*(int)pow(2,L*L)+j] += 0.5*the_mean_spin[0*L*L+k*L+L-1]+coe*0.5*the_mean_spin[1*L*L+k*L+L-1]*I;
                }else{
                    ;
                }
                //Sz
                if(i==j){
                    the_hamiltonian[i*(int)pow(2,L*L)+j] += coe*0.5*the_mean_spin[2*L*L+k*L+L-1];
                }
            }

            //the right-outside of the sublattice, particle 0*L+L-1 to (L-1)*L+L-1
            for(int k=0; k<L; k++){
                //the particle upon them is equal to 0*L to (L-1)*L
                //trans j to state
                num_to_state(state, j, L);
                //Sx and Sy
                coe = 1;
                if(state[k*L+L-1]==1){
                    coe = coe*(-1);
                }else{
                    ;
                }
                state[k*L+L-1] = state[k*L+L-1]-1==0? 0:1;
                //trans to num
                numj = state_to_num(state, L);
                if(i==numj){
                    the_hamiltonian[i*(int)pow(2,L*L)+j] += 0.5*the_mean_spin[0*L*L+k*L]+coe*0.5*the_mean_spin[1*L*L+k*L]*I;
                }else{
                    ;
                }
                //Sz
                if(i==j){
                    the_hamiltonian[i*(int)pow(2,L*L)+j] += coe*0.5*the_mean_spin[2*L*L+k*L];
                }
            }

            
        }
    }
    free(state);  
}

//renew the mean_spin and give the MAX error
double ren(double complex *mean_spin, double complex *g_eigenvetor, const int L){
    double E_MAX = 0;
    double complex x = 0+0*I;
    double complex y = 0+0*I;
    double complex z = 0+0*I;
    int num = 0;
    int coe = 1;
    int *state;
    state = (int*)malloc(L*L*sizeof(int));
    //遍历粒子
    for(int k=0; k<L*L; k++){
        x = 0+0*I;
        y = 0+0*I;
        z = 0+0*I;
        // 遍历态，系数在g_eigenvector中
        for(int i=0; i<(int)pow(2,L*L); i++){
            for(int j=0; j<(int)pow(2,L*L); j++){

                //Sx and Sy
                num_to_state(state, j, L);
                coe = 1;
                if(state[k]==1){
                    coe = coe*(-1);
                }else{;}
                state[k] = state[k]-1==0? 0:1;
                num = state_to_num(state, L);
                if(i==num){
                    x += conj(g_eigenvetor[i])*g_eigenvetor[j]*0.5;
                    y += conj(g_eigenvetor[i])*g_eigenvetor[j]*coe*0.5*I;
                }

                //Sz
                if(i==j){
                    z += conj(g_eigenvetor[i])*g_eigenvetor[j]*coe*0.5;
                }
            }
        }
        // find max error
        if(fabs( fabs(creal(mean_spin[0*L*L+k])) - fabs(creal(x)) )>E_MAX) E_MAX = fabs( fabs(creal(mean_spin[0*L*L+k])) - fabs(creal(x)) );
        if(fabs( fabs(creal(mean_spin[1*L*L+k])) - fabs(creal(y)) )>E_MAX) E_MAX = fabs( fabs(creal(mean_spin[1*L*L+k])) - fabs(creal(y)) );
        if(fabs( fabs(creal(mean_spin[2*L*L+k])) - fabs(creal(z)) )>E_MAX) E_MAX = fabs( fabs(creal(mean_spin[2*L*L+k])) - fabs(creal(z)) );
        //renew mean_spin
        mean_spin[0*L*L+k] = x;
        mean_spin[1*L*L+k] = y;
        mean_spin[2*L*L+k] = z;

    }
    free(state);
    return E_MAX;
}

//the magnet moment
void magnet_moment(double sx, double sy, double sz){
    double t = sx*sx+sy*sy+sz*sz;
    printf("the magnet moment M is:\n");
    printf("%3.4lf\n", sqrt(t));
}


#endif