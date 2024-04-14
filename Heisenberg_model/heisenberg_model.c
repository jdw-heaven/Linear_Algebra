//spin-1/2 AF Heisenberg_model
#include "include/header_file.h"

const int L = 2;    //the size of sublattice
const double error = 1e-5;

int main(void)
{
    //the Pauli matrix
    double complex Sx[4] = {0,0.5,0.5,0};
    double complex Sy[4] = {0,-0.5*I,0.5*I,0};
    double complex Sz[4] = {0.5,0,0,-0.5};

    //the expectation value of spin, use random number within -1/2 to 1/2
    init_genrand(time(NULL));
    double complex* mean_spin;
    mean_spin = (double complex*)malloc(3*L*L*sizeof(double complex));
    for(int i = 0; i < 3*L*L; i++){
        mean_spin[i] = (genrand_real1()-0.5)+0*I; // they're all real numbers
    }
    printf("the initial spin is :\n");
    m_cprint(mean_spin, 3, L*L);

    //the Hamiltonian
    double complex *Ham;
    Ham = (double complex*)malloc((int)pow(2,2*L*L)*sizeof(double complex));
    //double complex *Hamp;
    //Hamp = (double complex*)malloc((int)pow(2,2*L*L)*sizeof(double complex));
    double *eigenvalues;
    eigenvalues = (double *)malloc((int)pow(2,L*L)*sizeof(double));
    if (eigenvalues == NULL) {
        // 内存分配失败，打印错误消息并退出
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    double complex *g_eigenvector;
    g_eigenvector = (double complex*)malloc((int)pow(2,L*L)*sizeof(double complex));
    double E_MAX = 0;
    
    do{
        E_MAX = 0;  //make sure each loop with a zero E_MAX
        m_Ham(L, mean_spin, Ham);   //calculate hamiltonian
        printf("the hamiltonian is :\n");
        m_cprint(Ham, (int)pow(2,L*L), (int)pow(2,L*L));
        m_isermie(Ham, (int)pow(2,L*L));

        /*
        for(int i=0; i<(int)pow(2,2*L*L); i++){
            Hamp[i] = Ham[i];
        }*/
        //the next step is to get the ground-energy and it's eigen vector
        // use LAPACKE_zheev to calculate eigenvalues and eigenvectors
        // 特征向量将按列存储在eigenvector中
        int info = LAPACKE_zheev(LAPACK_ROW_MAJOR, 'V', 'U', (int)pow(2,L*L), Ham, (int)pow(2,L*L), eigenvalues);
        if (info != 0) {
            fprintf(stderr, "LAPACKE_zheev returned info=%d\n", info);
            free(Ham);
            free(eigenvalues);
            return 1;
        }

        // 很好的一点是特征值已经按从小到大排列好了，第一个元素就是基态能量。下面的部分将给出验证特征向量是按列排布的：
        for(int i=0; i<(int)pow(2,L*L); i++){
            g_eigenvector[i] = Ham[i*(int)pow(2,L*L)+0];
        }


        printf("g_eigenvector is :\n");
        m_cprint(g_eigenvector, (int)pow(2,L*L), 1);
        //printf("test the eigenvector is range by row or col:\n");
        //m_cprint(m_mul(Hamp, g_eigenvector, (int)pow(2,L*L), (int)pow(2,L*L), 1), (int)pow(2,L*L), 1);
        

        //renew the mean_spin




        
        printf("the eigenvalues are: \n");
        m_print(eigenvalues, 1, (int)pow(2,L*L));
        printf("the corresponding eigenvectors are (range as column): \n");
        m_cprint(Ham, (int)pow(2,L*L), (int)pow(2,L*L));
        //printf("is eigenvectors unitary?\n");
        //m_isUnitary(eigenvectors, (int)pow(2,2*L));

        
    }while( E_MAX>error );

    free(Ham);
    free(eigenvalues);
    free(g_eigenvector);
    free(mean_spin);
    return 0;
}