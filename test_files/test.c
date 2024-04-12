#include "../Heisenberg_model/include/header_file.h"
const int L = 2;    //the size of sublattice
 double complex Sx[4] = {0,0.5,0.5,0}; //the Pauli matrix
 double complex Sy[4] = {0,-0.5*I,0.5*I,0};
 double complex Sz[4] = {0.5,0,0,-0.5};


int main(void)
{
    double complex Sx[4] = {0,0.5,0.5,0}; //the Pauli matrix
    double complex Sy[4] = {0,-0.5*I,0.5*I,0};
    double complex Sz[4] = {0.5,0,0,-0.5};
    double complex SS[16] = {0.25,0,0,0,
                            0,-0.25,0.5,0,
                            0,0.5,-0.25,0,
                            0,0,0,0.25};

    m_cprint(Sx, 2, 2);
    m_cprint(Sy, 2, 2);
    m_cprint(Sz, 2, 2);
    double complex *S;
    S = (double complex*)malloc(4*4*sizeof(double complex));
    m_cprint(SS, 4, 4);
    m_KP(Sz, Sz, 2, 2, 2, 2, S);
    m_cprint(S, 4, 4);
    int x = 9;
    printf("%d\n", x%2);
    printf("%d\n", x);

    int state[4] = {1,1,1,1};
    int num = 10;
    printf("%d\n", state_to_num(state, 2));
    num_to_state(state, num, 2);
    printf("%d\n", state[0]);
    printf("%d\n", state[1]);
    printf("%d\n", state[2]);
    printf("%d\n", state[3]);
    
    int *state1;
    state1 = (int*)malloc(2*2*sizeof(int));
    printf("%d\n", (int)pow(2,2*2));

    free(state1);
    free(S);
    return 0;
}


