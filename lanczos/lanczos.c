/*
spin-1/2 XXY 反铁磁链基态能量求解
*/
#include "include/lanczos.h"

int main(void)
{
    FILE *data;
    data = fopen("draw/data/data12_13_2.txt", "w");
    init_genrand(time(NULL));

    for(double Delta = -2.0; Delta <2.01; Delta += 0.01){
        fprintf(data, "%lf %lf\n", Delta, m_lanczos(Delta));
    }


    // double Delta = -1.0;
    // printf("%lf\n", m_lanczos(Delta));
    fclose(data);
    return 0;
}