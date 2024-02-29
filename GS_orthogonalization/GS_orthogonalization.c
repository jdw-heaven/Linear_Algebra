#include"../m_header_file/m_complex.h"
#include"../m_header_file/random_num/mt19937ar-master/mt19937ar.c"

//归一化函数
void normalize(m_complex *a, int n){
    double a_mod = sqrt(c_mod(c_dot(a, a, n)));
    for(int i=0; i<n; i++){
        a[i].real = a[i].real/a_mod;
        a[i].image = a[i].image/a_mod;
    }

}
//主函数
void m_GS_Complex(void)
{
    //首先给出线性空间的维数
    int n = 0;
    printf("Please enter the demension of linear space:\n");
    scanf("%d", &n);
    //n个线性无关向量的生成
    m_complex x = {0, 0};//judge whether the arrow is linear independant
    m_complex* A;
    A = (m_complex *)malloc(n*n*sizeof(m_complex));
    do{
        for(int i=0; i<n*n; i++){
            A[i].real = 3*genrand_real1();//genrand_real1()生成[0,1]的随机数
            A[i].image = 3*genrand_real1();
        }
        x = c_det(A, n);
    }while(!c_mod(x));
    printf("//打印A看看\n");
    c_mprintf(A, n);//打印A看看
    //下面的代码可能有点难懂，我的目的是得到一个含有n个向量，每个向量都是n维复空间中的。
    struct m_complex_matrix{
        m_complex* arrow;
    };
    struct m_complex_matrix* m_complex_basis;
    m_complex_basis = (m_complex *)malloc(n*n*sizeof(m_complex));
    for(int i=0; i<n; i++){
        m_complex_basis[i].arrow = (m_complex *)malloc(n*sizeof(m_complex));
    }
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            m_complex_basis[i].arrow[j] = A[i*n+j];
        }
    }
    //printf("%lf",m_complex_basis[0].arrow[1].real);
    //正交化
    m_complex* Y;
    Y = (m_complex *)malloc(n*sizeof(m_complex));
    for(int i=1; i<n; i++){
        for(int l=0; l<n; l++){
            Y[l].real = 0;
            Y[l].image = 0;
        }
        for(int j=i-1; j>=0; j--){
            x = c_divi(c_dot(m_complex_basis[j].arrow, m_complex_basis[i].arrow, n), c_dot(m_complex_basis[j].arrow, m_complex_basis[j].arrow, n));
            for(int k=0; k<n; k++){
                Y[k] = c_add(Y[k], c_mul(x, m_complex_basis[j].arrow[k]));
            }
        }
        for(int k=0; k<n; k++){
            m_complex_basis[i].arrow[k] = c_sub(m_complex_basis[i].arrow[k], Y[k]);
        } 
    }
    //归一化
    for(int i=0; i<n; i++){
        normalize(m_complex_basis[i].arrow, n);
    }
    //按行打印基
    printf("//按行打印基\n");
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            printf("%4.3lf+%4.3lfi ", m_complex_basis[i].arrow[j].real, m_complex_basis[i].arrow[j].image);
        }printf("\n\n");
    }
    //检验基的正交归一性
    printf("检验基的正交归一性\n");
    //归一性
    for(int i=0; i<n; i++){
        printf("%lf   ", sqrt(c_mod(c_dot(m_complex_basis[i].arrow, m_complex_basis[i].arrow, n))));
    }
    printf("\n");
    //正交性
    for(int i=0; i<n; i++){
        for(int j=i+1; j<n; j++){
            printf("%lf   ", sqrt(c_mod(c_dot(m_complex_basis[i].arrow, m_complex_basis[j].arrow, n))));
        }
    }
    printf("\n");
    //坐标向量
    //任意整个向量
    printf("向量Y为:\n");
    for(int i=0; i<n; i++){
        Y[i].real = 3*genrand_real1();
        Y[i].image = 3*genrand_real1();
        printf("%4.3lf+%4.3lfi  ", Y[i].real, Y[i].image);
    }printf("\n");
    printf("Y的坐标为:\n");
    for(int i=0; i<n; i++){
        printf("%4.3lf+%4.3lfi\n", c_dot(Y, m_complex_basis[i].arrow, n).real, c_dot(Y, m_complex_basis[i].arrow, n).image);
    }
}
//实数也就是把上面的虚部整成0就可以了
void m_GS_Real(void)
{
    //首先给出线性空间的维数
    int n = 0;
    printf("Please enter the demension of linear space:\n");
    scanf("%d", &n);
    //n个线性无关向量的生成
    m_complex x = {0, 0};//judge whether the arrow is linear independant
    m_complex* A;
    A = (m_complex *)malloc(n*n*sizeof(m_complex));
    do{
        for(int i=0; i<n*n; i++){
            A[i].real = 3*genrand_real1();//genrand_real1()生成[0,1]的随机数
            A[i].image = 0;
        }
        x = c_det(A, n);
    }while(!c_mod(x));
    printf("//打印A看看\n");
    c_mprintf(A, n);//打印A看看
    //下面的代码可能有点难懂，我的目的是得到一个含有n个向量，每个向量都是n维复空间中的。
    struct m_complex_matrix{
        m_complex* arrow;
    };
    struct m_complex_matrix* m_complex_basis;
    m_complex_basis = (m_complex *)malloc(n*n*sizeof(m_complex));
    for(int i=0; i<n; i++){
        m_complex_basis[i].arrow = (m_complex *)malloc(n*sizeof(m_complex));
    }
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            m_complex_basis[i].arrow[j] = A[i*n+j];
        }
    }
    //printf("%lf",m_complex_basis[0].arrow[1].real);
    //正交化
    m_complex* Y;
    Y = (m_complex *)malloc(n*sizeof(m_complex));
    for(int i=1; i<n; i++){
        for(int l=0; l<n; l++){
            Y[l].real = 0;
            Y[l].image = 0;
        }
        for(int j=i-1; j>=0; j--){
            x = c_divi(c_dot(m_complex_basis[j].arrow, m_complex_basis[i].arrow, n), c_dot(m_complex_basis[j].arrow, m_complex_basis[j].arrow, n));
            for(int k=0; k<n; k++){
                Y[k] = c_add(Y[k], c_mul(x, m_complex_basis[j].arrow[k]));
            }
        }
        for(int k=0; k<n; k++){
            m_complex_basis[i].arrow[k] = c_sub(m_complex_basis[i].arrow[k], Y[k]);
        } 
    }
    //归一化
    for(int i=0; i<n; i++){
        normalize(m_complex_basis[i].arrow, n);
    }
    //按行打印基
    printf("//按行打印基\n");
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            printf("%4.3lf+%4.3lfi ", m_complex_basis[i].arrow[j].real, m_complex_basis[i].arrow[j].image);
        }printf("\n\n");
    }
    //检验基的正交归一性
    printf("检验基的正交归一性\n");
    //归一性
    for(int i=0; i<n; i++){
        printf("%lf   ", sqrt(c_mod(c_dot(m_complex_basis[i].arrow, m_complex_basis[i].arrow, n))));
    }
    printf("\n");
    //正交性
    for(int i=0; i<n; i++){
        for(int j=i+1; j<n; j++){
            printf("%lf   ", sqrt(c_mod(c_dot(m_complex_basis[i].arrow, m_complex_basis[j].arrow, n))));
        }
    }
    printf("\n");
    //坐标向量
    //任意整个向量
    printf("向量Y为:\n");
    for(int i=0; i<n; i++){
        Y[i].real = 3*genrand_real1();
        Y[i].image = 0;
        printf("%4.3lf+%4.3lfi  ", Y[i].real, Y[i].image);
    }printf("\n");
    printf("Y的坐标为:\n");
    for(int i=0; i<n; i++){
        printf("%4.3lf+%4.3lfi\n", c_dot(Y, m_complex_basis[i].arrow, n).real, c_dot(Y, m_complex_basis[i].arrow, n).image);
    }
}
void main(void)
{
    /* initializes mt[N] with a seed */
    init_genrand(time(NULL));
    m_GS_Complex();
    m_GS_Real();
}