/*
spin-1/2 XXY 反铁磁链基态能量求解
*/
#include <iostream>
#include <complex>
#include <torch/torch.h>
#include <cmath>

const int L = 12;   //单周期内的粒子数
const double error = 1e-12; //数值零

// 将 a 的第 b 位取反 ，最低位编号为 0
int flapBit(int a, int b) { return a ^ (1 << b); }
// 获取 a 的第 b 位，最低位编号为 0
int getBit(int a, int b) { return (a >> b) & 1; }
// H作用在f上
void m_Hf(torch::Tensor& v, torch::Tensor& f, double Delta){
    // v一定要是空向量
    int num = 0;
    for(int j=0; j<(1<<L); j++){
        for(int i=0; i<L; i++){
            if(getBit(j,i)==getBit(j,(i+1)%L)){
                v[j] += Delta*f[j]/4;
            }else{
                v[j] -= Delta*f[j]/4;
                num = flapBit(flapBit(j,(i+1)%L),i);
                v[num] += f[num]/2;
            }    
        }
    }
}
//lanczos iteration
double m_lanczos(double Delta){
    torch::Tensor alpha = torch::zeros({1<<L}, torch::kDouble);
    torch::Tensor beta = torch::zeros({1<<L}, torch::kDouble);

    torch::Tensor fr = torch::zeros(1<<L, torch::kDouble);
    torch::Tensor fi = torch::zeros(1<<L, torch::kDouble);
    torch::Tensor f = torch::complex(fr,fi);
    torch::Tensor v = torch::complex(fr,fi);
    torch::Tensor fp = torch::complex(fr,fi);
    
    f = torch::randn(1<<L, torch::kComplexDouble);
    f = f/torch::norm(f);

    int n = 0;
    do{
        n += 1;
        m_Hf(v,f,Delta);
        v = v - beta[n-1]*fp;
        alpha[n] = torch::einsum("i,i->", {torch::conj(f),v});
        v = v - alpha[n]*f;  //增加数值稳定性
        beta[n] = torch::norm(v);
        fp = f;
        f = v/beta[n];
        v = torch::complex(fr,fi);
    }while( beta[n].item<double>() > error );

    //求解本征值
    torch::Tensor L = torch::zeros({n,n}, torch::kDouble);
    for(int i=0; i<n; i++){
        L[i][i] = alpha[i+1];
    }
    for(int i=0; i<n-1; i++){
        L[i][i+1] = beta[i+1];
        L[i+1][i] = beta[i+1];
    }

    torch::Tensor e = torch::linalg_eigvalsh(L);
    //std::cout << e <<std::endl;
    double g_E = e[0].item<double>();

    return g_E;
}

int main(void)
{
    /*
    for(double Delta = -1.5; Delta <=1.5; Delta += 0.05){
        std::cout << m_lanczos(Delta) << std::endl;
    }
    */
    
    double Delta = 1.0;
    std::cout << m_lanczos(Delta) << std::endl;
    return 0;
}