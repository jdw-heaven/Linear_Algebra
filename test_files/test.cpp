/*
spin-1/2 XXY 反铁磁链基态能量求解
*/
#include <iostream>
#include <complex>
#include <torch/torch.h>
#include <cmath>

const int L = 2;   //单周期内的粒子数
const double error = 1e-12; //数值零

// 将 a 的第 b 位取反 ，最低位编号为 0
int flapBit(int a, int b) { return a ^ (1 << b); }
// 获取 a 的第 b 位，最低位编号为 0
int getBit(int a, int b) { return (a >> b) & 1; }
// H作用在f上
void m_Hf(torch::Tensor& v, torch::Tensor& f, double Delta){
    // v一定要是空向量
    int num = 0;
    for(int i=0; i<L; i++){
        for(int j=0; j<(1<<L); j++){
            if(getBit(j,i)==getBit(j,(i+1)%L)){
                v[j] += Delta*f[j]/4;
            }else{
                v[j] -= Delta*f[j]/4;
                num = flapBit(flapBit(j,(i+1)%L),i);
                //std::cout << num << std::endl;
                std::cout << f[num]/2 << std::endl;
                v[num] += f[num]/2;
            }    
        }
    }
}


int main(void)
{
    torch::Tensor fr = torch::tensor({1,2,3,4}, torch::kDouble);
    torch::Tensor fi = torch::zeros({1,2,3,4}, torch::kDouble);
    torch::Tensor f = torch::complex(fr,fi);

    auto v = torch::norm(f);


    std::cout << v << std::endl;


    return 0;
}
