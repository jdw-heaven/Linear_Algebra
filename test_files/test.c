#include <iostream>
#include <cmath>
#include <cstdlib>
#include <complex>
#include <torch/torch.h>

int main(void)
{
    /*the number that recognizes as zero*/
    const double jz = 1e-25;
    // 参数
    const int L = 4;    // 粒子数，本程序仅适用于偶数
    const double s = 1.0;   //tips: not 1/2, but 1.0/2
    const int D = 4;    // hilbert子空间维数
    const double tau = 0.01;    // 单位虚时

    // 自旋算符的矩阵表示
    /*单粒子自旋空间的维数*/
    const int d = static_cast<int>(std::round(2*s+1));
    // std::cout << d << std::endl;
    torch::Tensor Sz_real = torch::zeros({d,d}, torch::kDouble);
    torch::Tensor Sz_imag = torch::zeros({d,d}, torch::kDouble);
    torch::Tensor Sr_real = torch::zeros({d,d}, torch::kDouble);
    torch::Tensor Sr_imag = torch::zeros({d,d}, torch::kDouble);
    torch::Tensor Sl_real = torch::zeros({d,d}, torch::kDouble);
    torch::Tensor Sl_imag = torch::zeros({d,d}, torch::kDouble);

    /*自旋最大的为第0态*/
    double s_temp = s;
    for(int i = 0; i < d; i++){
        Sz_real[i][i] = s_temp;
        /*注意，对于Sr来说，s_temp需要-1*/
        if(i+1<d){
            Sr_real[i][i+1] = std::sqrt((s-s_temp+1)*(s+s_temp));
            Sl_real[i+1][i] = std::sqrt((s+s_temp)*(s-s_temp+1));
        }else{
            ;
        }
        s_temp -= 1;
    }
    // std::cout << Sz_real << std::endl;
    // std::cout << Sr_real << std::endl;
    // std::cout << Sl_real << std::endl;
    /*合并实虚部*/
    torch::Tensor Sz = torch::complex(Sz_real,Sz_imag);
    torch::Tensor Sr = torch::complex(Sr_real,Sr_imag);
    torch::Tensor Sl = torch::complex(Sl_real,Sl_imag);
    // std::cout << Sz.dtype() << std::endl;
    
    // h_i的矩阵表示
    torch::Tensor h = torch::kron(Sz,Sz) + 1.0/2*(torch::kron(Sr,Sl)+torch::kron(Sl,Sr));
    // std::cout << real(h) << std::endl;
    // std::cout << imag(h) << std::endl;

    // 虚时演化算符
    torch::Tensor ite = torch::zeros({d*d,d*d}, torch::kComplexDouble);
    for(int i=0; i<d*d; i++){
        for(int j=0; j<d*d; j++){
            if(std::abs(real(h[i][j]).item<double>())>jz){
                real(ite[i][j]) = torch::exp(-tau*real(h[i][j]));
            }else{
                ;
            }
        }
    }
    // std::cout << real(ite) << std::endl;
    // std::cout << imag(ite) << std::endl;

    // 随机生成周期性MPS
    std::vector<torch::Tensor> A;
    std::vector<torch::Tensor> Gamma;
    for(int i=0; i<L; i++){
        torch::Tensor a_r = torch::randn({D,d,D}, torch::kDouble);
        torch::Tensor a_i = torch::randn({D,d,D}, torch::kDouble);
        torch::Tensor a = torch::complex(a_r,a_i);
        A.push_back(a);

        torch::Tensor gamma = torch::randn({D}, torch::kDouble);
        Gamma.push_back(gamma);
    }

    // 进行虚时演化迭代，奇偶交替
    int times = 0;
    double g_E = -1.0;  //随便给个值
    double energy_error = 0.0;
    int j = 0;
    /*中转矩阵U、S、VH*/
    torch::Tensor U = torch::zeros({D*d,D*d}, torch::kComplexDouble);
    torch::Tensor S = torch::zeros({D*d}, torch::kDouble);
    torch::Tensor VH = torch::zeros({D*d,D*d}, torch::kComplexDouble);
    torch::Tensor M = torch::zeros({D*d,D*d}, torch::kComplexDouble);
    
    do{
        times += 1;
        /*奇数*/
        for(int i=0; i<L; i++,i++){
            j = i-1;
            if(i==0) j = L-1;
            auto K = torch::einsum("a,abj,j,jcd,d->abcd",{Gamma[j],A[i],Gamma[i],A[i+1],Gamma[i+1]});
            M = torch::einsum("ab,cbd->cad",{ite,K.reshape({D,d*d,D})}).reshape({D*d,d*D});
            std::tie(U, S, VH) = torch::linalg_svd(M, true);

            A[i] = U.slice(1,0,D).clone().reshape({D,d,D});
            A[i] = torch::einsum("a,abc->abd",{torch::reciprocal(Gamma[j]),A[i]});
            Gamma[i] = S.slice(0,0,D).clone();
            A[i+1] = VH.slice(0,0,D).clone().reshape({D,d,D});
            A[i+1] = torch::einsum("abc,c->abc",{A[i+1],torch::reciprocal(Gamma[i+1])});
            
            torch::Tensor E = torch::norm(S,2);
            Gamma[i] = Gamma[i]/E.item<double>();
            energy_error = std::abs(g_E-std::log(E.item<double>())/(-tau));
            g_E = std::log(E.item<double>())/(-tau);
        }
        /*偶数*/
        for(int i=1; i<L; i++,i++){
            auto K = torch::einsum("a,abj,j,jcd,d->abcd",{Gamma[i-1],A[i],Gamma[i],A[(i+1)%L],Gamma[(i+1)%L]});
            M = torch::einsum("ab,cbd->cad",{ite,K.reshape({D,d*d,D})}).reshape({D*d,d*D});
            std::tie(U, S, VH) = torch::linalg_svd(M, true);

            A[i] = U.slice(1,0,D).clone().reshape({D,d,D});
            A[i] = torch::einsum("a,abc->abd",{torch::reciprocal(Gamma[i-1]),A[i]});
            Gamma[i] = S.slice(0,0,D).clone();
            A[(i+1)%L] = VH.slice(0,0,D).clone().reshape({D,d,D});
            A[(i+1)%L] = torch::einsum("abc,c->abc",{A[(i+1)%L],torch::reciprocal(Gamma[(i+1)%L])});


            torch::Tensor E = torch::norm(S,2);
            Gamma[i] = Gamma[i]/E.item<double>();
            energy_error = std::abs(g_E-std::log(E.item<double>())/(-tau));
            g_E = std::log(E.item<double>())/(-tau);
        }
    }while(energy_error > jz);

    std::cout << "times is :" << times << std::endl;
    std::cout << "the ground energy is:" << g_E << std::endl;

    return 0;
}