#include <iostream>
#include <cmath> // 或者 #include <math.h>

int main() {
    double num = 1.5;
    int rounded = static_cast<int>(std::round(num));
    std::cout << rounded << std::endl; // 输出 2
    return 0;
}