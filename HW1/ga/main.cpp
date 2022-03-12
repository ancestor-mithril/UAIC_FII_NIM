#include <iostream>
#include <vector>
import GA.Cec22Specific;
// clang-format off
// import <iostream>;
// import <vector>;
// clang-format on

int main() {
    ga::functions::sanity_check();
    std::vector<double> a = {1.1, 1.2};
    std::vector<double> b = {0.1, 0.2};
    std::vector<std::vector<double>> c = {{1.1, 1.2}, {1.1, 1.2}};

    std::cout << 'a' << " = " << ga::functions::rastrigin_func(a, b, c, true, true)<< '\n';

    std::vector<double> x = {1.1, 1.2};
    std::vector<double> shift = {0.1, 0.2, 0.11, 0.21,0.12, 0.22};
    std::vector<std::vector<double>> rotate = {{1.1, 1.2}, {1.1, 1.2}, {1.1, 1.2}, {1.1, 1.2}, {1.1, 1.2}, {1.1, 1.2}};
    std::cout << "b = " << ga::functions::cf02(x, shift, rotate, true) << '\n';
    return 0;
}
