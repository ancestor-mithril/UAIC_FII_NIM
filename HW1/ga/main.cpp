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
    return 0;
}
