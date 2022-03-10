module;
#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <numeric>
module GA.Cec22Specific;
// clang-format off
// import <iostream>;
// import <algorithm>;
// import <functional>;
// import <cmath>;
// import <numeric>;
// clang-format on

namespace ga::functions {
namespace {

constexpr double PI = 3.1415926535897932384626433832795029;
constexpr double E = 2.7182818284590452353602874713526625;

void shiftfunc(std::vector<double>& x, const std::vector<double>& shift)
{
    // assuming x.size() == shift.size()
    std::transform(x.begin(), x.end(), shift.begin(), x.begin(),
                   std::minus<double>());
}

std::vector<double> rotatefunc(std::vector<double>& x,
                               const std::vector<std::vector<double>>& rotate)
{
    // assuming x.size() == rotate.size() == rotate[0].size()
    std::vector<double> aux(x.size(), 0.0);
    for (auto i = 0; i < x.size(); ++i) {
        for (auto j = 0; j < x.size(); ++j) {
            aux[i] += x[j] * rotate[i][j];
        }
    }
    return aux;
}

std::vector<double>
shift_rotate_transform(std::vector<double>& x, const std::vector<double>& shift,
                       const std::vector<std::vector<double>>& rotate,
                       double shift_rate, bool shift_flag, bool rotate_flag)
{
    if (shift_flag) {
        shiftfunc(x, shift);
    }
    std::transform(x.begin(), x.end(), x.begin(),
                   std::bind(std::multiplies<double>(), std::placeholders::_1,
                             shift_rate));
    if (rotate_flag) {
        return rotatefunc(x, rotate);
    }
    return x;
    // Kept here for historical reason
    //
    // if (shift_flag && rotate_flag) {
    //     shiftfunc(x, shift);
    //     std::transform(x.begin(), x.end(), x.begin(),
    //     std::bind(std::multiplies<double>(), std::placeholders::_1,
    //     shift_rate)); return rotatefunc(x, rotate);
    // }
    // if (shift_flag) {
    //     shiftfunc(x, shift);
    //     std::transform(x.begin(), x.end(), x.begin(),
    //     std::bind(std::multiplies<double>(), std::placeholders::_1,
    //     shift_rate)); return x;
    // }
    // if (rotate_flag) {
    //     std::transform(x.begin(), x.end(), x.begin(),
    //     std::bind(std::multiplies<double>(), std::placeholders::_1,
    //     shift_rate)); return rotatefunc(x, rotate);
    // }
    // std::transform(x.begin(), x.end(), x.begin(),
    // std::bind(std::multiplies<double>(), std::placeholders::_1,
    // shift_rate)); return x;
}

} // namespace

double ackley_func(std::vector<double>& x, const std::vector<double>& shift,
                   const std::vector<std::vector<double>>& rotate,
                   bool shift_flag, bool rotate_flag)
{
    double sum1 = 0.0;
    double sum2 = 0.0;

    const auto y =
        shift_rotate_transform(x, shift, rotate, 1.0, shift_flag, rotate_flag);
    for (auto i : y) {
        sum1 += i * i;
        sum2 += std::cos(2.0 * PI * i);
    }
    sum1 = -0.2 * std::sqrt(sum1 / y.size());
    sum2 /= y.size();
    return E - 20.0 * std::exp(sum1) - std::exp(sum2) + 20.0;
}

double bent_cigar_func(std::vector<double>& x, const std::vector<double>& shift,
                       const std::vector<std::vector<double>>& rotate,
                       bool shift_flag, bool rotate_flag)
{
    const auto y =
        shift_rotate_transform(x, shift, rotate, 1.0, shift_flag, rotate_flag);

    // assuming y.size() >= 1
    return std::accumulate(std::next(y.begin()), y.end(), y[0] * y[0],
                           [](auto f, auto elem) {
                               return f + elem * elem * 1000000.0;
                           }); // 1000000.0 = std::pow(10.0, 6.0)
    // TODO: std::move(f) vs f, which is better?
}

double discus_func(std::vector<double>& x, const std::vector<double>& shift,
                   const std::vector<std::vector<double>>& rotate,
                   bool shift_flag, bool rotate_flag)
{
    const auto y =
        shift_rotate_transform(x, shift, rotate, 1.0, shift_flag, rotate_flag);
    auto f = y[0] * y[0] * 1000000.0;
    // 1000000.0 = std::pow(10.0, 6.0)
    // assuming y.size() >= 1
    return std::accumulate(
        std::next(y.begin()), y.end(), y[0] * y[0] * 1000000.0,
        [](auto f, auto elem) { return std::move(f) + elem * elem; });
}

double ellips_func(std::vector<double>& x, const std::vector<double>& shift,
                   const std::vector<std::vector<double>>& rotate,
                   bool shift_flag, bool rotate_flag)
{
    const auto y =
        shift_rotate_transform(x, shift, rotate, 1.0, shift_flag, rotate_flag);
    auto i = 0;
    const auto n = y.size(); // TODO:  y.size() vs n, which is better?
    return std::accumulate(
        std::next(y.begin()), y.end(), 0.0, [&i, n](auto f, auto elem) {
            return std::move(f) + std::pow(10.0, 6.0 * (i++) / n) * elem * elem;
        });
}

double escaffer6_func(std::vector<double>& x, const std::vector<double>& shift,
                      const std::vector<std::vector<double>>& rotate,
                      bool shift_flag, bool rotate_flag)
{
    const auto y =
        shift_rotate_transform(x, shift, rotate, 1.0, shift_flag, rotate_flag);

    auto f = 0.0;
    for (auto i = 0; i < y.size() - 1; ++i) {
        const auto temp1 =
            std::sin(std::sqrt(y[i] * y[i] + y[i + 1] * y[i + 1]));
        const auto temp2 = 1.0 + 0.001 * (y[i] * y[i] + y[i + 1] * y[i + 1]);
        f += 0.5 + (temp1 * temp1 - 0.5) / (temp2 * temp2);
    }
    const auto temp1 =
        std::sin(std::sqrt(y[y.size() - 1] * y[y.size() - 1] + y[0] * y[0]));
    const auto temp2 =
        1.0 + 0.001 * (y[y.size() - 1] * y[y.size() - 1] + y[0] * y[0]);
    f += 0.5 + (temp1 * temp1 - 0.5) / (temp2 * temp2);
    return f;
}
double griewank_func(std::vector<double>& x, const std::vector<double>& shift,
                     const std::vector<std::vector<double>>& rotate,
                     bool shift_flag, bool rotate_flag)
{
    const auto y = shift_rotate_transform(x, shift, rotate, 600.0 / 100.0,
                                          shift_flag, rotate_flag);
    auto s = 0.0;
    auto p = 1.0;
    for (auto i = 0; i < y.size(); ++i) {
        s += y[i] * y[i];
        p *= std::cos(y[i] / std::sqrt(1.0 + i));
    }
    return 1.0 + s / 4000.0 - p;
}

double grie_rosen_func(std::vector<double>& x, const std::vector<double>& shift,
                       const std::vector<std::vector<double>>& rotate,
                       bool shift_flag, bool rotate_flag)
{
    auto y = shift_rotate_transform(x, shift, rotate, 5.0 / 100.0, shift_flag,
                                    rotate_flag);
    auto f = 0.0;
    y[0] += 1.0; // TODO: we could use const y and resolve 1.0 difference in
                 // calculus
    for (auto i = 0; i < y.size() - 1; ++i) {
        y[i + 1] += 1.0;
        const auto temp1 = y[i] * y[i] - y[i + 1];
        const auto temp2 = y[i] - 1.0;
        const auto temp = 100.0 * temp1 * temp1 + temp2 * temp2;
        f += (temp * temp) / 4000.0 - std::cos(temp) + 1.0;
    }
    const auto temp1 = y[y.size() - 1] * y[y.size() - 1] - y[0];
    const auto temp2 = y[y.size() - 1] - 1.0;
    const auto temp = 100.0 * temp1 * temp1 + temp2 * temp2;
    f += (temp * temp) / 4000.0 - std::cos(temp) + 1.0;
    return f;
}

double happycat_func(std::vector<double>& x, const std::vector<double>& shift,
                     const std::vector<std::vector<double>>& rotate,
                     bool shift_flag, bool rotate_flag)
{
    const auto y = shift_rotate_transform(x, shift, rotate, 5.0 / 100.0,
                                          shift_flag, rotate_flag);
    auto sum_y = 0.0;
    auto r2 = 0.0;
    std::for_each(y.begin(), y.end(), [&](auto elem) {
        const auto aux = elem - 1.0;
        sum_y += aux;
        r2 += aux * aux;
        // TODO: TEST aux vs elem - 1.0
    });
    return std::pow(std::abs(r2 - y.size()), 2 * 1.0 / 8.0) +
           (0.5 * r2 + sum_y) / y.size() + 0.5;
}

double hgbat_func(std::vector<double>& x, const std::vector<double>& shift,
                  const std::vector<std::vector<double>>& rotate,
                  bool shift_flag, bool rotate_flag)
{
    const auto y = shift_rotate_transform(x, shift, rotate, 5.0 / 100.0,
                                          shift_flag, rotate_flag);
    auto sum_y = 0.0;
    auto r2 = 0.0;
    std::for_each(y.begin(), y.end(), [&](auto elem) {
        const auto aux = elem - 1.0;
        sum_y += aux;
        r2 += aux * aux;
    });

    return std::pow(std::fabs(r2 * r2 - sum_y * sum_y), 2.0 * 1.0 / 4.0) +
           (0.5 * r2 + sum_y) / y.size() + 0.5;
}

double rosenbrock_func(std::vector<double>& x, const std::vector<double>& shift,
                       const std::vector<std::vector<double>>& rotate,
                       bool shift_flag, bool rotate_flag)
{
    const auto y = shift_rotate_transform(x, shift, rotate, 2.048 / 100.0,
                                          shift_flag, rotate_flag);
    auto f = 0.0;
    for (auto i = 0; i < y.size() - 1; ++i) {
        // TODO: aux vs not aux
        const auto aux = y[i] + 1.0;
        const auto temp1 = aux * aux - y[i + 1] - 1.0;
        f += 100.0 * temp1 * temp1 + y[i] * y[i];
    }
    // TODO: Make a test method to validate that all these methods are
    // equivalent to the python ones
    return f;
}

double rastrigin_func(std::vector<double>& x, const std::vector<double>& shift,
                      const std::vector<std::vector<double>>& rotate,
                      bool shift_flag, bool rotate_flag)
{
    const auto y = shift_rotate_transform(x, shift, rotate, 5.12 / 100.0,
                                          shift_flag, rotate_flag);
    return std::accumulate(y.begin(), y.end(), 0.0, [=](auto f, auto elem) {
        return f + elem * elem - 10.0 * std::cos(2.0 * PI * elem) + 10;
    });
}

double
schwefel_func(std::vector<double>& x, const std::vector<double>& shift,
              const std::vector<std::vector<double>>& rotate, bool shift_flag,
              bool rotate_flag)
{
    const auto y = shift_rotate_transform(x, shift, rotate, 1000.0 / 100.0,
                                          shift_flag, rotate_flag);
    return std::accumulate(y.begin(), y.end(), 0.0, [&](auto f, auto elem) {
        elem += 4.209687462275036e+002;
        if (elem > 500) {
            f -= (500.0 - std::fmod(elem, 500.0)) * 
            std::sin(std::sqrt(500.0 - std::fmod(elem, 500.0)));
            const auto temp = (elem - 500.0) / 100.0;
            f += temp * temp / y.size();
        } else if (elem < -500) {
            f -= (-500.0 + std::fmod(std::fabs(elem), 500.0)) * 
            std::sin(std::sqrt(500.0 - std::fmod(std::fabs(elem), 500.0)));
            const auto temp = (elem + 500.0) / 100.0;
            f += temp * temp / y.size();
        } else {
            f -= elem * std::sin(std::sqrt(std::fabs(elem)));
        }
        return f + 4.189828872724338e+002 * y.size();
    });
}

int sanity_check()
{
    std::cout << "da";
    return 1;
}

} // namespace ga::functions
