module;
#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <numeric>
#include <tuple>
#include <unordered_map>
module GA.Cec22Specific;
// clang-format off
// import <iostream>;
// import <algorithm>;
// import <functional>;
// import <cmath>;
// import <numeric>;
// import <unordered_map>;
// import <tuple>;
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

void shuffle(std::vector<double>& nums, const std::vector<std::size_t>& pos)
{
    // assuming that pos is a permutation of [0, .., nums.size() - 1]
    double aux[pos.size()];
    for (auto i = 0; i < pos.size(); ++i) {
        aux[pos[i]] = nums[i];
    }
    for (auto i = 0; i < pos.size(); ++i) {
        nums[i] = aux[i];
    }
}

template <std::size_t Size>
double composition_function_calculator(const std::vector<double>& x,
                                       const std::vector<double>& shift,
                                       const std::array<int, Size>& delta,
                                       const std::array<int, Size>& bias,
                                       const std::array<double, Size>& fit)
{
    throw "Not Implemented";
}

} // namespace

double ackley_func(std::vector<double>& x, const std::vector<double>& shift,
                   const std::vector<std::vector<double>>& rotate,
                   bool shift_flag, bool rotate_flag)
{
    const auto z =
        shift_rotate_transform(x, shift, rotate, 1.0, shift_flag, rotate_flag);
    auto sum1 = 0.0;
    auto sum2 = 0.0;
    for (auto i : z) {
        sum1 += i * i;
        sum2 += std::cos(2.0 * PI * i);
    }
    sum1 = -0.2 * std::sqrt(sum1 / z.size());
    sum2 /= z.size();
    return E - 20.0 * std::exp(sum1) - std::exp(sum2) + 20.0;
}

double bent_cigar_func(std::vector<double>& x, const std::vector<double>& shift,
                       const std::vector<std::vector<double>>& rotate,
                       bool shift_flag, bool rotate_flag)
{
    const auto z =
        shift_rotate_transform(x, shift, rotate, 1.0, shift_flag, rotate_flag);

    // assuming z.size() >= 1
    return std::accumulate(std::next(z.begin()), z.end(), z[0] * z[0],
                           [](auto f, auto elem) {
                               return f + elem * elem * 1000000.0;
                           }); // 1000000.0 = std::pow(10.0, 6.0)
    // TODO: std::move(f) vs f, which is better?
}

double discus_func(std::vector<double>& x, const std::vector<double>& shift,
                   const std::vector<std::vector<double>>& rotate,
                   bool shift_flag, bool rotate_flag)
{
    const auto z =
        shift_rotate_transform(x, shift, rotate, 1.0, shift_flag, rotate_flag);
    auto f = z[0] * z[0] * 1000000.0;
    // 1000000.0 = std::pow(10.0, 6.0)
    // assuming z.size() >= 1
    return std::accumulate(
        std::next(z.begin()), z.end(), z[0] * z[0] * 1000000.0,
        [](auto f, auto elem) { return f + elem * elem; });
}

double ellips_func(std::vector<double>& x, const std::vector<double>& shift,
                   const std::vector<std::vector<double>>& rotate,
                   bool shift_flag, bool rotate_flag)
{
    const auto z =
        shift_rotate_transform(x, shift, rotate, 1.0, shift_flag, rotate_flag);
    auto i = 0;
    return std::accumulate(
        std::next(z.begin()), z.end(), 0.0, [&i, n = z.size()](auto f, auto elem) {
            return std::move(f) + std::pow(10.0, 6.0 * (i++) / n) * elem * elem;
        });
}

double escaffer6_func(std::vector<double>& x, const std::vector<double>& shift,
                      const std::vector<std::vector<double>>& rotate,
                      bool shift_flag, bool rotate_flag)
{
    const auto z =
        shift_rotate_transform(x, shift, rotate, 1.0, shift_flag, rotate_flag);

    auto f = 0.0;
    for (auto i = 0; i < z.size() - 1; ++i) {
        const auto temp1 =
            std::sin(std::sqrt(z[i] * z[i] + z[i + 1] * z[i + 1]));
        const auto temp2 = 1.0 + 0.001 * (z[i] * z[i] + z[i + 1] * z[i + 1]);
        f += 0.5 + (temp1 * temp1 - 0.5) / (temp2 * temp2);
    }
    const auto last = *std::prev(z.end()); // vs z[z.size() - 1]
    const auto temp1 =
        std::sin(std::sqrt(last * last + z[0] * z[0]));
    const auto temp2 =
        1.0 + 0.001 * (last * last + z[0] * z[0]);
    f += 0.5 + (temp1 * temp1 - 0.5) / (temp2 * temp2);
    return f;
}
double griewank_func(std::vector<double>& x, const std::vector<double>& shift,
                     const std::vector<std::vector<double>>& rotate,
                     bool shift_flag, bool rotate_flag)
{
    const auto z = shift_rotate_transform(x, shift, rotate, 600.0 / 100.0,
                                          shift_flag, rotate_flag);
    auto s = 0.0;
    auto p = 1.0;
    for (auto i = 0; i < z.size(); ++i) {
        s += z[i] * z[i];
        p *= std::cos(z[i] / std::sqrt(1.0 + i));
    }
    return 1.0 + s / 4000.0 - p;
}

double grie_rosen_func(std::vector<double>& x, const std::vector<double>& shift,
                       const std::vector<std::vector<double>>& rotate,
                       bool shift_flag, bool rotate_flag)
{
    auto z = shift_rotate_transform(x, shift, rotate, 5.0 / 100.0, shift_flag,
                                    rotate_flag);
    auto f = 0.0;
    z[0] += 1.0; // TODO: we could use const z and resolve 1.0 difference in
                 // calculus
    for (auto i = 0; i < z.size() - 1; ++i) {
        z[i + 1] += 1.0;
        const auto temp1 = z[i] * z[i] - z[i + 1];
        const auto temp2 = z[i] - 1.0;
        const auto temp = 100.0 * temp1 * temp1 + temp2 * temp2;
        f += (temp * temp) / 4000.0 - std::cos(temp) + 1.0;
    }
    const auto temp1 = z[z.size() - 1] * z[z.size() - 1] - z[0];
    const auto temp2 = z[z.size() - 1] - 1.0;
    const auto temp = 100.0 * temp1 * temp1 + temp2 * temp2;
    f += (temp * temp) / 4000.0 - std::cos(temp) + 1.0;
    return f;
}

double happycat_func(std::vector<double>& x, const std::vector<double>& shift,
                     const std::vector<std::vector<double>>& rotate,
                     bool shift_flag, bool rotate_flag)
{
    const auto z = shift_rotate_transform(x, shift, rotate, 5.0 / 100.0,
                                          shift_flag, rotate_flag);
    auto sum_y = 0.0;
    auto r2 = 0.0;
    std::for_each(z.begin(), z.end(), [&](auto elem) {
        const auto aux = elem - 1.0;
        sum_y += aux;
        r2 += aux * aux;
        // TODO: TEST aux vs elem - 1.0
    });
    return std::pow(std::abs(r2 - z.size()), 2 * 1.0 / 8.0) +
           (0.5 * r2 + sum_y) / z.size() + 0.5;
}

double hgbat_func(std::vector<double>& x, const std::vector<double>& shift,
                  const std::vector<std::vector<double>>& rotate,
                  bool shift_flag, bool rotate_flag)
{
    const auto z = shift_rotate_transform(x, shift, rotate, 5.0 / 100.0,
                                          shift_flag, rotate_flag);
    auto sum_y = 0.0;
    auto r2 = 0.0;
    std::for_each(z.begin(), z.end(), [&](auto elem) {
        const auto aux = elem - 1.0;
        sum_y += aux;
        r2 += aux * aux;
    });

    return std::pow(std::fabs(r2 * r2 - sum_y * sum_y), 2.0 * 1.0 / 4.0) +
           (0.5 * r2 + sum_y) / z.size() + 0.5;
}

double rosenbrock_func(std::vector<double>& x, const std::vector<double>& shift,
                       const std::vector<std::vector<double>>& rotate,
                       bool shift_flag, bool rotate_flag)
{
    const auto z = shift_rotate_transform(x, shift, rotate, 2.048 / 100.0,
                                          shift_flag, rotate_flag);
    auto f = 0.0;
    for (auto i = 0; i < z.size() - 1; ++i) {
        // TODO: aux vs not aux
        const auto aux = z[i] + 1.0;
        const auto temp1 = aux * aux - z[i + 1] - 1.0;
        f += 100.0 * temp1 * temp1 + z[i] * z[i];
    }
    // TODO: Make a test method to validate that all these methods are
    // equivalent to the python ones
    return f;
}

double rastrigin_func(std::vector<double>& x, const std::vector<double>& shift,
                      const std::vector<std::vector<double>>& rotate,
                      bool shift_flag, bool rotate_flag)
{
    const auto z = shift_rotate_transform(x, shift, rotate, 5.12 / 100.0,
                                          shift_flag, rotate_flag);
    return std::accumulate(z.begin(), z.end(), 0.0, [=](auto f, auto elem) {
        return f + elem * elem - 10.0 * std::cos(2.0 * PI * elem) + 10;
    });
}

double schwefel_func(std::vector<double>& x, const std::vector<double>& shift,
                     const std::vector<std::vector<double>>& rotate,
                     bool shift_flag, bool rotate_flag)
{
    const auto z = shift_rotate_transform(x, shift, rotate, 1000.0 / 100.0,
                                          shift_flag, rotate_flag);
    return std::accumulate(
        z.begin(), z.end(), 0.0, [n = z.size()](auto f, auto elem) {
            elem += 4.209687462275036e+002;
            if (elem > 500) {
                f -= (500.0 - std::fmod(elem, 500.0)) *
                     std::sin(std::sqrt(500.0 - std::fmod(elem, 500.0)));
                const auto temp = (elem - 500.0) / 100.0;
                f += temp * temp / n;
            } else if (elem < -500) {
                f -= (-500.0 + std::fmod(std::fabs(elem), 500.0)) *
                     std::sin(
                         std::sqrt(500.0 - std::fmod(std::fabs(elem), 500.0)));
                const auto temp = (elem + 500.0) / 100.0;
                f += temp * temp / n;
            } else {
                f -= elem * std::sin(std::sqrt(std::fabs(elem)));
            }
            return f + 4.189828872724338e+002 * n;
        });
}

double
schaffer_F7_func(std::vector<double>& x, const std::vector<double>& shift,
                 const std::vector<std::vector<double>>& rotate,
                 bool shift_flag, bool rotate_flag)
{
    const auto z =
        shift_rotate_transform(x, shift, rotate, 1.0, shift_flag, rotate_flag);
    // shaffer_F7_func is tottaly different from original, which is broken
    // source:
    // https://www.sciencedirect.com/topics/computer-science/benchmark-function
    // TODO: check if results are ok
    auto f = 0.0;
    for (auto i = 0; i < z.size() - 1; ++i) {
        const auto si = std::sqrt(z[i] * z[i] + z[i + 1] * z[i + 1]);
        const auto temp = std::sin(50.0 * std::pow(si, 0.2));
        f += si + si * temp * temp;
    }
    return f * f / z.size() / z.size();
}

double
step_rastrigin_func(std::vector<double>& x, const std::vector<double>& shift,
                    const std::vector<std::vector<double>>& rotate,
                    bool shift_flag, bool rotate_flag)
{
    // ??? is exactly step rastrigin if shift_flag and rotate_flag
    // TODO: maybe remove shift and rotate flag and always shift and rotate
    // see h01
    return rastrigin_func(x, shift, rotate, shift_flag, rotate_flag);
}

double levy_func(std::vector<double>& x, const std::vector<double>& shift,
                 const std::vector<std::vector<double>>& rotate,
                 bool shift_flag, bool rotate_flag)
{
    const auto z =
        shift_rotate_transform(x, shift, rotate, 1.0, shift_flag, rotate_flag);
    const auto w = [](auto elem) { return 1.0 + (elem - 1.0) / 4.0; };
    const auto term1 = std::sin(PI * w(z[1]));
    const auto term2 = std::accumulate(
        z.begin(), std::prev(z.end()), 0.0, [&w](auto f, auto elem) {
            const auto wi = w(elem);
            const auto temp = std::sin(PI * wi + 1.0);
            return f + (wi - 1.0) * (wi - 1.0) * (1.0 + 10.0 * temp * temp);
        });
    const auto last = w(z[z.size() - 1]);
    const auto temp = std::sin(2.0 * PI * last);
    const auto term3 = (last - 1.0) * (last - 1.0) * (1.0 + temp * temp);
    return term1 * term1 + term2 + term3;
}

double zakharov_func(std::vector<double>& x, const std::vector<double>& shift,
                     const std::vector<std::vector<double>>& rotate,
                     bool shift_flag, bool rotate_flag)
{
    const auto z =
        shift_rotate_transform(x, shift, rotate, 1.0, shift_flag, rotate_flag);
    auto sum1 = 0.0;
    auto sum2 = 0.0;
    for (auto i = 0; i < z.size(); ++i) {
        sum1 += z[i] * z[i];
        sum2 += 0.5 * i * z[i];
    }
    return sum1 + sum2 * sum2 + sum2 * sum2 * sum2 * sum2;
}

double katsuura_func(std::vector<double>& x, const std::vector<double>& shift,
                     const std::vector<std::vector<double>>& rotate,
                     bool shift_flag, bool rotate_flag)
{
    const auto z = shift_rotate_transform(x, shift, rotate, 5.0 / 100.0,
                                          shift_flag, rotate_flag);
    auto f = 1.0;
    const auto temp3 = std::pow(static_cast<double>(z.size()), 1.2);
    for (auto i = 0; i < z.size(); ++i) {
        auto temp = 0.0;
        for (auto j = 1; j < 33; ++j) {
            const auto temp1 = std::pow(2.0, j);
            const auto temp2 = temp1 * z[i];
            temp += std::fabs(temp2 - std::floor(temp2 + 0.5)) / temp1;
        }
        f *= std::pow(1.0 + (i + 1) * temp, 10.0 / temp3);
    }
    const auto temp1 = 10.0 / z.size() / z.size();
    return f * temp1 - temp1;
}

double
hf01(std::vector<double>& x, const std::vector<double>& shift,
     const std::vector<std::vector<double>>& rotate,
     const std::vector<std::size_t>& indices, bool shift_flag, bool rotate_flag)
{
    // [0.4, 0.4, 0.2]
    constexpr auto N = 3;
    auto z =
        shift_rotate_transform(x, shift, rotate, 1.0, shift_flag, rotate_flag);
    shuffle(z, indices);
    const auto range = std::ceil(0.4 * z.size());
    const auto end1 = std::next(z.begin(), range); // 0.4
    const auto begin2 = std::next(end1);
    const auto end2 = std::next(begin2, range); // another 0.4
    const auto begin3 = std::next(end2);

    std::vector<double> z1{z.begin(), end1};
    std::vector<double> z2{begin2, end2};
    std::vector<double> z3{begin3, z.end()};
    return bent_cigar_func(z1, shift, rotate, false, false) +
           hgbat_func(z2, shift, rotate, false, false) +
           rastrigin_func(z3, shift, rotate, false, false);
}

double
hf02(std::vector<double>& x, const std::vector<double>& shift,
     const std::vector<std::vector<double>>& rotate,
     const std::vector<std::size_t>& indices, bool shift_flag, bool rotate_flag)
{
    throw "Not Implemented";
}

double
hf03(std::vector<double>& x, const std::vector<double>& shift,
     const std::vector<std::vector<double>>& rotate,
     const std::vector<std::size_t>& indices, bool shift_flag, bool rotate_flag)
{
    throw "Not Implemented";
}

double cf01(std::vector<double>& x, const std::vector<double>& shift,
            const std::vector<std::vector<double>>& rotate, bool shift_flag,
            bool rotate_flag)
{
    // asumming that shift.size and rotate.size is (at least) 5 * x.size
    const auto shift_end_1 = std::next(shift.begin(), x.size());
    const auto shift_begin_2 = std::next(shift_end_1);
    const auto shift_end_2 = std::next(shift_begin_2, x.size());
    const auto shift_begin_3 = std::next(shift_end_2);
    const auto shift_end_3 = std::next(shift_begin_3, x.size());
    const auto shift_begin_4 = std::next(shift_end_3);
    const auto shift_end_4 = std::next(shift_begin_4, x.size());
    const auto shift_begin_5 = std::next(shift_end_4);
    const auto rotate_end_1 = std::next(rotate.begin(), x.size());
    const auto rotate_begin_2 = std::next(rotate_end_1);
    const auto rotate_end_2 = std::next(rotate_begin_2, x.size());
    const auto rotate_begin_3 = std::next(rotate_end_2);
    const auto rotate_end_3 = std::next(rotate_begin_3, x.size());
    const auto rotate_begin_4 = std::next(rotate_end_3);
    const auto rotate_end_4 = std::next(rotate_begin_4, x.size());
    const auto rotate_begin_5 = std::next(rotate_end_4);
    constexpr auto N = 5;
    const std::array<double, N> fit{
        rosenbrock_func(x, {shift.begin(), shift_end_1},
                        {rotate.begin(), rotate_end_1}, true, rotate_flag),
        ellips_func(x, {shift_begin_2, shift_end_2},
                    {rotate_begin_2, rotate_end_2}, true, rotate_flag) *
            1e-6,
        bent_cigar_func(x, {shift_begin_3, shift_end_3},
                        {rotate_begin_3, rotate_end_3}, true, rotate_flag) *
            1e-6,
        discus_func(x, {shift_begin_4, shift_end_4},
                    {rotate_begin_4, rotate_end_4}, true, rotate_flag) *
            1e-6,
        ellips_func(x, {shift_begin_5, shift.end()},
                    {rotate_begin_5, rotate.end()}, true, false) *
            1e-6, // ?? why false
    };
    const std::array<int, N> delta{10, 20, 30, 40, 50};
    const std::array<int, N> bias{0, 200, 300, 100, 400};
    return composition_function_calculator<N>(x, shift, delta, bias, fit);
}

double cf02(std::vector<double>& x, const std::vector<double>& shift,
            const std::vector<std::vector<double>>& rotate, bool shift_flag,
            bool rotate_flag)
{
    // asumming that shift.size and rotate.size is (at least) 3 * x.size
    const auto shift_end_1 = std::next(shift.begin(), x.size());
    const auto shift_begin_2 = std::next(shift_end_1);
    const auto shift_end_2 = std::next(shift_begin_2, x.size());
    const auto shift_begin_3 = std::next(shift_end_2);
    const auto rotate_end_1 = std::next(rotate.begin(), x.size());
    const auto rotate_begin_2 = std::next(rotate_end_1);
    const auto rotate_end_2 = std::next(rotate_begin_2, x.size());
    const auto rotate_begin_3 = std::next(rotate_end_2);
    constexpr auto N = 3;
    const std::array<double, N> fit{
        schwefel_func(x, {shift.begin(), shift_end_1},
                      {rotate.begin(), rotate_end_1}, true,
                      false), // ?? why false
        rastrigin_func(x, {shift_begin_2, shift_end_2},
                       {rotate_begin_2, rotate_end_2}, true, rotate_flag),
        hgbat_func(x, {shift_begin_3, shift.end()},
                   {rotate_begin_3, rotate.end()}, true, rotate_flag),
    };
    const std::array<int, N> delta{20, 10, 10};
    const std::array<int, N> bias{0, 200, 100};
    return composition_function_calculator<N>(x, shift, delta, bias, fit);
}

double cf03(std::vector<double>& x, const std::vector<double>& shift,
            const std::vector<std::vector<double>>& rotate, bool shift_flag,
            bool rotate_flag)
{
    // asumming that shift.size and rotate.size is (at least) 3 * x.size
    throw "Not Implemented";
}

double cf04(std::vector<double>& x, const std::vector<double>& shift,
            const std::vector<std::vector<double>>& rotate, bool shift_flag,
            bool rotate_flag)
{
    // asumming that shift.size and rotate.size is (at least) 3 * x.size
    throw "Not Implemented";
}

int sanity_check()
{
    std::cout << "da";
    return 1;
}

} // namespace ga::functions
