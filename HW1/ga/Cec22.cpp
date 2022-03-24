#include "Cec22.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <numeric>
#include <tuple>
#include <unordered_map>

namespace ga::functions {
namespace {

constexpr double PI = 3.1415926535897932384626433832795029;
constexpr double E = 2.7182818284590452353602874713526625;

void shiftfunc(std::vector<double>& x, const std::vector<double>& shift)
{
    // assuming x.size() == shift.size()
    std::transform(x.begin(), x.end(), shift.begin(), x.begin(),
                   std::minus<double>());
    // TODO: Check what nethods can be vectorized
}

void rotatefunc(std::vector<double>& x, std::vector<double>& aux,
                const std::vector<std::vector<double>>& rotate)
{
    // assuming x.size() == rotate.size() == rotate[0].size() == aux.size()
    std::fill(aux.begin(), aux.end(), 0.0);
    for (std::size_t i = 0; i < x.size(); ++i) {
        for (std::size_t j = 0; j < x.size(); ++j) {
            aux[i] += x[j] * rotate[i][j];
        }
    }
    x.swap(aux);
}

void shift_rotate_transform(std::vector<double>& x, std::vector<double>& aux,
                            const std::vector<double>& shift,
                            const std::vector<std::vector<double>>& rotate,
                            double shift_rate, bool shift_flag,
                            bool rotate_flag)
{
    // template bool and if constexpr
    if (shift_flag) [[likely]] {
        shiftfunc(x, shift);
    }
    std::transform(x.begin(), x.end(), x.begin(),
                   std::bind(std::multiplies<double>(), std::placeholders::_1,
                             shift_rate));
    if (rotate_flag) [[likely]] {
        rotatefunc(x, aux, rotate);
    }
}

void apply_permutation(std::vector<double>& nums,
                       const std::vector<std::size_t>& indices)
{
    std::vector<double> aux;
    aux.reserve(nums.size());
    std::transform(indices.begin(), indices.end(), std::back_inserter(aux),
                   [&](auto i) { return std::move(nums[i]); });
    nums = std::move(aux);
}

template <std::size_t Size>
double composition_function_calculator(const std::vector<double>& x,
                                       const std::vector<double>& shift,
                                       const std::array<int, Size>& delta,
                                       const std::array<double, Size>& fit)
{
    auto w_max = 0.0;
    auto w_sum = 0.0;
    std::array<double, Size> w{};
    for (std::size_t i = 0; i < Size; ++i) {
        for (std::size_t j = 0; j < x.size(); ++j) {
            const auto temp = x[j] - shift[i * x.size() + j];
            w[i] += temp * temp;
        }

        // else will happen only when x is shift
        if (w[i] != 0.0) [[likely]] {
            w[i] = std::sqrt(1.0 / w[i]) *
                   std::exp(-w[i] / 2.0 / x.size() / delta[i] / delta[i]);
        } else [[unlikely]] {
            w[i] = 1.0e99; // INF
        }

        if (w[i] > w_max) {
            w_max = w[i];
        }
        w_sum += w[i];
    }

    // This does not happen if there's any w[i] >= 0
    if (w_max == 0.0) [[unlikely]] {
        return std::accumulate(
            fit.begin(), fit.end(), 0.0,
            [](auto f, auto elem) { return f + elem / Size; });
    }

    // we could calculate w_sum in the previous iteration
    return std::inner_product(
        w.begin(), w.end(), fit.begin(), 0.0, std::plus<double>(),
        [=](auto w, auto fit) { return w / w_sum * fit; });
}

} // namespace

// TODO: visit bellow
double ackley_func(std::vector<double>& x, std::vector<double>& aux,
                   const std::vector<double>& shift,
                   const std::vector<std::vector<double>>& rotate,
                   bool shift_flag, bool rotate_flag)
{
    shift_rotate_transform(x, aux, shift, rotate, 1.0, shift_flag, rotate_flag);
    auto sum1 = 0.0;
    auto sum2 = 0.0;
    for (auto i : x) {
        sum1 += i * i;
        sum2 += std::cos(2.0 * PI * i);
    }
    sum1 = -0.2 * std::sqrt(sum1 / x.size());
    sum2 /= x.size();
    return E - 20.0 * std::exp(sum1) - std::exp(sum2) + 20.0;
}

double bent_cigar_func(std::vector<double>& x, std::vector<double>& aux,
                       const std::vector<double>& shift,
                       const std::vector<std::vector<double>>& rotate,
                       bool shift_flag, bool rotate_flag)
{
    shift_rotate_transform(x, aux, shift, rotate, 1.0, shift_flag, rotate_flag);

    // assuming x.size() >= 1
    return std::accumulate(std::next(x.begin()), x.end(), x[0] * x[0],
                           [](auto f, auto elem) {
                               return f + elem * elem * 1000000.0;
                           }); // 1000000.0 = std::pow(10.0, 6.0)
    // TODO: accumulate vs reduce
}

double discus_func(std::vector<double>& x, std::vector<double>& aux,
                   const std::vector<double>& shift,
                   const std::vector<std::vector<double>>& rotate,
                   bool shift_flag, bool rotate_flag)
{
    shift_rotate_transform(x, aux, shift, rotate, 1.0, shift_flag, rotate_flag);
    // 1000000.0 = std::pow(10.0, 6.0)
    // assuming x.size() >= 1
    return std::accumulate(std::next(x.begin()), x.end(),
                           x[0] * x[0] * 1000000.0,
                           [](auto f, auto elem) { return f + elem * elem; });
}

double ellips_func(std::vector<double>& x, std::vector<double>& aux,
                   const std::vector<double>& shift,
                   const std::vector<std::vector<double>>& rotate,
                   bool shift_flag, bool rotate_flag)
{
    shift_rotate_transform(x, aux, shift, rotate, 1.0, shift_flag, rotate_flag);
    auto i = 0;
    // here we need to accumulate because reduce does not maintain order
    return std::accumulate(
        std::next(x.begin()), x.end(), 0.0,
        [&i, n = x.size()](auto f, auto elem) {
            return std::move(f) + std::pow(10.0, 6.0 * (i++) / n) * elem * elem;
        });
}

double escaffer6_func(std::vector<double>& x, std::vector<double>& aux,
                      const std::vector<double>& shift,
                      const std::vector<std::vector<double>>& rotate,
                      bool shift_flag, bool rotate_flag)
{
    shift_rotate_transform(x, aux, shift, rotate, 1.0, shift_flag, rotate_flag);

    auto f = 0.0;
    // assuming x.size() is not 0
    for (std::size_t i = 0; i < x.size() - 1; ++i) {
        const auto temp1 =
            std::sin(std::sqrt(x[i] * x[i] + x[i + 1] * x[i + 1]));
        const auto temp2 = 1.0 + 0.001 * (x[i] * x[i] + x[i + 1] * x[i + 1]);
        f += 0.5 + (temp1 * temp1 - 0.5) / (temp2 * temp2);
    }
    const auto last = *std::prev(x.end()); // vs x[x.size() - 1]
    const auto temp1 = std::sin(std::sqrt(last * last + x[0] * x[0]));
    const auto temp2 = 1.0 + 0.001 * (last * last + x[0] * x[0]);
    f += 0.5 + (temp1 * temp1 - 0.5) / (temp2 * temp2);
    return f;
}
double griewank_func(std::vector<double>& x, std::vector<double>& aux,
                     const std::vector<double>& shift,
                     const std::vector<std::vector<double>>& rotate,
                     bool shift_flag, bool rotate_flag)
{
    shift_rotate_transform(x, aux, shift, rotate, 600.0 / 100.0, shift_flag,
                           rotate_flag);
    auto s = 0.0;
    auto p = 1.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        s += x[i] * x[i];
        p *= std::cos(x[i] / std::sqrt(1.0 + i));
    }
    return 1.0 + s / 4000.0 - p;
}

double grie_rosen_func(std::vector<double>& x, std::vector<double>& aux,
                       const std::vector<double>& shift,
                       const std::vector<std::vector<double>>& rotate,
                       bool shift_flag, bool rotate_flag)
{
    shift_rotate_transform(x, aux, shift, rotate, 5.0 / 100.0, shift_flag,
                           rotate_flag);
    auto f = 0.0;
    x[0] += 1.0;
    // assuming x.size() is not 0
    for (std::size_t i = 0; i < x.size() - 1; ++i) {
        x[i + 1] += 1.0;
        const auto temp1 = x[i] * x[i] - x[i + 1];
        const auto temp2 = x[i] - 1.0;
        const auto temp = 100.0 * temp1 * temp1 + temp2 * temp2;
        f += (temp * temp) / 4000.0 - std::cos(temp) + 1.0;
    }
    const auto last = *std::prev(x.end());
    const auto temp1 = last * last - x[0];
    const auto temp2 = last - 1.0;
    const auto temp = 100.0 * temp1 * temp1 + temp2 * temp2;
    f += (temp * temp) / 4000.0 - std::cos(temp) + 1.0;
    return f;
}

double happycat_func(std::vector<double>& x, std::vector<double>& aux,
                     const std::vector<double>& shift,
                     const std::vector<std::vector<double>>& rotate,
                     bool shift_flag, bool rotate_flag)
{
    shift_rotate_transform(x, aux, shift, rotate, 5.0 / 100.0, shift_flag,
                           rotate_flag);
    auto sum_y = 0.0;
    auto r2 = 0.0;
    std::for_each(x.begin(), x.end(), [&](auto elem) {
        const auto aux = elem - 1.0;
        sum_y += aux;
        r2 += aux * aux;
    });
    return std::pow(std::abs(r2 - x.size()), 2 * 1.0 / 8.0) +
           (0.5 * r2 + sum_y) / x.size() + 0.5;
}

double hgbat_func(std::vector<double>& x, std::vector<double>& aux,
                  const std::vector<double>& shift,
                  const std::vector<std::vector<double>>& rotate,
                  bool shift_flag, bool rotate_flag)
{
    shift_rotate_transform(x, aux, shift, rotate, 5.0 / 100.0, shift_flag,
                           rotate_flag);
    auto sum_y = 0.0;
    auto r2 = 0.0;
    std::for_each(x.begin(), x.end(), [&](auto elem) {
        const auto aux = elem - 1.0;
        sum_y += aux;
        r2 += aux * aux;
    });

    return std::pow(std::fabs(r2 * r2 - sum_y * sum_y), 2.0 * 1.0 / 4.0) +
           (0.5 * r2 + sum_y) / x.size() + 0.5;
}

double rosenbrock_func(std::vector<double>& x, std::vector<double>& aux,
                       const std::vector<double>& shift,
                       const std::vector<std::vector<double>>& rotate,
                       bool shift_flag, bool rotate_flag)
{
    shift_rotate_transform(x, aux, shift, rotate, 2.048 / 100.0, shift_flag,
                           rotate_flag);
    auto f = 0.0;
    // assuming x.size() is not 0
    for (std::size_t i = 0; i < x.size() - 1; ++i) {
        const auto aux = x[i] + 1.0;
        const auto temp1 = aux * aux - x[i + 1] - 1.0;
        f += 100.0 * temp1 * temp1 + x[i] * x[i];
    }
    // TODO: Make a test method to validate that all these methods are
    // equivalent to the python ones
    return f;
}

double rastrigin_func(std::vector<double>& x, std::vector<double>& aux,
                      const std::vector<double>& shift,
                      const std::vector<std::vector<double>>& rotate,
                      bool shift_flag, bool rotate_flag)
{
    shift_rotate_transform(x, aux, shift, rotate, 5.12 / 100.0, shift_flag,
                           rotate_flag);
    return std::accumulate(
        x.begin(), x.end(), 10.0 * x.size(), [=](auto f, auto elem) {
            return f + elem * elem - 10.0 * std::cos(2.0 * PI * elem);
        });
}

double schwefel_func(std::vector<double>& x, std::vector<double>& aux,
                     const std::vector<double>& shift,
                     const std::vector<std::vector<double>>& rotate,
                     bool shift_flag, bool rotate_flag)
{
    shift_rotate_transform(x, aux, shift, rotate, 1000.0 / 100.0, shift_flag,
                           rotate_flag);
    return std::accumulate(
        x.begin(), x.end(), 0.0, [n = x.size()](auto f, auto elem) {
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

double schaffer_F7_func(std::vector<double>& x, std::vector<double>& aux,
                        const std::vector<double>& shift,
                        const std::vector<std::vector<double>>& rotate,
                        bool shift_flag, bool rotate_flag)
{
    shift_rotate_transform(x, aux, shift, rotate, 1.0, shift_flag, rotate_flag);
    // shaffer_F7_func is tottaly different from original, which is broken
    // source:
    // https://www.sciencedirect.com/topics/computer-science/benchmark-function
    // TODO: check if results are ok
    auto f = 0.0;
    // assuming x.size() is not 0
    for (std::size_t i = 0; i < x.size() - 1; ++i) {
        const auto si = std::sqrt(x[i] * x[i] + x[i + 1] * x[i + 1]);
        const auto temp = std::sin(50.0 * std::pow(si, 0.2));
        f += si + si * temp * temp;
    }
    return f * f / x.size() / x.size();
}

double step_rastrigin_func(std::vector<double>& x, std::vector<double>& aux,
                           const std::vector<double>& shift,
                           const std::vector<std::vector<double>>& rotate,
                           bool shift_flag, bool rotate_flag)
{
    // ??? is exactly step rastrigin if shift_flag and rotate_flag
    // TODO: maybe remove shift and rotate flag and always shift and rotate
    // see h01
    return rastrigin_func(x, aux, shift, rotate, shift_flag, rotate_flag);
}

double levy_func(std::vector<double>& x, std::vector<double>& aux,
                 const std::vector<double>& shift,
                 const std::vector<std::vector<double>>& rotate,
                 bool shift_flag, bool rotate_flag)
{
    shift_rotate_transform(x, aux, shift, rotate, 1.0, shift_flag, rotate_flag);

    const auto w = [](auto elem) { return 1.0 + (elem - 1.0) / 4.0; };

    const auto term1 = std::sin(PI * w(x[0]));
    const auto term2 = std::accumulate(
        x.begin(), std::prev(x.end()), 0.0, [&w](auto f, auto elem) {
            const auto wi = w(elem);
            const auto temp = std::sin(PI * wi + 1.0);
            return f + (wi - 1.0) * (wi - 1.0) * (1.0 + 10.0 * temp * temp);
        });
    const auto last = w(*std::prev(x.end()));
    const auto temp = std::sin(2.0 * PI * last);
    const auto term3 = (last - 1.0) * (last - 1.0) * (1.0 + temp * temp);

    return term1 * term1 + term2 + term3;
}

double zakharov_func(std::vector<double>& x, std::vector<double>& aux,
                     const std::vector<double>& shift,
                     const std::vector<std::vector<double>>& rotate,
                     bool shift_flag, bool rotate_flag)
{
    shift_rotate_transform(x, aux, shift, rotate, 1.0, shift_flag, rotate_flag);
    auto sum1 = 0.0;
    auto sum2 = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        sum1 += x[i] * x[i];
        sum2 += 0.5 * i * x[i];
    }
    return sum1 + sum2 * sum2 + sum2 * sum2 * sum2 * sum2;
}

double katsuura_func(std::vector<double>& x, std::vector<double>& aux,
                     const std::vector<double>& shift,
                     const std::vector<std::vector<double>>& rotate,
                     bool shift_flag, bool rotate_flag)
{
    shift_rotate_transform(x, aux, shift, rotate, 5.0 / 100.0, shift_flag,
                           rotate_flag);

    auto f = 1.0;
    const auto size = x.size();
    const auto temp3 = std::pow(static_cast<double>(size), 1.2);
    for (std::size_t i = 0; i < size; ++i) {
        auto temp = 0.0;
        for (auto j = 1; j < 33; ++j) {
            const auto temp1 = std::pow(2.0, j);
            const auto temp2 = temp1 * x[i];
            temp += std::fabs(temp2 - std::floor(temp2 + 0.5)) / temp1;
        }
        f *= std::pow(1.0 + (i + 1) * temp, 10.0 / temp3);
    }
    const auto temp1 = 10.0 / size / size;
    return f * temp1 - temp1;
}

double
hf01(std::vector<double>& x, std::vector<double>& aux,
     const std::vector<double>& shift,
     const std::vector<std::vector<double>>& rotate,
     const std::vector<std::size_t>& indices, bool shift_flag, bool rotate_flag)
{
    // [0.4, 0.4, 0.2]
    shift_rotate_transform(x, aux, shift, rotate, 1.0, shift_flag, rotate_flag);
    apply_permutation(x, indices);

    const auto range = std::ceil(0.4 * x.size());
    const auto margin_1 = std::next(x.begin(), range); // 0.4
    const auto margin_2 = std::next(margin_1, range);  // another 0.4
    const auto aux_margin_1 = std::next(aux.begin(), range);
    const auto aux_margin_2 = std::next(margin_1, range);

    std::vector<double> z1{x.begin(), margin_1};
    std::vector<double> z2{margin_1, margin_2};
    std::vector<double> z3{margin_2, x.end()};
    std::vector<double> aux1(z1.size(), 0.0);
    std::vector<double> aux2(z2.size(), 0.0);
    std::vector<double> aux3(z3.size(), 0.0);

    return bent_cigar_func(z1, aux1, shift, rotate, false, false) +
           hgbat_func(z2, aux2, shift, rotate, false, false) +
           rastrigin_func(z3, aux3, shift, rotate, false, false);
}

double
hf02(std::vector<double>& x, std::vector<double>& aux,
     const std::vector<double>& shift,
     const std::vector<std::vector<double>>& rotate,
     const std::vector<std::size_t>& indices, bool shift_flag, bool rotate_flag)
{
    throw "Not Implemented";
}

double
hf03(std::vector<double>& x, std::vector<double>& aux,
     const std::vector<double>& shift,
     const std::vector<std::vector<double>>& rotate,
     const std::vector<std::size_t>& indices, bool shift_flag, bool rotate_flag)
{
    throw "Not Implemented";
}

double cf01(std::vector<double>& x, std::vector<double>& aux,
            const std::vector<double>& shift,
            const std::vector<std::vector<double>>& rotate, bool rotate_flag)
{
    // asumming that shift.size and rotate.size is (at least) 5 * x.size
    const auto shift_margin_1 = std::next(shift.begin(), x.size());
    const auto shift_margin_2 = std::next(shift_margin_1, x.size());
    const auto shift_margin_3 = std::next(shift_margin_2, x.size());
    const auto shift_margin_4 = std::next(shift_margin_3, x.size());

    const auto rotate_margin_1 = std::next(rotate.begin(), x.size());
    const auto rotate_margin_2 = std::next(rotate_margin_1, x.size());
    const auto rotate_margin_3 = std::next(rotate_margin_2, x.size());
    const auto rotate_margin_4 = std::next(rotate_margin_3, x.size());

    constexpr auto N = 5;
    // fit is function result * lambda + bias
    // lambda is 1, 1e-6, 1e-6, 1e-6, 1e-6
    // bias is 0, 200, 300, 100, 400
    const std::array<double, N> fit{
        rosenbrock_func(x, aux, {shift.begin(), shift_margin_1},
                        {rotate.begin(), rotate_margin_1}, true, rotate_flag),
        ellips_func(x, aux, {shift_margin_1, shift_margin_2},
                    {rotate_margin_1, rotate_margin_2}, true, rotate_flag) *
                1e-6 +
            200,
        bent_cigar_func(x, aux, {shift_margin_2, shift_margin_3},
                        {rotate_margin_2, rotate_margin_3}, true, rotate_flag) *
                1e-6 +
            300,
        discus_func(x, aux, {shift_margin_3, shift_margin_4},
                    {rotate_margin_3, rotate_margin_4}, true, rotate_flag) *
                1e-6 +
            100,
        ellips_func(x, aux, {shift_margin_4, shift.end()},
                    {rotate_margin_4, rotate.end()}, true, false) *
                1e-6 +
            400, // ?? why false
    };
    const std::array<int, N> delta{10, 20, 30, 40, 50};
    return composition_function_calculator<N>(x, shift, delta, fit);
}

double cf02(std::vector<double>& x, std::vector<double>& aux,
            const std::vector<double>& shift,
            const std::vector<std::vector<double>>& rotate, bool rotate_flag)
{
    // asumming that shift.size and rotate.size is (at least) 3 * x.size
    const auto shift_margin_1 = std::next(shift.begin(), x.size());
    const auto shift_margin_2 = std::next(shift_margin_1, x.size());

    const auto rotate_margin_1 = std::next(rotate.begin(), x.size());
    const auto rotate_margin_2 = std::next(rotate_margin_1, x.size());

    constexpr auto N = 3;
    // fit is function result * lambda + bias
    // lambda is 1, 1, 1
    // bias is 0, 200, 100
    // TODO: use iterators, do not create new vectors
    const std::array<double, N> fit{
        schwefel_func(x, aux, {shift.begin(), shift_margin_1},
                      {rotate.begin(), rotate_margin_1}, true,
                      false), // ?? why false
        rastrigin_func(x, aux, {shift_margin_1, shift_margin_2},
                       {rotate_margin_1, rotate_margin_2}, true, rotate_flag) +
            200,
        hgbat_func(x, aux, {shift_margin_2, shift.end()},
                   {rotate_margin_2, rotate.end()}, true, rotate_flag) +
            100,
    };
    const std::array<int, N> delta{20, 10, 10};
    return composition_function_calculator<N>(x, shift, delta, fit);
}

double cf03(std::vector<double>& x, std::vector<double>& aux,
            const std::vector<double>& shift,
            const std::vector<std::vector<double>>& rotate, bool rotate_flag)
{
    // asumming that shift.size and rotate.size is (at least) 3 * x.size
    throw "Not Implemented";
}

double cf04(std::vector<double>& x, std::vector<double>& aux,
            const std::vector<double>& shift,
            const std::vector<std::vector<double>>& rotate, bool rotate_flag)
{
    // asumming that shift.size and rotate.size is (at least) 3 * x.size
    throw "Not Implemented";
}

int sanity_check()
{
    std::cout << "da";
    std::vector<double> a = {1.1, 1.2};
    std::vector<double> aux = {1.1, 1.2};
    std::vector<double> b = {0.1, 0.2};
    std::vector<std::vector<double>> c = {{1.1, 1.2}, {1.1, 1.2}};

    std::cout << 'a' << " = "
              << ga::functions::rastrigin_func(a, aux, b, c, true, true)
              << '\n';

    std::vector<double> x = {1.1, 1.2};
    std::vector<double> shift = {0.1, 0.2, 0.11, 0.21, 0.12, 0.22};
    std::vector<std::vector<double>> rotate = {
        {1.1, 1.2}, {1.1, 1.2}, {1.1, 1.2}, {1.1, 1.2}, {1.1, 1.2}, {1.1, 1.2}};
    std::cout << "b = " << ga::functions::cf02(x, aux, shift, rotate, true)
              << '\n';
    return 1;
}

} // namespace ga::functions
