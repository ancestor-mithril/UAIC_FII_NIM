#include "Cec22.h"

#include <algorithm>
#include <cmath>
#include <execution>
#include <functional>
#include <iostream>
#include <numeric>

namespace cec22 {
namespace {

constexpr double PI = 3.1415926535897932384626433832795029;
constexpr double E = 2.7182818284590452353602874713526625;

using vector_begin = std::vector<double>::const_iterator;
using matrix_begin = std::vector<std::vector<double>>::const_iterator;

struct VectorRange {
    std::vector<double>::iterator begin;
    std::vector<double>::iterator end;
};

// void print_vec(const std::vector<double>& x)
// {
//     for (auto i : x) {
//         std::cout << i << ' ';
//     }
//     std::cout << '\n';
// }

void shiftfunc(const std::vector<double>& x, std::vector<double>& aux,
               const vector_begin shiftBegin)
{
    // assuming x.size() == shift.size()
    std::transform(std::execution::unseq, x.begin(), x.end(), shiftBegin,
                   aux.begin(), std::minus<double>());
}

void rotatefunc(std::vector<double>& vec, const matrix_begin rotateBegin)
{
    // assuming vec.size() == rotate.size() == rotate[0].size() == aux.size()
    const auto n = vec.size();
    auto aux = std::vector<double>(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            aux[i] += vec[j] * rotateBegin[i][j];
            // accesing iterator as if it were a matrix
        }
    }
    std::swap(aux, vec); // swapping result to back to vec
}

// const VectorRange is misleading. Value of iterators can change, the position
// won't
void justShift(const VectorRange& x, double shiftRate)
{
    std::transform(
        x.begin, x.end, x.begin,
        std::bind(std::multiplies<double>(), std::placeholders::_1, shiftRate));
}

void shift(const std::vector<double>& x, std::vector<double>& aux,
           const vector_begin shiftBegin, bool shiftFlag)
{
    if (shiftFlag) [[likely]] {
        shiftfunc(x, aux, shiftBegin);
    } else [[unlikely]] {
        std::copy(x.begin(), x.end(), aux.begin());
    }
}

// TODO: move function into utils header and make it generic with constraints
// (and constexpr)
void multiplyVectorWithScalar(std::vector<double>& vec, double scalar)
{
    std::transform(
        std::execution::unseq, vec.begin(), vec.end(), vec.begin(),
        std::bind(std::multiplies<double>(), std::placeholders::_1, scalar));
}

void rotate(std::vector<double>& vec, const matrix_begin rotateBegin,
            bool rotateFlag)
{
    if (rotateFlag) [[likely]] {
        rotatefunc(vec, rotateBegin);
    }
}

void shiftRotateTransform(const std::vector<double>& x,
                          std::vector<double>& aux,
                          const vector_begin shiftBegin,
                          const matrix_begin rotateBegin, double shiftRate,
                          bool shiftFlag, bool rotateFlag)
{
    // shift rotate transform with shift rate != 1.0
    shift(x, aux, shiftBegin, shiftFlag);
    multiplyVectorWithScalar(aux, shiftRate);
    rotate(aux, rotateBegin, rotateFlag);
}

void shiftRotateTransform(const std::vector<double>& x,
                          std::vector<double>& aux,
                          const vector_begin shiftBegin,
                          const matrix_begin rotateBegin, bool shiftFlag,
                          bool rotateFlag)
{
    // shift rotate transform with shiftRate = 1.0
    shift(x, aux, shiftBegin, shiftFlag);
    rotate(aux, rotateBegin, rotateFlag);
}

void shiftRotateTransform(const std::vector<double>& x,
                          std::vector<double>& aux,
                          const std::vector<double>& shift,
                          const std::vector<std::vector<double>>& rotate,
                          double shiftRate, bool shiftFlag, bool rotateFlag)
{
    shiftRotateTransform(x, aux, shift.cbegin(), rotate.cbegin(), shiftRate,
                         shiftFlag, rotateFlag);
}

void shiftRotateTransform(const std::vector<double>& x,
                          std::vector<double>& aux,
                          const std::vector<double>& shift,
                          const std::vector<std::vector<double>>& rotate,
                          bool shiftFlag, bool rotateFlag)
{
    shiftRotateTransform(x, aux, shift.cbegin(), rotate.cbegin(), shiftFlag,
                         rotateFlag);
}

void applyPermutation(std::vector<double>& nums,
                      const std::vector<std::size_t>& indices)
{
    auto aux = std::vector<double>{};
    aux.reserve(indices.size());

    // assuming indices.size() == nums.size()
    std::transform(indices.begin(), indices.end(), std::back_inserter(aux),
                   [&](auto i) { return nums[i]; });

    // moving result of permutation to nums
    std::swap(aux, nums);
}

template <std::size_t Size>
double compositionFunctionCalculator(const std::vector<double>& x,
                                     const std::vector<double>& shift,
                                     const std::array<int, Size>& delta,
                                     const std::array<double, Size>& fit)
{
    const auto n = x.size();
    auto w_max = 0.0;
    auto w_sum = 0.0;
    std::array<double, Size> w{0.0};
    for (std::size_t i = 0; i < Size; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            const auto temp = x[j] - shift[i * n + j];
            w[i] += temp * temp;
        }

        // else will happen only when x is shift
        if (w[i] != 0.0) [[likely]] {
            w[i] = std::sqrt(1.0 / w[i]) *
                   std::exp(-w[i] / 2.0 / n / delta[i] / delta[i]);
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
        return std::accumulate(fit.begin(), fit.end(), 0.0,
                               [](auto f, auto elem) { return f + elem; }) /
               Size;
    }

    return std::inner_product(
        w.begin(), w.end(), fit.begin(), 0.0, std::plus<double>(),
        [=](auto w, auto fit) { return w / w_sum * fit; });
}

} // namespace

double do_ackley_func(const VectorRange& x)
{
    auto sum1 = 0.0;
    auto sum2 = 0.0;
    std::for_each(x.begin, x.end, [&](auto elem) {
        sum1 += elem * elem;
        sum2 += std::cos(2.0 * PI * elem);
    });
    const auto n = std::distance(x.begin, x.end);
    sum1 = -0.2 * std::sqrt(sum1 / n);
    sum2 /= n;
    return E - 20.0 * std::exp(sum1) - std::exp(sum2) + 20.0;
}

double ackley_func(const std::vector<double>& x, std::vector<double>& aux,
                   const std::vector<double>& shift,
                   const std::vector<std::vector<double>>& rotate,
                   bool shiftFlag, bool rotateFlag)
{
    // shift rate is 1.0
    shiftRotateTransform(x, aux, shift.cbegin(), rotate.cbegin(), shiftFlag,
                         rotateFlag);

    return do_ackley_func({aux.begin(), aux.end()});
}

double ackley_func(const VectorRange& x)
{
    // justShift(x, 1.0) // not needed

    return do_ackley_func(x);
}

double do_bent_cigar_func(const VectorRange& x)
{
    // assuming x.size() >= 1
    return std::reduce(std::next(x.begin), x.end, (*x.begin) * (*x.begin),
                       [](auto f, auto elem) {
                           return f + elem * elem * 1000000.0;
                       }); // 1000000.0 = std::pow(10.0, 6.0)
}

double bent_cigar_func(const VectorRange& x)
{
    // justShift(x, 1.0); // not needed

    return do_bent_cigar_func(x);
}

double
bent_cigar_func(const std::vector<double>& x, std::vector<double>& aux,
                const vector_begin shiftBegin, const matrix_begin rotateBegin,
                bool shiftFlag, bool rotateFlag)
{
    // shift rate is 1.0
    shiftRotateTransform(x, aux, shiftBegin, rotateBegin, shiftFlag,
                         rotateFlag);

    return do_bent_cigar_func({aux.begin(), aux.end()});
}

double bent_cigar_func(const std::vector<double>& x, std::vector<double>& aux,
                       const std::vector<double>& shift,
                       const std::vector<std::vector<double>>& rotate,
                       bool shiftFlag, bool rotateFlag)
{
    return bent_cigar_func(x, aux, shift.cbegin(), rotate.cbegin(), shiftFlag,
                           rotateFlag);
}

double
discus_func(const std::vector<double>& x, std::vector<double>& aux,
            const vector_begin shiftBegin, const matrix_begin rotateBegin,
            bool shiftFlag, bool rotateFlag)
{
    // shift rate is 1.0
    shiftRotateTransform(x, aux, shiftBegin, rotateBegin, shiftFlag,
                         rotateFlag);
    // 1000000.0 = std::pow(10.0, 6.0)
    // assuming x.size() >= 1
    return std::reduce(std::next(aux.begin()), aux.end(),
                       aux[0] * aux[0] * 1000000.0,
                       [](auto f, auto elem) { return f + elem * elem; });
}

double discus_func(const std::vector<double>& x, std::vector<double>& aux,
                   const std::vector<double>& shift,
                   const std::vector<std::vector<double>>& rotate,
                   bool shiftFlag, bool rotateFlag)
{
    return discus_func(x, aux, shift.cbegin(), rotate.cbegin(), shiftFlag,
                       rotateFlag);
}

double
ellips_func(const std::vector<double>& x, std::vector<double>& aux,
            const vector_begin shiftBegin, const matrix_begin rotateBegin,
            bool shiftFlag, bool rotateFlag)
{
    // shift rate is 1.0
    shiftRotateTransform(x, aux, shiftBegin, rotateBegin, shiftFlag,
                         rotateFlag);
    // here we need to accumulate because reduce does not maintain order
    return std::accumulate(
        aux.begin(), aux.end(), 0.0,
        [i = 0, n = aux.size() - 1.0](auto f, auto elem) mutable {
            return std::move(f) + std::pow(10.0, 6.0 * (i++) / n) * elem * elem;
        });
}

double ellips_func(const std::vector<double>& x, std::vector<double>& aux,
                   const std::vector<double>& shift,
                   const std::vector<std::vector<double>>& rotate,
                   bool shiftFlag, bool rotateFlag)
{
    return ellips_func(x, aux, shift.cbegin(), rotate.cbegin(), shiftFlag,
                       rotateFlag);
}

double
escaffer6_func(const std::vector<double>& x, std::vector<double>& aux,
               const vector_begin shiftBegin, const matrix_begin rotateBegin,
               bool shiftFlag, bool rotateFlag)
{
    // shift rate is 1.0
    shiftRotateTransform(x, aux, shiftBegin, rotateBegin, shiftFlag,
                         rotateFlag);

    const auto n = aux.size();
    auto f = 0.0;
    for (std::size_t i = 0; i < n - 1; ++i) {
        const auto xi = aux[i] * aux[i];
        const auto xinext = aux[i + 1] * aux[i + 1];
        const auto temp1 = std::sin(std::sqrt(xi + xinext));
        const auto temp2 = 1.0 + 0.001 * (xi + xinext);
        f += 0.5 + (temp1 * temp1 - 0.5) / (temp2 * temp2);
    }

    const auto first = aux[0] * aux[0];
    // assuming x.size() is not 0
    const auto last = aux[n - 1] * aux[n - 1];
    const auto temp1 = std::sin(std::sqrt(last + first));
    const auto temp2 = 1.0 + 0.001 * (last + first);
    f += 0.5 + (temp1 * temp1 - 0.5) / (temp2 * temp2);
    //
    return f;
}

double escaffer6_func(const std::vector<double>& x, std::vector<double>& aux,
                      const std::vector<double>& shift,
                      const std::vector<std::vector<double>>& rotate,
                      bool shiftFlag, bool rotateFlag)
{
    return escaffer6_func(x, aux, shift.cbegin(), rotate.cbegin(), shiftFlag,
                          rotateFlag);
}

double
griewank_func(const std::vector<double>& x, std::vector<double>& aux,
              const vector_begin shiftBegin, const matrix_begin rotateBegin,
              bool shiftFlag, bool rotateFlag)
{
    shiftRotateTransform(x, aux, shiftBegin, rotateBegin, 600.0 / 100.0,
                         shiftFlag, rotateFlag);

    auto s = 0.0;
    auto p = 1.0;
    for (std::size_t i = 0, n = aux.size(); i < n; ++i) {
        s += aux[i] * aux[i];
        p *= std::cos(aux[i] / std::sqrt(1.0 + i));
    }
    return 1.0 + s / 4000.0 - p;
}

double griewank_func(const std::vector<double>& x, std::vector<double>& aux,
                     const std::vector<double>& shift,
                     const std::vector<std::vector<double>>& rotate,
                     bool shiftFlag, bool rotateFlag)
{
    return griewank_func(x, aux, shift.cbegin(), rotate.cbegin(), shiftFlag,
                         rotateFlag);
}

double do_grie_rosen_func(const VectorRange& x)
{
    // assuming x.size() is not 0
    const auto first = *x.begin + 1.0;
    const auto end = std::prev(x.end);
    auto f = 0.0;

    for (auto it = x.begin; it != end;) {
        const auto current = *it + 1.0;
        const auto next = *(++it) + 1.0;
        const auto temp1 = current * current - next;
        const auto temp2 = current - 1.0;
        const auto temp = 100.0 * temp1 * temp1 + temp2 * temp2;
        f += (temp * temp) / 4000.0 - std::cos(temp) + 1.0;
    }

    const auto last = *end + 1.0;
    const auto temp1 = last * last - first;
    const auto temp2 = last - 1.0;
    const auto temp = 100.0 * temp1 * temp1 + temp2 * temp2;
    f += (temp * temp) / 4000.0 - std::cos(temp) + 1.0;
    return f;
}

double grie_rosen_func(const std::vector<double>& x, std::vector<double>& aux,
                       const std::vector<double>& shift,
                       const std::vector<std::vector<double>>& rotate,
                       bool shiftFlag, bool rotateFlag)
{
    shiftRotateTransform(x, aux, shift.cbegin(), rotate.cbegin(), 5.0 / 100.0,
                         shiftFlag, rotateFlag);

    return do_grie_rosen_func({aux.begin(), aux.end()});
}

double grie_rosen_func(const VectorRange& x)
{
    justShift(x, 5.0 / 100.0);

    return do_grie_rosen_func(x);
}

double do_happycat_func(const VectorRange& x)
{
    auto sum_y = 0.0;
    auto r2 = 0.0;
    std::for_each(x.begin, x.end, [&](auto elem) {
        const auto aux = elem - 1.0;
        sum_y += aux;
        r2 += aux * aux;
    });

    const auto n = std::distance(x.begin, x.end);
    return std::pow(std::abs(r2 - n), 2 * 1.0 / 8.0) + (0.5 * r2 + sum_y) / n +
           0.5;
}

double happycat_func(const std::vector<double>& x, std::vector<double>& aux,
                     const std::vector<double>& shift,
                     const std::vector<std::vector<double>>& rotate,
                     bool shiftFlag, bool rotateFlag)
{
    shiftRotateTransform(x, aux, shift.cbegin(), rotate.cbegin(), 5.0 / 100.0,
                         shiftFlag, rotateFlag);

    return do_happycat_func({aux.begin(), aux.end()});
}

double happycat_func(const VectorRange& x)
{
    justShift(x, 5.0 / 100.0);

    return do_happycat_func(x);
}

double do_hgbat_func(const VectorRange& x)
{
    auto sum_y = 0.0;
    auto r2 = 0.0;
    std::for_each(x.begin, x.end, [&](auto elem) {
        const auto aux = elem - 1.0;
        sum_y += aux;
        r2 += aux * aux;
    });

    return std::pow(std::fabs(r2 * r2 - sum_y * sum_y), 2.0 * 1.0 / 4.0) +
           (0.5 * r2 + sum_y) / std::distance(x.begin, x.end) + 0.5;
}

double hgbat_func(const std::vector<double>& x, std::vector<double>& aux,
                  const vector_begin shiftBegin, const matrix_begin rotateBegin,
                  bool shiftFlag, bool rotateFlag)
{
    shiftRotateTransform(x, aux, shiftBegin, rotateBegin, 5.0 / 100.0,
                         shiftFlag, rotateFlag);

    return do_hgbat_func({aux.begin(), aux.end()});
}

double hgbat_func(const VectorRange& x)
{
    justShift(x, 5.0 / 100.0);

    return do_hgbat_func(x);
}

double hgbat_func(const std::vector<double>& x, std::vector<double>& aux,
                  const std::vector<double>& shift,
                  const std::vector<std::vector<double>>& rotate,
                  bool shiftFlag, bool rotateFlag)
{
    return hgbat_func(x, aux, shift.cbegin(), rotate.cbegin(), shiftFlag,
                      rotateFlag);
}

double
rosenbrock_func(const std::vector<double>& x, std::vector<double>& aux,
                const vector_begin shiftBegin, const matrix_begin rotateBegin,
                bool shiftFlag, bool rotateFlag)
{
    shiftRotateTransform(x, aux, shiftBegin, rotateBegin, 2.048 / 100.0,
                         shiftFlag, rotateFlag);
    auto f = 0.0;
    for (std::size_t i = 0, n = aux.size() - 1; i < n; ++i) {
        const auto temp = aux[i] + 1.0;
        const auto temp1 = temp * temp - aux[i + 1] - 1.0;
        f += 100.0 * temp1 * temp1 + aux[i] * aux[i];
    }
    return f;
}

double rosenbrock_func(const std::vector<double>& x, std::vector<double>& aux,
                       const std::vector<double>& shift,
                       const std::vector<std::vector<double>>& rotate,
                       bool shiftFlag, bool rotateFlag)
{
    return rosenbrock_func(x, aux, shift.cbegin(), rotate.cbegin(), shiftFlag,
                           rotateFlag);
}

double do_rastrigin_func(const VectorRange& x)
{
    return std::reduce(x.begin, x.end, std::distance(x.begin, x.end) * 10.0,
                       [=](auto f, auto elem) {
                           return f + elem * elem -
                                  10.0 * std::cos(2.0 * PI * elem);
                       });
}

double
rastrigin_func(const std::vector<double>& x, std::vector<double>& aux,
               const vector_begin shiftBegin, const matrix_begin rotateBegin,
               bool shiftFlag, bool rotateFlag)
{
    shiftRotateTransform(x, aux, shiftBegin, rotateBegin, 5.12 / 100.0,
                         shiftFlag, rotateFlag);

    return do_rastrigin_func({aux.begin(), aux.end()});
}

double rastrigin_func(const VectorRange& x)
{
    justShift(x, 5.12 / 100.0);

    return do_rastrigin_func(x);
}

double rastrigin_func(const std::vector<double>& x, std::vector<double>& aux,
                      const std::vector<double>& shift,
                      const std::vector<std::vector<double>>& rotate,
                      bool shiftFlag, bool rotateFlag)
{
    return rastrigin_func(x, aux, shift.cbegin(), rotate.cbegin(), shiftFlag,
                          rotateFlag);
}

double do_schwefel_func(const VectorRange& x)
{
    const auto n = std::distance(x.begin, x.end);
    return 4.189828872724338e+002 * n +
           std::reduce(x.begin, x.end, 0.0, [=](auto f, auto elem) {
               const auto xi = elem + 4.209687462275036e+002;
               if (xi > 500.0) {
                   const auto temp1 =
                       (500.0 - std::fmod(xi, 500)) *
                       std::sin(std::sqrt(500.0 - std::fmod(xi, 500)));
                   const auto temp2 = (xi - 500.0) / 100.0;
                   return f - temp1 + temp2 * temp2 / n;
               } else if (xi < -500.0) {
                   const auto temp1 =
                       (-500.0 + std::fmod(std::fabs(xi), 500)) *
                       std::sin(
                           std::sqrt(500.0 - std::fmod(std::fabs(xi), 500)));
                   const auto temp2 = (xi + 500.0) / 100.0;
                   return f - temp1 + temp2 * temp2 / n;
               } // else
               return f - xi * std::sin(std::sqrt(std::fabs(xi)));
           });
}

double
schwefel_func(const std::vector<double>& x, std::vector<double>& aux,
              const vector_begin shiftBegin, const matrix_begin rotateBegin,
              bool shiftFlag, bool rotateFlag)
{
    shiftRotateTransform(x, aux, shiftBegin, rotateBegin, 1000.0 / 100.0,
                         shiftFlag, rotateFlag);

    return do_schwefel_func({aux.begin(), aux.end()});
}

double schwefel_func(const VectorRange& x)
{
    justShift(x, 1000.0 / 100.0);

    return do_schwefel_func(x);
}

double schwefel_func(const std::vector<double>& x, std::vector<double>& aux,
                     const std::vector<double>& shift,
                     const std::vector<std::vector<double>>& rotate,
                     bool shiftFlag, bool rotateFlag)
{
    return schwefel_func(x, aux, shift.cbegin(), rotate.cbegin(), shiftFlag,
                         rotateFlag);
}

double do_schaffer_F7_func(const VectorRange& x)
{
    const auto end = std::prev(x.end);
    auto f = 0.0;

    for (auto it = x.begin; it != end;) {
        const auto temp1 = *it;
        std::advance(it, 1);
        const auto temp2 = *it;
        const auto si = std::sqrt(temp1 * temp1 + temp2 * temp2);
        const auto temp = std::sin(50.0 * std::pow(si, 0.2));
        const auto sqrtsi = std::sqrt(si);
        f += sqrtsi + sqrtsi * temp * temp;
    }

    const auto n = std::distance(x.begin, x.end) - 1;
    f = f * f / n / n;
    return f;
}

double schaffer_F7_func(const std::vector<double>& x, std::vector<double>& aux,
                        const std::vector<double>& shift,
                        const std::vector<std::vector<double>>& rotate,
                        bool shiftFlag, [[maybe_unused]] bool rotateFlag)
{
    // schaffer_F7_func is wrong, it's not rotated
    // shift rate is 1.0
    shiftRotateTransform(x, aux, shift, rotate, shiftFlag, false);

    return do_schaffer_F7_func({aux.begin(), aux.end()});
}

double schaffer_F7_func(const VectorRange& x)
{
    // justShift(x, 1.0) // not needed

    return do_schaffer_F7_func(x);
}

double
step_rastrigin_func(const std::vector<double>& x, std::vector<double>& aux,
                    const std::vector<double>& shift,
                    const std::vector<std::vector<double>>& rotate,
                    bool shiftFlag, bool rotateFlag)
{
    return rastrigin_func(x, aux, shift, rotate, shiftFlag, rotateFlag);
}

double levy_func(const std::vector<double>& x, std::vector<double>& aux,
                 const std::vector<double>& shift,
                 const std::vector<std::vector<double>>& rotate, bool shiftFlag,
                 bool rotateFlag)
{
    // min is 1.49966e-32, close to 0
    // shift rate is 1.0
    shiftRotateTransform(x, aux, shift, rotate, shiftFlag, rotateFlag);

    // Correct is (elem - 1.0), but it does not provide a minumum close to 0
    const auto w = [](auto elem) { return 1.0 + (elem - 0.0) / 4.0; };

    const auto term1 = std::sin(PI * w(aux[0]));
    const auto term2 = std::accumulate(
        aux.begin(), std::prev(aux.end()), 0.0, [&w](auto f, auto elem) {
            const auto wi = w(elem);
            const auto temp = std::sin(PI * wi + 1.0);
            return f + (wi - 1.0) * (wi - 1.0) * (1.0 + 10.0 * temp * temp);
        });
    const auto last = w(*std::prev(aux.end()));
    const auto temp = std::sin(2.0 * PI * last);
    const auto term3 = (last - 1.0) * (last - 1.0) * (1.0 + temp * temp);
    return term1 * term1 + term2 + term3;
}

double zakharov_func(const std::vector<double>& x, std::vector<double>& aux,
                     const std::vector<double>& shift,
                     const std::vector<std::vector<double>>& rotate,
                     bool shiftFlag, bool rotateFlag)
{
    // shift rate is 1.0
    shiftRotateTransform(x, aux, shift, rotate, shiftFlag, rotateFlag);

    auto sum1 = 0.0;
    auto sum2 = 0.0;
    for (std::size_t i = 0, n = aux.size(); i < n; ++i) {
        sum1 += aux[i] * aux[i];
        sum2 += 0.5 * i * aux[i];
    }
    return sum1 + sum2 * sum2 + sum2 * sum2 * sum2 * sum2;
}

double do_katsuura_func(const VectorRange& x)
{
    const auto size = std::distance(x.begin, x.end);
    const auto temp3 = std::pow(static_cast<double>(size), 1.2);
    auto f = 1.0;
    for (auto i = 0; i < size; ++i) {
        auto temp = 0.0;
        const auto xi = *std::next(x.begin, i);
        for (auto j = 1; j < 33; ++j) {
            const auto temp1 = std::pow(2.0, j);
            const auto temp2 = temp1 * xi;
            temp += std::fabs(temp2 - std::floor(temp2 + 0.5)) / temp1;
        }
        f *= std::pow(1.0 + (i + 1) * temp, 10.0 / temp3);
    }
    const auto temp1 = 10.0 / size / size;
    return f * temp1 - temp1;
}

double katsuura_func(const std::vector<double>& x, std::vector<double>& aux,
                     const std::vector<double>& shift,
                     const std::vector<std::vector<double>>& rotate,
                     bool shiftFlag, bool rotateFlag)
{
    shiftRotateTransform(x, aux, shift, rotate, 5.0 / 100.0, shiftFlag,
                         rotateFlag);

    return do_katsuura_func({aux.begin(), aux.end()});
}

double katsuura_func(const VectorRange& x)
{
    justShift(x, 5.0 / 100.0);

    return do_katsuura_func(x);
}

double
hf01(const std::vector<double>& x, std::vector<double>& aux,
     const std::vector<double>& shift,
     const std::vector<std::vector<double>>& rotate,
     const std::vector<std::size_t>& indices, bool shiftFlag, bool rotateFlag)
{
    // [0.4, 0.4, 0.2]
    // shift rate is 1.0
    shiftRotateTransform(x, aux, shift, rotate, shiftFlag, rotateFlag);
    applyPermutation(aux, indices);

    const auto limit = std::ceil(0.4 * aux.size());
    const auto margin_1 = std::next(aux.begin(), limit); // 0.4
    const auto margin_2 = std::next(margin_1, limit);    // another 0.4

    const auto range1 = VectorRange{aux.begin(), margin_1};
    const auto range2 = VectorRange{margin_1, margin_2};
    const auto range3 = VectorRange{margin_2, aux.end()};

    // these methods do not shift and rotate
    return bent_cigar_func(range1) + hgbat_func(range2) +
           rastrigin_func(range3);
}

double
hf02(const std::vector<double>& x, std::vector<double>& aux,
     const std::vector<double>& shift,
     const std::vector<std::vector<double>>& rotate,
     const std::vector<std::size_t>& indices, bool shiftFlag, bool rotateFlag)
{
    // [0.1, 0.2, 0.2, 0.2, 0.1, 0.2]
    // shift rate is 1.0
    shiftRotateTransform(x, aux, shift, rotate, shiftFlag, rotateFlag);
    applyPermutation(aux, indices);

    const auto limit1 = std::ceil(0.1 * aux.size());
    const auto limit2 = std::ceil(0.2 * aux.size());

    const auto margin1 = std::next(aux.begin(), limit1); // 0.1
    const auto margin2 = std::next(margin1, limit2);     // 0.2
    const auto margin3 = std::next(margin2, limit2);     // 0.2
    const auto margin4 = std::next(margin3, limit2);     // 0.2
    const auto margin5 = std::next(margin4, limit1);     // 0.1
    // remaining: 0.2

    const auto range1 = VectorRange{aux.begin(), margin1};
    const auto range2 = VectorRange{margin1, margin2};
    const auto range3 = VectorRange{margin2, margin3};
    const auto range4 = VectorRange{margin3, margin4};
    const auto range5 = VectorRange{margin4, margin5};
    // const auto range6 = VectorRange{margin5, x.end()}; // CORRECT
    auto copy = std::vector<double>{aux.begin(), aux.begin() + 2};
    const auto range6 = VectorRange{copy.begin(), copy.end()}; // WRONG

    return hgbat_func(range1) + katsuura_func(range2) + ackley_func(range3) +
           rastrigin_func(range4) + schwefel_func(range5) +
           schaffer_F7_func(range6);
}

double
hf03(const std::vector<double>& x, std::vector<double>& aux,
     const std::vector<double>& shift,
     const std::vector<std::vector<double>>& rotate,
     const std::vector<std::size_t>& indices, bool shiftFlag, bool rotateFlag)
{
    //  [0.3, 0.2, 0.2, 0.1, 0.2]
    // shift rate is 1.0
    shiftRotateTransform(x, aux, shift, rotate, shiftFlag, rotateFlag);
    applyPermutation(aux, indices);

    const auto limit1 = std::ceil(0.1 * aux.size());
    const auto limit2 = std::ceil(0.2 * aux.size());
    const auto limit3 = std::ceil(0.3 * aux.size());

    const auto margin1 = std::next(aux.begin(), limit3); // 0.3
    const auto margin2 = std::next(margin1, limit2);     // 0.2
    const auto margin3 = std::next(margin2, limit2);     // 0.2
    const auto margin4 = std::next(margin3, limit1);     // 0.1
    // remaining: 0.2

    const auto range1 = VectorRange{aux.begin(), margin1};
    const auto range2 = VectorRange{margin1, margin2};
    const auto range3 = VectorRange{margin2, margin3};
    const auto range4 = VectorRange{margin3, margin4};
    const auto range5 = VectorRange{margin4, aux.end()};

    return katsuura_func(range1) + happycat_func(range2) +
           grie_rosen_func(range3) + schwefel_func(range4) +
           ackley_func(range5);
}

double cf01(const std::vector<double>& x, std::vector<double>& aux,
            const std::vector<double>& shift,
            const std::vector<std::vector<double>>& rotate, bool rotateFlag)
{
    // asumming that shift.size and rotate.size is (at least) 5 * x.size
    const auto size = x.size();
    const auto shift_margin_1 = std::next(shift.cbegin(), size);
    const auto shift_margin_2 = std::next(shift_margin_1, size);
    const auto shift_margin_3 = std::next(shift_margin_2, size);
    const auto shift_margin_4 = std::next(shift_margin_3, size);

    const auto rotate_margin_1 = std::next(rotate.cbegin(), size);
    const auto rotate_margin_2 = std::next(rotate_margin_1, size);
    const auto rotate_margin_3 = std::next(rotate_margin_2, size);
    const auto rotate_margin_4 = std::next(rotate_margin_3, size);

    constexpr auto N = 5;
    // fit is function result * lambda + bias
    // lambda is 1, 1e-6, 1e-6, 1e-6, 1e-6
    // lambda in their implementation is 1, 1e-6, 1e-26, 1e-6, 1e-6
    // bias is 0, 200, 300, 100, 400

    const std::array<double, N> fit{
        rosenbrock_func(x, aux, shift.cbegin(), rotate.cbegin(), true,
                        rotateFlag),
        ellips_func(x, aux, shift_margin_1, rotate_margin_1, true, rotateFlag) *
                1e-6 +
            200,
        bent_cigar_func(x, aux, shift_margin_2, rotate_margin_2, true,
                        rotateFlag) *
                1e-26 +
            300,
        discus_func(x, aux, shift_margin_3, rotate_margin_3, true, rotateFlag) *
                1e-6 +
            100,
        ellips_func(x, aux, shift_margin_4, rotate_margin_4, true, false) *
                1e-6 +
            400, // ?? why false
    };

    const std::array<int, N> delta{10, 20, 30, 40, 50};
    return compositionFunctionCalculator<N>(x, shift, delta, fit);
}

double cf02(const std::vector<double>& x, std::vector<double>& aux,
            const std::vector<double>& shift,
            const std::vector<std::vector<double>>& rotate, bool rotateFlag)
{
    // asumming that shift.size and rotate.size is (at least) 3 * x.size
    const auto size = x.size();
    const auto shift_margin_1 = std::next(shift.begin(), size);
    const auto shift_margin_2 = std::next(shift_margin_1, size);

    const auto rotate_margin_1 = std::next(rotate.begin(), size);
    const auto rotate_margin_2 = std::next(rotate_margin_1, size);

    constexpr auto N = 3;
    // fit is function result * lambda + bias
    // lambda is 1, 1, 1
    // bias is 0, 200, 100
    const std::array<double, N> fit{
        schwefel_func(x, aux, shift.cbegin(), rotate.cbegin(), true,
                      false), // ?? why false
        rastrigin_func(x, aux, shift_margin_1, rotate_margin_1, true,
                       rotateFlag) +
            200,
        hgbat_func(x, aux, shift_margin_2, rotate_margin_2, true, rotateFlag) +
            100,
    };
    const std::array<int, N> delta{20, 10, 10};
    return compositionFunctionCalculator<N>(x, shift, delta, fit);
}

double cf03(const std::vector<double>& x, std::vector<double>& aux,
            const std::vector<double>& shift,
            const std::vector<std::vector<double>>& rotate, bool rotateFlag)
{
    // asumming that shift.size and rotate.size is (at least) 5 * x.size
    const auto size = x.size();
    const auto shift_margin_1 = std::next(shift.cbegin(), size);
    const auto shift_margin_2 = std::next(shift_margin_1, size);
    const auto shift_margin_3 = std::next(shift_margin_2, size);
    const auto shift_margin_4 = std::next(shift_margin_3, size);

    const auto rotate_margin_1 = std::next(rotate.cbegin(), size);
    const auto rotate_margin_2 = std::next(rotate_margin_1, size);
    const auto rotate_margin_3 = std::next(rotate_margin_2, size);
    const auto rotate_margin_4 = std::next(rotate_margin_3, size);

    constexpr auto N = 5;
    // fit is function result * lambda + bias
    // lambda = [1e-26, 10, 1e-6, 10, 5e-4]
    // their lambda is [10000/2e+7, 1, 1000 / 100, 1,  10000 / 1e+3]
    // bias = [0, 200, 300, 400, 200]
    const std::array<double, N> fit{
        escaffer6_func(x, aux, shift.cbegin(), rotate.cbegin(), true,
                       rotateFlag) *
            (10000.0 / 2e+7),
        schwefel_func(x, aux, shift_margin_1, rotate_margin_1, true,
                      rotateFlag) *
                1.0 +
            200,
        griewank_func(x, aux, shift_margin_2, rotate_margin_2, true,
                      rotateFlag) *
                (1000 / 100) +
            300,
        rosenbrock_func(x, aux, shift_margin_3, rotate_margin_3, true,
                        rotateFlag) *
                1 +
            400,
        rastrigin_func(x, aux, shift_margin_4, rotate_margin_4, true,
                       rotateFlag) *
                (10000 / 1e+3) +
            200,

    };
    const std::array<int, N> delta{20, 20, 30, 30, 20};
    return compositionFunctionCalculator<N>(x, shift, delta, fit);
}

double cf04(const std::vector<double>& x, std::vector<double>& aux,
            const std::vector<double>& shift,
            const std::vector<std::vector<double>>& rotate, bool rotateFlag)
{
    // asumming that shift.size and rotate.size is (at least) 6 * x.size
    const auto size = x.size();
    const auto shift_margin_1 = std::next(shift.cbegin(), size);
    const auto shift_margin_2 = std::next(shift_margin_1, size);
    const auto shift_margin_3 = std::next(shift_margin_2, size);
    const auto shift_margin_4 = std::next(shift_margin_3, size);
    const auto shift_margin_5 = std::next(shift_margin_4, size);

    const auto rotate_margin_1 = std::next(rotate.cbegin(), size);
    const auto rotate_margin_2 = std::next(rotate_margin_1, size);
    const auto rotate_margin_3 = std::next(rotate_margin_2, size);
    const auto rotate_margin_4 = std::next(rotate_margin_3, size);
    const auto rotate_margin_5 = std::next(rotate_margin_4, size);

    constexpr auto N = 6;
    // fit is function result * lambda + bias
    // lambda = [10, 10, 2.5, 1e−26, 1e-6, 5e-4]
    // bias = [0, 300, 500, 100, 400, 200]
    const std::array<double, N> fit{
        hgbat_func(x, aux, shift.cbegin(), rotate.cbegin(), true, rotateFlag) *
            10,
        rastrigin_func(x, aux, shift_margin_1, rotate_margin_1, true,
                       rotateFlag) *
                10 +
            300,
        schwefel_func(x, aux, shift_margin_2, rotate_margin_2, true,
                      rotateFlag) *
                2.5 +
            500,
        bent_cigar_func(x, aux, shift_margin_3, rotate_margin_3, true,
                        rotateFlag) *
                1e-26 +
            100,
        ellips_func(x, aux, shift_margin_4, rotate_margin_4, true, rotateFlag) *
                1e-6 +
            400,
        escaffer6_func(x, aux, shift_margin_5, rotate_margin_5, true,
                       rotateFlag) *
                5e-4 +
            200,

    };
    const std::array<int, N> delta{10, 20, 30, 40, 50, 60};
    return compositionFunctionCalculator<N>(x, shift, delta, fit);
}

int sanity_check()
{
    auto x = std::vector<double>(10, 0.0);
    auto aux = std::vector<double>(10, 0.0);
    const auto shift = std::vector<double>(60, 0.0);
    const auto rotate = std::vector<std::vector<double>>{};
    const auto indices = std::vector<std::size_t>{0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    std::cout << "zakharov_func = "
              << zakharov_func(x, aux, shift, rotate, false, false)
              << std::endl;
    std::fill(x.begin(), x.end(), 0.0);
    std::cout << "rosenbrock_func = "
              << rosenbrock_func(x, aux, shift, rotate, false, false)
              << std::endl;
    std::fill(x.begin(), x.end(), 0.0);
    std::cout << "schaffer_F7_func = "
              << schaffer_F7_func(x, aux, shift, rotate, false, false)
              << std::endl;
    std::fill(x.begin(), x.end(), 0.0);
    std::cout << "rastrigin_func = "
              << rastrigin_func(x, aux, shift, rotate, false, false)
              << std::endl;
    std::fill(x.begin(), x.end(), 0.0);
    std::cout << "levy_func = "
              << levy_func(x, aux, shift, rotate, false, false) << std::endl;
    std::fill(x.begin(), x.end(), 0.0);
    std::cout << "hf01 = " << hf01(x, aux, shift, rotate, indices, false, false)
              << std::endl;
    std::fill(x.begin(), x.end(), 0.0);
    std::cout << "hf02 = " << hf02(x, aux, shift, rotate, indices, false, false)
              << std::endl;
    std::fill(x.begin(), x.end(), 0.0);
    std::cout << "hf03 = " << hf03(x, aux, shift, rotate, indices, false, false)
              << std::endl;

    // Composition functions do not seem to have their minimum in 0
    std::fill(x.begin(), x.end(), 0.0);
    std::cout << "cf01 = " << cf01(x, aux, shift, rotate, false) << std::endl;
    std::fill(x.begin(), x.end(), 0.0);
    std::cout << "cf02 = " << cf02(x, aux, shift, rotate, false) << std::endl;
    std::fill(x.begin(), x.end(), 0.0);
    std::cout << "cf03 = " << cf03(x, aux, shift, rotate, false) << std::endl;
    std::fill(x.begin(), x.end(), 0.0);
    std::cout << "cf04 = " << cf04(x, aux, shift, rotate, false) << std::endl;
    return 1;
}

} // namespace cec22
