#include "Cec22.h"

#include <algorithm>
#include <cmath>
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

void shiftfunc(std::vector<double>& x, const vector_begin shiftBegin)
{
    // assuming x.size() == shift.size()
    std::transform(x.begin(), x.end(), shiftBegin, x.begin(),
                   std::minus<double>());
    // TODO: Check what nethods can be vectorized
}

void rotatefunc(std::vector<double>& x, std::vector<double>& aux,
                const matrix_begin rotateBegin)
{
    // assuming x.size() == rotate.size() == rotate[0].size() <= aux.size()
    const auto n = x.size();
    const auto auxBegin = aux.begin();
    std::fill(auxBegin, std::next(auxBegin, n), 0.0); // setting aux with 0.0
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            aux[i] +=
                x[j] * (*std::next(rotateBegin,
                                   i))[j]; // rotate[i][j] but with iterators
        }
    }
    x.swap(aux);
}

// const VectorRange is misleading. Value of iterators can change, the position
// won't
void justShift(const VectorRange& x, double shiftRate)
{
    std::transform(
        x.begin, x.end, x.begin,
        std::bind(std::multiplies<double>(), std::placeholders::_1, shiftRate));
}

void shiftRotateTransform(std::vector<double>& x, std::vector<double>& aux,
                          const vector_begin shiftBegin,
                          const matrix_begin rotateBegin, double shiftRate,
                          bool shiftFlag, bool rotateFlag)
{
    // template bool and if constexpr
    if (shiftFlag) [[likely]] {
        shiftfunc(x, shiftBegin);
    }

    // should we check if shift rate is 1.0? Or should we provide a different
    // overload without shift rate?
    std::transform(
        x.begin(), x.end(), x.begin(),
        std::bind(std::multiplies<double>(), std::placeholders::_1, shiftRate));

    if (rotateFlag) [[likely]] {
        rotatefunc(x, aux, rotateBegin);
    }
}

void shiftRotateTransform(std::vector<double>& x, std::vector<double>& aux,
                          const std::vector<double>& shift,
                          const std::vector<std::vector<double>>& rotate,
                          double shiftRate, bool shiftFlag, bool rotateFlag)
{
    shiftRotateTransform(x, aux, shift.cbegin(), rotate.cbegin(), shiftRate,
                         shiftFlag, rotateFlag);
}

void applyPermutation(std::vector<double>& nums, std::vector<double>& aux,
                      const std::vector<std::size_t>& indices)
{
    // assuming indices.size() == nums.size() == aux.size()
    std::transform(indices.begin(), indices.end(), aux.begin(),
                   [&](auto i) { return nums[i]; });
    nums.swap(aux);
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

double ackley_func(std::vector<double>& x, std::vector<double>& aux,
                   const std::vector<double>& shift,
                   const std::vector<std::vector<double>>& rotate,
                   bool shiftFlag, bool rotateFlag)
{
    shiftRotateTransform(x, aux, shift.cbegin(), rotate.cbegin(), 1.0,
                         shiftFlag, rotateFlag);

    return do_ackley_func({x.begin(), x.end()});
}

double ackley_func(const VectorRange& x)
{
    // justShift(x, 1.0) // not needed

    return do_ackley_func(x);
}

double do_bent_cigar_func(const VectorRange& x)
{
    // assuming x.size() >= 1
    // TODO: first
    return std::accumulate(std::next(x.begin), x.end, (*x.begin) * (*x.begin),
                           [](auto f, auto elem) {
                               return f + elem * elem * 1000000.0;
                           }); // 1000000.0 = std::pow(10.0, 6.0)
    // TODO: accumulate vs reduce?
}

double bent_cigar_func(const VectorRange& x)
{
    // justShift(x, 1.0); // not needed

    return do_bent_cigar_func(x);
}

double
bent_cigar_func(std::vector<double>& x, std::vector<double>& aux,
                const vector_begin shiftBegin, const matrix_begin rotateBegin,
                bool shiftFlag, bool rotateFlag)
{
    shiftRotateTransform(x, aux, shiftBegin, rotateBegin, 1.0, shiftFlag,
                         rotateFlag);

    return do_bent_cigar_func({x.begin(), x.end()});
}

double bent_cigar_func(std::vector<double>& x, std::vector<double>& aux,
                       const std::vector<double>& shift,
                       const std::vector<std::vector<double>>& rotate,
                       bool shiftFlag, bool rotateFlag)
{
    return bent_cigar_func(x, aux, shift.cbegin(), rotate.cbegin(), shiftFlag,
                           rotateFlag);
}

double
discus_func(std::vector<double>& x, std::vector<double>& aux,
            const vector_begin shiftBegin, const matrix_begin rotateBegin,
            bool shiftFlag, bool rotateFlag)
{
    shiftRotateTransform(x, aux, shiftBegin, rotateBegin, 1.0, shiftFlag,
                         rotateFlag);
    // 1000000.0 = std::pow(10.0, 6.0)
    // assuming x.size() >= 1
    return std::accumulate(std::next(x.begin()), x.end(),
                           x[0] * x[0] * 1000000.0,
                           [](auto f, auto elem) { return f + elem * elem; });
}

double discus_func(std::vector<double>& x, std::vector<double>& aux,
                   const std::vector<double>& shift,
                   const std::vector<std::vector<double>>& rotate,
                   bool shiftFlag, bool rotateFlag)
{
    return discus_func(x, aux, shift.cbegin(), rotate.cbegin(), shiftFlag,
                       rotateFlag);
}

double
ellips_func(std::vector<double>& x, std::vector<double>& aux,
            const vector_begin shiftBegin, const matrix_begin rotateBegin,
            bool shiftFlag, bool rotateFlag)
{
    shiftRotateTransform(x, aux, shiftBegin, rotateBegin, 1.0, shiftFlag,
                         rotateFlag);

    auto i = 0;
    // here we need to accumulate because reduce does not maintain order
    return std::accumulate(
        std::next(x.begin()), x.end(), 0.0,
        [&i, n = x.size()](auto f, auto elem) {
            return std::move(f) + std::pow(10.0, 6.0 * (i++) / n) * elem * elem;
        });
}

double ellips_func(std::vector<double>& x, std::vector<double>& aux,
                   const std::vector<double>& shift,
                   const std::vector<std::vector<double>>& rotate,
                   bool shiftFlag, bool rotateFlag)
{
    return ellips_func(x, aux, shift.cbegin(), rotate.cbegin(), shiftFlag,
                       rotateFlag);
}

double
escaffer6_func(std::vector<double>& x, std::vector<double>& aux,
               const vector_begin shiftBegin, const matrix_begin rotateBegin,
               bool shiftFlag, bool rotateFlag)
{
    shiftRotateTransform(x, aux, shiftBegin, rotateBegin, 1.0, shiftFlag,
                         rotateFlag);

    const auto n = x.size();
    auto f = 0.0;

    for (std::size_t i = 0; i < n - 1; ++i) {
        const auto xi = x[i] * x[i];
        const auto xinext = x[i + 1] * x[i + 1];
        // TODO: test in godbolt xi&xinext against code without them
        const auto temp1 = std::sin(std::sqrt(xi + xinext));
        const auto temp2 = 1.0 + 0.001 * (xi + xinext);
        f += 0.5 + (temp1 * temp1 - 0.5) / (temp2 * temp2);
    }

    const auto first = x[0] * x[0];
    // assuming x.size() is not 0
    const auto last = x[n - 1] * x[n - 1];
    const auto temp1 = std::sin(std::sqrt(last + first));
    const auto temp2 = 1.0 + 0.001 * (last + first);
    f += 0.5 + (temp1 * temp1 - 0.5) / (temp2 * temp2);
    //
    return f;
}

double escaffer6_func(std::vector<double>& x, std::vector<double>& aux,
                      const std::vector<double>& shift,
                      const std::vector<std::vector<double>>& rotate,
                      bool shiftFlag, bool rotateFlag)
{
    return escaffer6_func(x, aux, shift.cbegin(), rotate.cbegin(), shiftFlag,
                          rotateFlag);
}

double
griewank_func(std::vector<double>& x, std::vector<double>& aux,
              const vector_begin shiftBegin, const matrix_begin rotateBegin,
              bool shiftFlag, bool rotateFlag)
{
    shiftRotateTransform(x, aux, shiftBegin, rotateBegin, 600.0 / 100.0,
                         shiftFlag, rotateFlag);

    // TODO: n
    auto s = 0.0;
    auto p = 1.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        s += x[i] * x[i];
        p *= std::cos(x[i] / std::sqrt(1.0 + i));
    }
    return 1.0 + s / 4000.0 - p;
}

double griewank_func(std::vector<double>& x, std::vector<double>& aux,
                     const std::vector<double>& shift,
                     const std::vector<std::vector<double>>& rotate,
                     bool shiftFlag, bool rotateFlag)
{
    return griewank_func(x, aux, shift.cbegin(), rotate.cbegin(), shiftFlag,
                         rotateFlag);
}

double do_grie_rosen_func(const VectorRange& x)
{
    // TODO: Check impl
    auto f = 0.0;
    auto& first = *x.begin;
    first += 1.0;

    for (auto it = x.begin; it != x.end;) {
        const auto current = *it;
        auto& next = *(++it);
        next += 1.0;
        const auto temp1 = current * current - next;
        const auto temp2 = current - 1.0;
        const auto temp = 100.0 * temp1 * temp1 + temp2 * temp2;
        f += (temp * temp) / 4000.0 - std::cos(temp) + 1.0;
    }

    // assuming x.size() is not 0
    const auto last = *std::prev(x.end);
    const auto temp1 = last * last - first;
    const auto temp2 = last - 1.0;
    const auto temp = 100.0 * temp1 * temp1 + temp2 * temp2;
    f += (temp * temp) / 4000.0 - std::cos(temp) + 1.0;
    return f;
}

double grie_rosen_func(std::vector<double>& x, std::vector<double>& aux,
                       const std::vector<double>& shift,
                       const std::vector<std::vector<double>>& rotate,
                       bool shiftFlag, bool rotateFlag)
{
    shiftRotateTransform(x, aux, shift.cbegin(), rotate.cbegin(), 5.0 / 100.0,
                         shiftFlag, rotateFlag);

    return do_grie_rosen_func({x.begin(), x.end()});
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

double happycat_func(std::vector<double>& x, std::vector<double>& aux,
                     const std::vector<double>& shift,
                     const std::vector<std::vector<double>>& rotate,
                     bool shiftFlag, bool rotateFlag)
{
    shiftRotateTransform(x, aux, shift.cbegin(), rotate.cbegin(), 5.0 / 100.0,
                         shiftFlag, rotateFlag);

    return do_happycat_func({x.begin(), x.end()});
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

double hgbat_func(std::vector<double>& x, std::vector<double>& aux,
                  const vector_begin shiftBegin, const matrix_begin rotateBegin,
                  bool shiftFlag, bool rotateFlag)
{
    shiftRotateTransform(x, aux, shiftBegin, rotateBegin, 5.0 / 100.0,
                         shiftFlag, rotateFlag);

    return do_hgbat_func({x.begin(), x.end()});
}

double hgbat_func(const VectorRange& x)
{
    justShift(x, 5.0 / 100.0);

    return do_hgbat_func(x);
}

double hgbat_func(std::vector<double>& x, std::vector<double>& aux,
                  const std::vector<double>& shift,
                  const std::vector<std::vector<double>>& rotate,
                  bool shiftFlag, bool rotateFlag)
{
    return hgbat_func(x, aux, shift.cbegin(), rotate.cbegin(), shiftFlag,
                      rotateFlag);
}

double
rosenbrock_func(std::vector<double>& x, std::vector<double>& aux,
                const vector_begin shiftBegin, const matrix_begin rotateBegin,
                bool shiftFlag, bool rotateFlag)
{
    shiftRotateTransform(x, aux, shiftBegin, rotateBegin, 2.048 / 100.0,
                         shiftFlag, rotateFlag);

    auto f = 0.0;
    for (std::size_t i = 0; i < x.size() - 1; ++i) {
        const auto aux = x[i] + 1.0;
        const auto temp1 = aux * aux - x[i + 1] - 1.0;
        f += 100.0 * temp1 * temp1 + x[i] * x[i];
    }
    return f;
}

double rosenbrock_func(std::vector<double>& x, std::vector<double>& aux,
                       const std::vector<double>& shift,
                       const std::vector<std::vector<double>>& rotate,
                       bool shiftFlag, bool rotateFlag)
{
    return rosenbrock_func(x, aux, shift.cbegin(), rotate.cbegin(), shiftFlag,
                           rotateFlag);
}

double do_rastrigin_func(const VectorRange& x)
{
    return std::accumulate(
        x.begin, x.end, std::distance(x.begin, x.end), [=](auto f, auto elem) {
            return f + elem * elem - 10.0 * std::cos(2.0 * PI * elem);
        });
}

double
rastrigin_func(std::vector<double>& x, std::vector<double>& aux,
               const vector_begin shiftBegin, const matrix_begin rotateBegin,
               bool shiftFlag, bool rotateFlag)
{
    shiftRotateTransform(x, aux, shiftBegin, rotateBegin, 5.12 / 100.0,
                         shiftFlag, rotateFlag);

    return do_rastrigin_func({x.begin(), x.end()});
}

double rastrigin_func(const VectorRange& x)
{
    justShift(x, 5.12 / 100.0);

    return do_rastrigin_func(x);
}

double rastrigin_func(std::vector<double>& x, std::vector<double>& aux,
                      const std::vector<double>& shift,
                      const std::vector<std::vector<double>>& rotate,
                      bool shiftFlag, bool rotateFlag)
{
    return rastrigin_func(x, aux, shift.cbegin(), rotate.cbegin(), shiftFlag,
                          rotateFlag);
}

double do_schwefel_func(const VectorRange& x)
{
    return std::accumulate(
        x.begin, x.end, 0.0,
        [n = std::distance(x.begin, x.end)](auto f, auto elem) {
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
schwefel_func(std::vector<double>& x, std::vector<double>& aux,
              const vector_begin shiftBegin, const matrix_begin rotateBegin,
              bool shiftFlag, bool rotateFlag)
{
    shiftRotateTransform(x, aux, shiftBegin, rotateBegin, 1000.0 / 100.0,
                         shiftFlag, rotateFlag);

    return do_schwefel_func({x.begin(), x.end()});
}

double schwefel_func(const VectorRange& x)
{
    justShift(x, 1000.0 / 100.0);

    return do_schwefel_func(x);
}

double schwefel_func(std::vector<double>& x, std::vector<double>& aux,
                     const std::vector<double>& shift,
                     const std::vector<std::vector<double>>& rotate,
                     bool shiftFlag, bool rotateFlag)
{
    return schwefel_func(x, aux, shift.cbegin(), rotate.cbegin(), shiftFlag,
                         rotateFlag);
}

double do_schaffer_F7_func(const VectorRange& x)
{
    // shaffer_F7_func is tottaly different from original, which is broken
    // source:
    // https://www.sciencedirect.com/topics/computer-science/benchmark-function
    // TODO: check if results are ok and minimum is zero. If not, use the
    // provided version

    auto f = 0.0;

    for (auto it = x.begin; it != x.end;) {
        const auto temp1 = *it;
        std::advance(it, 1); // TODO: vs ++it
        const auto temp2 = *it;

        const auto si = std::sqrt(temp1 * temp1 + temp2 * temp2);
        const auto temp = std::sin(50.0 * std::pow(si, 0.2));
        f += si + si * temp * temp;
    }

    const auto n = std::distance(x.begin, x.end);
    return f * f / n / n;
}

double schaffer_F7_func(std::vector<double>& x, std::vector<double>& aux,
                        const std::vector<double>& shift,
                        const std::vector<std::vector<double>>& rotate,
                        bool shiftFlag, bool rotateFlag)
{
    shiftRotateTransform(x, aux, shift, rotate, 1.0, shiftFlag, rotateFlag);

    return do_schaffer_F7_func({x.begin(), x.end()});
}

double schaffer_F7_func(const VectorRange& x)
{
    // justShift(x, 1.0) // not needed

    return do_schaffer_F7_func(x);
}

double step_rastrigin_func(std::vector<double>& x, std::vector<double>& aux,
                           const std::vector<double>& shift,
                           const std::vector<std::vector<double>>& rotate,
                           bool shiftFlag, bool rotateFlag)
{
    // ??? is exactly step rastrigin
    // TODO: maybe remove shift and rotate flag and always shift and rotate
    // see h01
    return rastrigin_func(x, aux, shift, rotate, shiftFlag, rotateFlag);
}

double levy_func(std::vector<double>& x, std::vector<double>& aux,
                 const std::vector<double>& shift,
                 const std::vector<std::vector<double>>& rotate, bool shiftFlag,
                 bool rotateFlag)
{
    shiftRotateTransform(x, aux, shift, rotate, 1.0, shiftFlag, rotateFlag);

    // TODO: Different than provided implementation. Check if min is 0,
    // otherwise change
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
                     bool shiftFlag, bool rotateFlag)
{
    shiftRotateTransform(x, aux, shift, rotate, 1.0, shiftFlag, rotateFlag);

    auto sum1 = 0.0;
    auto sum2 = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        sum1 += x[i] * x[i];
        sum2 += 0.5 * i * x[i];
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
        // TODO: Check implementation
        for (auto j = 1; j < 33; ++j) {
            const auto temp1 = std::pow(2.0, j);
            const auto temp2 = temp1 * (*std::next(x.begin, i));
            temp += std::fabs(temp2 - std::floor(temp2 + 0.5)) / temp1;
        }
        f *= std::pow(1.0 + (i + 1) * temp, 10.0 / temp3);
    }
    const auto temp1 = 10.0 / size / size;
    return f * temp1 - temp1;
}

double katsuura_func(std::vector<double>& x, std::vector<double>& aux,
                     const std::vector<double>& shift,
                     const std::vector<std::vector<double>>& rotate,
                     bool shiftFlag, bool rotateFlag)
{
    shiftRotateTransform(x, aux, shift, rotate, 5.0 / 100.0, shiftFlag,
                         rotateFlag);

    return do_katsuura_func({x.begin(), x.end()});
}

double katsuura_func(const VectorRange& x)
{
    justShift(x, 5.0 / 100.0);

    return do_katsuura_func(x);
}

double
hf01(std::vector<double>& x, std::vector<double>& aux,
     const std::vector<double>& shift,
     const std::vector<std::vector<double>>& rotate,
     const std::vector<std::size_t>& indices, bool shiftFlag, bool rotateFlag)
{
    // [0.4, 0.4, 0.2]
    shiftRotateTransform(x, aux, shift, rotate, 1.0, shiftFlag, rotateFlag);
    applyPermutation(x, aux, indices);

    const auto limit = std::ceil(0.4 * x.size());
    const auto margin_1 = std::next(x.begin(), limit); // 0.4
    const auto margin_2 = std::next(margin_1, limit);  // another 0.4

    const auto range1 = VectorRange{x.begin(), margin_1};
    const auto range2 = VectorRange{margin_1, margin_2};
    const auto range3 = VectorRange{margin_2, x.end()};

    // these methods do not shift and rotate
    return bent_cigar_func(range1) + hgbat_func(range2) +
           rastrigin_func(range3);
}

double
hf02(std::vector<double>& x, std::vector<double>& aux,
     const std::vector<double>& shift,
     const std::vector<std::vector<double>>& rotate,
     const std::vector<std::size_t>& indices, bool shiftFlag, bool rotateFlag)
{
    // [0.1, 0.2, 0.2, 0.2, 0.1, 0.2]
    // ??? omega = [10,20,30]. Check what's this and remove
    shiftRotateTransform(x, aux, shift, rotate, 1.0, shiftFlag, rotateFlag);
    applyPermutation(x, aux, indices);

    // TODO: n
    const auto limit1 = std::ceil(0.1 * x.size());
    const auto limit2 = std::ceil(0.2 * x.size());

    const auto margin1 = std::next(x.begin(), limit1); // 0.1
    const auto margin2 = std::next(margin1, limit2);   // 0.2
    const auto margin3 = std::next(margin2, limit2);   // 0.2
    const auto margin4 = std::next(margin3, limit2);   // 0.2
    const auto margin5 = std::next(margin4, limit1);   // 0.1
    // remaining: 0.2

    const auto range1 = VectorRange{x.begin(), margin1};
    const auto range2 = VectorRange{margin1, margin2};
    const auto range3 = VectorRange{margin2, margin3};
    const auto range4 = VectorRange{margin3, margin4};
    const auto range5 = VectorRange{margin4, margin5};
    const auto range6 = VectorRange{margin5, x.end()};

    return hgbat_func(range1) + katsuura_func(range2) + ackley_func(range3) +
           rastrigin_func(range4) + schwefel_func(range5) +
           schaffer_F7_func(range6);
}

double
hf03(std::vector<double>& x, std::vector<double>& aux,
     const std::vector<double>& shift,
     const std::vector<std::vector<double>>& rotate,
     const std::vector<std::size_t>& indices, bool shiftFlag, bool rotateFlag)
{
    //  [0.3, 0.2, 0.2, 0.1, 0.2]
    shiftRotateTransform(x, aux, shift, rotate, 1.0, shiftFlag, rotateFlag);
    applyPermutation(x, aux, indices);

    const auto limit1 = std::ceil(0.1 * x.size());
    const auto limit2 = std::ceil(0.2 * x.size());
    const auto limit3 = std::ceil(0.3 * x.size());

    const auto margin1 = std::next(x.begin(), limit3); // 0.3
    const auto margin2 = std::next(margin1, limit2);   // 0.2
    const auto margin3 = std::next(margin2, limit2);   // 0.2
    const auto margin4 = std::next(margin3, limit1);   // 0.1
    // remaining: 0.2

    const auto range1 = VectorRange{x.begin(), margin1};
    const auto range2 = VectorRange{margin1, margin2};
    const auto range3 = VectorRange{margin2, margin3};
    const auto range4 = VectorRange{margin3, margin4};
    const auto range5 = VectorRange{margin4, x.end()};

    return katsuura_func(range1) + happycat_func(range2) +
           grie_rosen_func(range3) + schwefel_func(range4) +
           ackley_func(range5);
}

double cf01(std::vector<double>& x, std::vector<double>& aux,
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
    // bias is 0, 200, 300, 100, 400
    const std::array<double, N> fit{
        rosenbrock_func(x, aux, shift.cbegin(), rotate.cbegin(), true,
                        rotateFlag),
        ellips_func(x, aux, shift_margin_1, rotate_margin_1, true, rotateFlag) *
                1e-6 +
            200,
        bent_cigar_func(x, aux, shift_margin_2, rotate_margin_2, true,
                        rotateFlag) *
                1e-6 +
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

double cf02(std::vector<double>& x, std::vector<double>& aux,
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

double cf03(std::vector<double>& x, std::vector<double>& aux,
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
    // bias = [0, 200, 300, 400, 200]
    const std::array<double, N> fit{
        escaffer6_func(x, aux, shift.cbegin(), rotate.cbegin(), true,
                       rotateFlag) *
            1e-26,
        schwefel_func(x, aux, shift_margin_1, rotate_margin_1, true,
                      rotateFlag) *
                10 +
            200,
        griewank_func(x, aux, shift_margin_2, rotate_margin_2, true,
                      rotateFlag) *
                1e-6 +
            300,
        rosenbrock_func(x, aux, shift_margin_3, rotate_margin_3, true,
                        rotateFlag) *
                10 +
            400,
        rastrigin_func(x, aux, shift_margin_4, rotate_margin_4, true,
                       rotateFlag) *
                5e-4 +
            200,

    };
    const std::array<int, N> delta{20, 20, 30, 30, 20};
    return compositionFunctionCalculator<N>(x, shift, delta, fit);
}

double cf04(std::vector<double>& x, std::vector<double>& aux,
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
    return 1;
}

} // namespace cec22
