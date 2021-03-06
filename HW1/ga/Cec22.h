#pragma once
#include <vector>

namespace cec22 {
// Minimum for these functions is 0 in unshifted and not rotated [0]n
// ?? except for composition functions

// Ackley's
double ackley_func(const std::vector<double>& x, std::vector<double>& aux,
                   const std::vector<double>& shift,
                   const std::vector<std::vector<double>>& rotate,
                   bool shift_flag, bool rotate_flag);

// Bent Cigar
double bent_cigar_func(const std::vector<double>& x, std::vector<double>& aux,
                       const std::vector<double>& shift,
                       const std::vector<std::vector<double>>& rotate,
                       bool shift_flag, bool rotate_flag);

// Discus
double discus_func(const std::vector<double>& x, std::vector<double>& aux,
                   const std::vector<double>& shift,
                   const std::vector<std::vector<double>>& rotate,
                   bool shift_flag, bool rotate_flag);

// High Conditioned Elliptic Function
double ellips_func(const std::vector<double>& x, std::vector<double>& aux,
                   const std::vector<double>& shift,
                   const std::vector<std::vector<double>>& rotate,
                   bool shift_flag, bool rotate_flag);

// Expanded Schaffer’s f6
double escaffer6_func(const std::vector<double>& x, std::vector<double>& aux,
                      const std::vector<double>& shift,
                      const std::vector<std::vector<double>>& rotate,
                      bool shift_flag, bool rotate_flag);

// Griewank's
double griewank_func(const std::vector<double>& x, std::vector<double>& aux,
                     const std::vector<double>& shift,
                     const std::vector<std::vector<double>>& rotate,
                     bool shift_flag, bool rotate_flag);

// Griewank-Rosenbrock
double grie_rosen_func(const std::vector<double>& x, std::vector<double>& aux,
                       const std::vector<double>& shift,
                       const std::vector<std::vector<double>>& rotate,
                       bool shift_flag, bool rotate_flag);

// HappyCat
double happycat_func(const std::vector<double>& x, std::vector<double>& aux,
                     const std::vector<double>& shift,
                     const std::vector<std::vector<double>>& rotate,
                     bool shift_flag, bool rotate_flag);

// HGBat
double hgbat_func(const std::vector<double>& x, std::vector<double>& aux,
                  const std::vector<double>& shift,
                  const std::vector<std::vector<double>>& rotate,
                  bool shift_flag, bool rotate_flag);

// Rosenbrock's
double rosenbrock_func(const std::vector<double>& x, std::vector<double>& aux,
                       const std::vector<double>& shift,
                       const std::vector<std::vector<double>>& rotate,
                       bool shift_flag, bool rotate_flag);

// Rastrigin's
double rastrigin_func(const std::vector<double>& x, std::vector<double>& aux,
                      const std::vector<double>& shift,
                      const std::vector<std::vector<double>>& rotate,
                      bool shift_flag, bool rotate_flag);

// Schwefel's
double schwefel_func(const std::vector<double>& x, std::vector<double>& aux,
                     const std::vector<double>& shift,
                     const std::vector<std::vector<double>>& rotate,
                     bool shift_flag, bool rotate_flag);

// Schaffer’s F7
double schaffer_F7_func(const std::vector<double>& x, std::vector<double>& aux,
                        const std::vector<double>& shift,
                        const std::vector<std::vector<double>>& rotate,
                        bool shift_flag, bool rotate_flag);

// Noncontinuous Rastrigin's
double
step_rastrigin_func(const std::vector<double>& x, std::vector<double>& aux,
                    const std::vector<double>& shift,
                    const std::vector<std::vector<double>>& rotate,
                    bool shift_flag, bool rotate_flag);

// Levy
double levy_func(const std::vector<double>& x, std::vector<double>& aux,
                 const std::vector<double>& shift,
                 const std::vector<std::vector<double>>& rotate,
                 bool shift_flag, bool rotate_flag);

// ZAKHAROV
double zakharov_func(const std::vector<double>& x, std::vector<double>& aux,
                     const std::vector<double>& shift,
                     const std::vector<std::vector<double>>& rotate,
                     bool shift_flag, bool rotate_flag);

// Katsuura
double katsuura_func(const std::vector<double>& x, std::vector<double>& aux,
                     const std::vector<double>& shift,
                     const std::vector<std::vector<double>>& rotate,
                     bool shift_flag, bool rotate_flag);

// Hybrid Function 1
double hf01(const std::vector<double>& x, std::vector<double>& aux,
            const std::vector<double>& shift,
            const std::vector<std::vector<double>>& rotate,
            const std::vector<std::size_t>& indices, bool shift_flag,
            bool rotate_flag);

// Hybrid Function 2
double hf02(const std::vector<double>& x, std::vector<double>& aux,
            const std::vector<double>& shift,
            const std::vector<std::vector<double>>& rotate,
            const std::vector<std::size_t>& indices, bool shift_flag,
            bool rotate_flag);

// Hybrid Function 3
double hf03(const std::vector<double>& x, std::vector<double>& aux,
            const std::vector<double>& shift,
            const std::vector<std::vector<double>>& rotate,
            const std::vector<std::size_t>& indices, bool shift_flag,
            bool rotate_flag);

// Composition Function 1
double cf01(const std::vector<double>& x, std::vector<double>& aux,
            const std::vector<double>& shift,
            const std::vector<std::vector<double>>& rotate, bool rotate_flag);

// Composition Function 1
double cf02(const std::vector<double>& x, std::vector<double>& aux,
            const std::vector<double>& shift,
            const std::vector<std::vector<double>>& rotate, bool rotate_flag);

// Composition Function 1
double cf03(const std::vector<double>& x, std::vector<double>& aux,
            const std::vector<double>& shift,
            const std::vector<std::vector<double>>& rotate, bool rotate_flag);

// Composition Function 1
double cf04(const std::vector<double>& x, std::vector<double>& aux,
            const std::vector<double>& shift,
            const std::vector<std::vector<double>>& rotate, bool rotate_flag);

int sanity_check();

} // namespace cec22