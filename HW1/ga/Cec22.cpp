module;
#include <vector>
export module GA.Cec22Specific;
// clang-format off
// import <vector>;
// clang-format on

namespace ga::functions {

// Ackley's
export double
ackley_func(std::vector<double>& x, const std::vector<double>& shift,
            const std::vector<std::vector<double>>& rotate, bool shift_flag,
            bool rotate_flag);

// Bent Cigar
export double
bent_cigar_func(std::vector<double>& x, const std::vector<double>& shift,
                const std::vector<std::vector<double>>& rotate, bool shift_flag,
                bool rotate_flag);

// Discus
export double
discus_func(std::vector<double>& x, const std::vector<double>& shift,
            const std::vector<std::vector<double>>& rotate, bool shift_flag,
            bool rotate_flag);

// High Conditioned Elliptic Function
export double
ellips_func(std::vector<double>& x, const std::vector<double>& shift,
            const std::vector<std::vector<double>>& rotate, bool shift_flag,
            bool rotate_flag);

// Expanded Schaffer’s f6
export double
escaffer6_func(std::vector<double>& x, const std::vector<double>& shift,
               const std::vector<std::vector<double>>& rotate, bool shift_flag,
               bool rotate_flag);

// Griewank's
export double
griewank_func(std::vector<double>& x, const std::vector<double>& shift,
              const std::vector<std::vector<double>>& rotate, bool shift_flag,
              bool rotate_flag);

// Griewank-Rosenbrock
export double
grie_rosen_func(std::vector<double>& x, const std::vector<double>& shift,
                const std::vector<std::vector<double>>& rotate, bool shift_flag,
                bool rotate_flag);

// HappyCat
export double
happycat_func(std::vector<double>& x, const std::vector<double>& shift,
              const std::vector<std::vector<double>>& rotate, bool shift_flag,
              bool rotate_flag);

// HGBat
export double
hgbat_func(std::vector<double>& x, const std::vector<double>& shift,
           const std::vector<std::vector<double>>& rotate, bool shift_flag,
           bool rotate_flag);

// Rosenbrock's
export double
rosenbrock_func(std::vector<double>& x, const std::vector<double>& shift,
                const std::vector<std::vector<double>>& rotate, bool shift_flag,
                bool rotate_flag);

// Rastrigin's
export double
rastrigin_func(std::vector<double>& x, const std::vector<double>& shift,
               const std::vector<std::vector<double>>& rotate, bool shift_flag,
               bool rotate_flag);

// Schwefel's
export double
schwefel_func(std::vector<double>& x, const std::vector<double>& shift,
              const std::vector<std::vector<double>>& rotate, bool shift_flag,
              bool rotate_flag);

// Schaffer’s F7
export double
schaffer_F7_func(std::vector<double>& x, const std::vector<double>& shift,
                 const std::vector<std::vector<double>>& rotate,
                 bool shift_flag, bool rotate_flag);

// Noncontinuous Rastrigin's
export double
step_rastrigin_func(std::vector<double>& x, const std::vector<double>& shift,
                    const std::vector<std::vector<double>>& rotate,
                    bool shift_flag, bool rotate_flag);

// Levy
export double
levy_func(std::vector<double>& x, const std::vector<double>& shift,
          const std::vector<std::vector<double>>& rotate, bool shift_flag,
          bool rotate_flag);

// ZAKHAROV
export double
zakharov_func(std::vector<double>& x, const std::vector<double>& shift,
              const std::vector<std::vector<double>>& rotate, bool shift_flag,
              bool rotate_flag);

// Katsuura
export double
katsuura_func(std::vector<double>& x, const std::vector<double>& shift,
              const std::vector<std::vector<double>>& rotate, bool shift_flag,
              bool rotate_flag);

// Hybrid Function 1
export double hf01(std::vector<double>& x, const std::vector<double>& shift,
                   const std::vector<std::vector<double>>& rotate,
                   const std::vector<std::size_t>& indices, bool shift_flag,
                   bool rotate_flag);

// Hybrid Function 2
export double hf02(std::vector<double>& x, const std::vector<double>& shift,
                   const std::vector<std::vector<double>>& rotate,
                   const std::vector<std::size_t>& indices, bool shift_flag,
                   bool rotate_flag);

// Hybrid Function 3
export double hf03(std::vector<double>& x, const std::vector<double>& shift,
                   const std::vector<std::vector<double>>& rotate,
                   const std::vector<std::size_t>& indices, bool shift_flag,
                   bool rotate_flag);

// Composition Function 1
export double cf01(std::vector<double>& x, const std::vector<double>& shift,
                   const std::vector<std::vector<double>>& rotate,
                   bool shift_flag, bool rotate_flag);

// Composition Function 1
export double cf02(std::vector<double>& x, const std::vector<double>& shift,
                   const std::vector<std::vector<double>>& rotate,
                   bool shift_flag, bool rotate_flag);

// Composition Function 1
export double cf03(std::vector<double>& x, const std::vector<double>& shift,
                   const std::vector<std::vector<double>>& rotate,
                   bool shift_flag, bool rotate_flag);

// Composition Function 1
export double cf04(std::vector<double>& x, const std::vector<double>& shift,
                   const std::vector<std::vector<double>>& rotate,
                   bool shift_flag, bool rotate_flag);

export int sanity_check();

} // namespace ga::functions