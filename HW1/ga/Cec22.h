#include <vector>

namespace ga::functions {

// Ackley's
double ackley_func(std::vector<double>& x, const std::vector<double>& shift,
                   const std::vector<std::vector<double>>& rotate,
                   bool shift_flag, bool rotate_flag);

// Bent Cigar
double bent_cigar_func(std::vector<double>& x, const std::vector<double>& shift,
                       const std::vector<std::vector<double>>& rotate,
                       bool shift_flag, bool rotate_flag);

// Discus
double discus_func(std::vector<double>& x, const std::vector<double>& shift,
                   const std::vector<std::vector<double>>& rotate,
                   bool shift_flag, bool rotate_flag);

// High Conditioned Elliptic Function
double ellips_func(std::vector<double>& x, const std::vector<double>& shift,
                   const std::vector<std::vector<double>>& rotate,
                   bool shift_flag, bool rotate_flag);

// Expanded Schaffer’s f6
double escaffer6_func(std::vector<double>& x, const std::vector<double>& shift,
                      const std::vector<std::vector<double>>& rotate,
                      bool shift_flag, bool rotate_flag);

// Griewank's
double griewank_func(std::vector<double>& x, const std::vector<double>& shift,
                     const std::vector<std::vector<double>>& rotate,
                     bool shift_flag, bool rotate_flag);

// Griewank-Rosenbrock
double grie_rosen_func(std::vector<double>& x, const std::vector<double>& shift,
                       const std::vector<std::vector<double>>& rotate,
                       bool shift_flag, bool rotate_flag);

// HappyCat
double happycat_func(std::vector<double>& x, const std::vector<double>& shift,
                     const std::vector<std::vector<double>>& rotate,
                     bool shift_flag, bool rotate_flag);

// HGBat
double hgbat_func(std::vector<double>& x, const std::vector<double>& shift,
                  const std::vector<std::vector<double>>& rotate,
                  bool shift_flag, bool rotate_flag);

// Rosenbrock's
double rosenbrock_func(std::vector<double>& x, const std::vector<double>& shift,
                       const std::vector<std::vector<double>>& rotate,
                       bool shift_flag, bool rotate_flag);

// Rastrigin's
double rastrigin_func(std::vector<double>& x, const std::vector<double>& shift,
                      const std::vector<std::vector<double>>& rotate,
                      bool shift_flag, bool rotate_flag);

// Schwefel's
double schwefel_func(std::vector<double>& x, const std::vector<double>& shift,
                     const std::vector<std::vector<double>>& rotate,
                     bool shift_flag, bool rotate_flag);

// Schaffer’s F7
double
schaffer_F7_func(std::vector<double>& x, const std::vector<double>& shift,
                 const std::vector<std::vector<double>>& rotate,
                 bool shift_flag, bool rotate_flag);

// Noncontinuous Rastrigin's
double
step_rastrigin_func(std::vector<double>& x, const std::vector<double>& shift,
                    const std::vector<std::vector<double>>& rotate,
                    bool shift_flag, bool rotate_flag);

// Levy
double levy_func(std::vector<double>& x, const std::vector<double>& shift,
                 const std::vector<std::vector<double>>& rotate,
                 bool shift_flag, bool rotate_flag);

// ZAKHAROV
double zakharov_func(std::vector<double>& x, const std::vector<double>& shift,
                     const std::vector<std::vector<double>>& rotate,
                     bool shift_flag, bool rotate_flag);

// Katsuura
double katsuura_func(std::vector<double>& x, const std::vector<double>& shift,
                     const std::vector<std::vector<double>>& rotate,
                     bool shift_flag, bool rotate_flag);

// Hybrid Function 1
double hf01(std::vector<double>& x, const std::vector<double>& shift,
            const std::vector<std::vector<double>>& rotate,
            const std::vector<std::size_t>& indices, bool shift_flag,
            bool rotate_flag);

// Hybrid Function 2
double hf02(std::vector<double>& x, const std::vector<double>& shift,
            const std::vector<std::vector<double>>& rotate,
            const std::vector<std::size_t>& indices, bool shift_flag,
            bool rotate_flag);

// Hybrid Function 3
double hf03(std::vector<double>& x, const std::vector<double>& shift,
            const std::vector<std::vector<double>>& rotate,
            const std::vector<std::size_t>& indices, bool shift_flag,
            bool rotate_flag);

// Composition Function 1
double cf01(std::vector<double>& x, const std::vector<double>& shift,
            const std::vector<std::vector<double>>& rotate, bool rotate_flag);

// Composition Function 1
double cf02(std::vector<double>& x, const std::vector<double>& shift,
            const std::vector<std::vector<double>>& rotate, bool rotate_flag);

// Composition Function 1
double cf03(std::vector<double>& x, const std::vector<double>& shift,
            const std::vector<std::vector<double>>& rotate, bool rotate_flag);

// Composition Function 1
double cf04(std::vector<double>& x, const std::vector<double>& shift,
            const std::vector<std::vector<double>>& rotate, bool rotate_flag);

int sanity_check();

} // namespace ga::functions