#include <random>
#include <string>

namespace ga {

class GeneticAlgorithm
{
  public:
    // TODO: add a function manager to call the functions
    GeneticAlgorithm() = default;
    GeneticAlgorithm(double crossoverProbability, double mutationProbability,
                     double hypermutationRate, double elitesPercentage,
                     double selectionPressure, double encodingChangeRate,
                     int iterations, int populationSize,
                     int stepsToHypermutation, int maxNoImprovementSteps);
    void sanityCheck();

  private: // TODO: find good values
    double crossoverProbability = 0.3;
    double hypermutationRate = 0.1;
    double elitesPercentage = 0.04;
    double selectionPressure = 1.01;
    double encodingChangeRate = 0.1;
    double mutationProbability = 0.005;
    int iterations = 15; // TODO: might not need iterations here
    int populationSize = 100;
    int stepsToHypermutation = 10;
    int maxNoImprovementSteps = 100;
    std::random_device seed;
    std::mt19937_64 gen{seed()};
    std::bernoulli_distribution randomBool;
    std::uniform_real_distribution<double> randomDouble{0.0, 1.0};
};

} // namespace ga
