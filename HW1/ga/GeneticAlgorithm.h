#include <random>
#include <string>

namespace ga {

// TODO: Maybe use template for population size
class GeneticAlgorithm
{
  public:
    // TODO: add a function manager to call the functions
    GeneticAlgorithm(double crossoverProbability, double mutationProbability,
                     double hypermutationRate, double elitesPercentage,
                     double selectionPressure, double encodingChangeRate,
                     int maxSteps, int populationSize, int dimensions,
                     int stepsToHypermutation, int maxNoImprovementSteps);
    void sanityCheck() const;
    void run();

  private:
    void randomizePopulation();

    // const?
    std::random_device seed;
    std::mt19937_64 gen{seed()};
    std::bernoulli_distribution randomBool;
    std::uniform_real_distribution<double> randomDouble{0.0, 1.0};

    // TODO: find good values
    // TODO: use const where we should use const
    double crossoverProbability = 0.3;  // can this modify?
    double mutationProbability = 0.005; // can this modify?
    double hypermutationRate = 0.1;     // can this modify?
    double elitesPercentage = 0.04;     // can this modify?
    double selectionPressure = 1.01;    // can this modify?
    double encodingChangeRate = 0.1;    // can this modify?

    const int maxSteps = 200'000;
    const int populationSize = 100; // should we always use 100?
    const int dimensions = 10;
    const int bitsPerChromozome;
    const int stepsToHypermutation = 10;
    const int maxNoImprovementSteps = 100;

    int elitesNumber = 4; // can this modify?

    // TODO: test against char, might be faster because std::vector<double> is a
    // space efficient but not time efficient
    std::vector<std::vector<bool>> population;
};

GeneticAlgorithm getDefault();

} // namespace ga
