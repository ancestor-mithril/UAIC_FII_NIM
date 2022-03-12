#include "GeneticAlgorithm.h"

#include <iostream>

namespace ga {

// GeneticAlgorithm::GeneticAlgorithm(
//     double crossoverProbability, double mutationProbability,
//     double hypermutationRate, double elitesPercentage, double
//     selectionPressure, double encodingChangeRate, int iterations, int
//     populationSize, int stepsToHypermutation, int maxNoImprovementSteps) :
//     crossoverProbability{crossoverProbability},
//       mutationProbability{mutationProbability},
//       hypermutationRate{hypermutationRate},
//       elitesPercentage{elitesPercentage},
//       selectionPressure{selectionPressure},
//       encodingChangeRate{encodingChangeRate}, iterations{iterations},
//       stepsToHypermutation{stepsToHypermutation}, maxNoImprovementSteps{
//                                                       maxNoImprovementSteps}
// {
// }

void GeneticAlgorithm::sanityCheck()
{
    std::cout << "GeneticAlgorithm::sanityCheck" << '\n';
}

} // namespace ga
