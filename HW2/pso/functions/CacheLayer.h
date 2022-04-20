#pragma once

#include <map>
#include <vector>

namespace function_layer::cache_layer {

// TODO: use templates and concepts
// TODO: use static polymorphism

struct VectorCache {
    std::vector<std::vector<double>> cache;
    std::vector<double> values;
};

struct MapLayer {
    // map because we need lower_bound
    std::map<std::vector<double>, double> cache;
    // this will be slow. Key is huge, but value is small. It would be more
    // efficient to use 2 sorted std::vector.

    // TODO: implement caching strategy (keep-recent) to reduce size
    // TODO: Add Custom comparator which maintains uniqueness and takes distance
    // into account.
    // Solution feasible in a domain in which a1 >= a2 >= .. >= an: L2 norm
    // L2 norm of such vectors is both unique and offers sorting by
    // closeness

    // TODO: Another strategy is to use cosine similarity
};
} // namespace function_layer::cache_layer
