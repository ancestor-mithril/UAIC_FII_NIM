#pragma once

#include <map>
#include <vector>

namespace function_layer::cache_layer {

struct CacheLayer {
  public:
    // map because we need lower_bound
    std::map<std::vector<double>, double> cache;
    // this will be slow

    // TODO: implement caching strategy (keep-recent) to reduce size
    // TODO: Add Custom comparator which maintains uniqueness and takes distance
    // into account.
    // Solution feasible in a domain in which a1 >= a2 >= .. >= an: L2 norm
    // L2 norm of such vectors is both unique (i hope) and offers sorting by
    // closeness
};
} // namespace function_layer::cache_layer
