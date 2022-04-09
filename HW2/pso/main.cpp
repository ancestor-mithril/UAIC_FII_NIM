#include "cec22/Cec22.h"

#include <iostream>

int main([[maybe_unused]] int argc, [[maybe_unused]] char* argv[])
{
    std::cout << cec22::sanity_check() << '\n';
    return 0;
}
