#include "cec22/Cec22.h"

#include <iostream>

int main(int argc, char* argv[])
{
    std::cout << cec22::sanity_check() << '\n';
    return 0;
}