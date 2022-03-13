cmake_minimum_required(VERSION 3.23)

project(
  "GA"
  VERSION 1.0.0
  LANGUAGES CXX
)

set(OPT O3)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_compile_options(-std=c++20 
                    -${OPT}
                    # -fmodules-ts
                    -Wall 
                    -Wextra 
                    -Wpedantic)

add_custom_target(std_modules ALL
    COMMAND ${CMAKE_COMMAND} -E echo "Building standard library modules"
    # COMMAND g++ -fmodules-ts -std=c++20 -c -xc++-system-header string
    # COMMAND g++ -fmodules-ts -std=c++20 -c -xc++-system-header iostream
    # COMMAND g++ -fmodules-ts -std=c++20 -${OPT} -c -xc++-system-header vector
    # COMMAND g++ -fmodules-ts -std=c++20 -${OPT} -c -xc++-system-header algorithm
    # COMMAND g++ -fmodules-ts -std=c++20 -${OPT} -c -xc++-system-header functional
    # COMMAND g++ -fmodules-ts -std=c++20 -${OPT} -c -xc++-system-header cmath
    # COMMAND g++ -fmodules-ts -std=c++20 -${OPT} -c -xc++-system-header numeric
    # COMMAND g++ -fmodules-ts -std=c++20 -${OPT} -c -xc++-system-header unordered_map
    # COMMAND g++ -fmodules-ts -std=c++20 -${OPT} -c -xc++-system-header tuple
    # COMMAND g++ -fmodules-ts -std=c++20 -${OPT} -c -xc++-system-header random
)

function (prepare_for_module TGT)
    set_property(TARGET ${TGT} PROPERTY CXX_STANDARD 20)
    set_property(TARGET ${TGT} PROPERTY CXX_EXTENSIONS OFF)
    add_dependencies(${TGT} std_modules)
endfunction()


set (TARGET program)
set (SOURCES ga/main.cpp)
set (MODULES 
    ga/Cec22.cpp
    # ga/Cec22Impl.cpp
    ga/GeneticAlgorithm.cpp
    ga/FunctionManager.cpp
    # ga/GeneticAlgorithmImpl.cpp
    )

# Setup program modules object library
set (MODULE_TARGET module-lib)
add_library(${MODULE_TARGET} OBJECT ${MODULES})
# prepare_for_module(${MODULE_TARGET})

# Setup executable
add_executable(${TARGET} ${SOURCES})
# prepare_for_module(${TARGET})

# Add modules to application using object library
target_link_libraries(${TARGET} PRIVATE ${MODULE_TARGET})
# Test end
