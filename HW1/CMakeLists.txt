cmake_minimum_required(VERSION 3.23)

project(
  "GA"
  VERSION 1.0.0
  LANGUAGES CXX
)

set(OPT O3)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20 -${OPT} -fmodules-ts")
add_compile_options(-std=c++20 -${OPT} -fmodules-ts)

add_custom_target(std_modules ALL
    COMMAND ${CMAKE_COMMAND} -E echo "Building standard library modules"
    COMMAND g++ -fmodules-ts -std=c++20 -c -xc++-system-header string
    COMMAND g++ -fmodules-ts -std=c++20 -c -xc++-system-header iostream
    COMMAND g++ -fmodules-ts -std=c++20 -${OPT} -c -xc++-system-header vector
    COMMAND g++ -fmodules-ts -std=c++20 -${OPT} -c -xc++-system-header algorithm
    COMMAND g++ -fmodules-ts -std=c++20 -${OPT} -c -xc++-system-header functional
    COMMAND g++ -fmodules-ts -std=c++20 -${OPT} -c -xc++-system-header cmath
    COMMAND g++ -fmodules-ts -std=c++20 -${OPT} -c -xc++-system-header numeric
)

function (prepare_for_module TGT)
    target_compile_options(${TGT} PUBLIC -fmodules-ts)
    set_property(TARGET ${TGT} PROPERTY CXX_STANDARD 20)
    set_property(TARGET ${TGT} PROPERTY CXX_EXTENSIONS OFF)
    add_dependencies(${TGT} std_modules)
endfunction()

# Testing C++20:
set (TARGET test)
set (SOURCES ga/main.cpp)
set (MODULES 
    ga/Cec22.cpp
    ga/Cec22Impl.cpp)

# Setup program modules object library
set (MODULE_TARGET module-lib)
add_library(${MODULE_TARGET} OBJECT ${MODULES})
prepare_for_module(${MODULE_TARGET})

# Setup executable
add_executable(${TARGET} ${SOURCES})
target_link_libraries(${TARGET} ${Boost_LIBRARIES})
prepare_for_module(${TARGET})

# Add modules to application using object library
target_link_libraries(${TARGET} PRIVATE ${MODULE_TARGET})
# Test end
