cmake_minimum_required(VERSION 3.23)
project(main)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(PREBUILT_MODULE_PATH ${CMAKE_BINARY_DIR}/modules)

function(add_module name)
    file(MAKE_DIRECTORY ${PREBUILT_MODULE_PATH})
    add_custom_target(${name}.pcm
            COMMAND
                ${CMAKE_CXX_COMPILER}
                -std=c++20
                -Mmodules
                -fmodules-ts
                -c
                ${CMAKE_CURRENT_SOURCE_DIR}/${ARGN}
                -emit-module-interface
                -o ${PREBUILT_MODULE_PATH}/${name}.pcm

            )
endfunction()


add_compile_options(-Mmodules)
add_compile_options(-fmodules-ts)
add_compile_options(-fprebuilt-module-path=${PREBUILT_MODULE_PATH})

add_module(helloworld helloworld.cpp)

add_executable(main main.cpp)
add_dependencies(main helloworld.pcm)