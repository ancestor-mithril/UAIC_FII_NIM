Compiler: gcc-11.2.0
On linux please change the fifth line in the Makefile from g++ to g++-11

Use `make rel1` to compile all
Use `make exp` to compile and run experiments

Main arguments:
 * no argument: Runs the test function
 * FunctionName: Runs 30 runs using the function provided it exists (see make exp for examples)

Don't use `make rel2` or `make debug` on linux because it uses CMake with MinGW Makefiles


