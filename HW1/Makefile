BUILDDIR=build
APP=program
OPT=O3
CMAKE_CXX_FLAGS=-std=c++20 -${OPT} -Wall -Wextra -Wpedantic

cxx:
	export CXX=g++

format:
	clang-format -i ga/Cec22.h
	clang-format -i ga/Cec22.cpp
	clang-format -i ga/GeneticAlgorithm.h
	clang-format -i ga/GeneticAlgorithm.cpp
	clang-format -i ga/FunctionManager.h
	clang-format -i ga/FunctionManager.cpp
	clang-format -i ga/main.cpp

builddir:
	if [[ -d "${BUILDDIR}" ]]; then echo "Directory ${BUILDDIR} already exists"; else echo "Creating ${BUILDDIR}"; mkdir ${BUILDDIR}; fi;

clean:
	rm -rf ${BUILDDIR}

rel1: clean builddir cxx # We need clean because building std modules comes with gcm.cache errors
	cd ${BUILDDIR} \
	&& g++ ${CMAKE_CXX_FLAGS} -c ../ga/Cec22.cpp \
	&& g++ ${CMAKE_CXX_FLAGS} -c ../ga/GeneticAlgorithm.cpp \
	&& g++ ${CMAKE_CXX_FLAGS} -c ../ga/FunctionManager.cpp \
	&& g++ ${CMAKE_CXX_FLAGS} -c ../ga/main.cpp \
	&& g++ ${CMAKE_CXX_FLAGS} Cec22.o FunctionManager.o GeneticAlgorithm.o main.o -o ${APP}.exe

rel2: clean builddir cxx
	cmake \
		-D CMAKE_BUILD_TYPE=Release \
		-D CMAKE_CXX_COMPILER=g++ \
		-G "MinGW Makefiles"  \
		-S . \
		-B ${BUILDDIR}
	cmake --build ${BUILDDIR}

debug: clean builddir cxx # We need clean because building std modules comes with gcm.cache errors
	cmake \
		-D CMAKE_BUILD_TYPE=Debug \
		-D CMAKE_CXX_COMPILER=g++ \
		-G "MinGW Makefiles"  \
		-S . \
		-B ${BUILDDIR}
	cmake --build ${BUILDDIR}

r1: rel1
	./${BUILDDIR}/${APP}.exe 

r2: rel2
	./${BUILDDIR}/${APP}.exe 

gdb: debug
	gdb ./${BUILDDIR}/${APP}.exe

all: clean format run
