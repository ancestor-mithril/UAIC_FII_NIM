BUILDDIR=build
APP=program
OPT=Ofast # Check O2 against Ofast.
CMAKE_CXX_FLAGS=-std=c++20 -${OPT} -Wall -Wextra -Wpedantic

cxx:
	export CXX=g++

format:
	clang-format -i ga/Constexpr.h
	clang-format -i ga/Constants.h
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

rel1: builddir cxx
	cd ${BUILDDIR} \
	&& g++ ${CMAKE_CXX_FLAGS} -c ../ga/Cec22.cpp \
	&& g++ ${CMAKE_CXX_FLAGS} -c ../ga/GeneticAlgorithm.cpp \
	&& g++ ${CMAKE_CXX_FLAGS} -c ../ga/FunctionManager.cpp \
	&& g++ ${CMAKE_CXX_FLAGS} -c ../ga/main.cpp \
	&& g++ ${CMAKE_CXX_FLAGS} Cec22.o FunctionManager.o GeneticAlgorithm.o main.o -o ${APP}.exe

main: builddir cxx  # debug only
	cd ${BUILDDIR} \
	&& g++ ${CMAKE_CXX_FLAGS} -c ../ga/main.cpp \
	&& g++ ${CMAKE_CXX_FLAGS}  main.o -o ${APP}.exe
	./${BUILDDIR}/${APP}.exe 

rel2: builddir cxx
	cmake \
		-D CMAKE_BUILD_TYPE=Release \
		-D CMAKE_CXX_COMPILER=g++ \
		-G "MinGW Makefiles"  \
		-S . \
		-B ${BUILDDIR}
	cmake --build ${BUILDDIR}

debug: builddir cxx
	cmake \
		-D CMAKE_BUILD_TYPE=Debug \
		-D CMAKE_CXX_COMPILER=g++ \
		-G "MinGW Makefiles"  \
		-S . \
		-B ${BUILDDIR}
	cmake --build ${BUILDDIR}

run:
	./${BUILDDIR}/${APP}.exe 

r1: rel1 run

r2: rel2 run

gdb: debug
	gdb ./${BUILDDIR}/${APP}.exe

all: clean format run
