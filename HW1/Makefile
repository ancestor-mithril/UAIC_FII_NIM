BUILDDIR=build
APP=test
OPT=O3
CMAKE_CXX_FLAGS=-std=c++20 -${OPT} -fmodules-ts -Wall -Wextra -Wpedantic

cxx:
	export CXX=g++

format:
	clang-format -i ga/Cec22Impl.cpp
	clang-format -i ga/Cec22.cpp
	# clang-format -i ga/main.cpp

builddir:
	if [[ -d "${BUILDDIR}" ]]; then echo "Directory ${BUILDDIR} already exists"; else echo "Creating ${BUILDDIR}"; mkdir ${BUILDDIR}; fi;

clean:
	rm -rf ${BUILDDIR}

rel1: clean builddir cxx # We need clean because building std modules comes with gcm.cache errors
	cd ${BUILDDIR} \
	&& g++ ${CMAKE_CXX_FLAGS} -c -xc++-system-header iostream \
	&& g++ ${CMAKE_CXX_FLAGS} -c -xc++-system-header vector \
	&& g++ ${CMAKE_CXX_FLAGS} -c -xc++-system-header algorithm \
	&& g++ ${CMAKE_CXX_FLAGS} -c -xc++-system-header functional \
	&& g++ ${CMAKE_CXX_FLAGS} -c -xc++-system-header cmath \
	&& g++ ${CMAKE_CXX_FLAGS} -c -xc++-system-header numeric \
	&& g++ ${CMAKE_CXX_FLAGS} -c -xc++-system-header unordered_map \
	&& g++ ${CMAKE_CXX_FLAGS} -c -xc++-system-header tuple \
	&& g++ ${CMAKE_CXX_FLAGS} -c ../ga/Cec22.cpp ../ga/Cec22Impl.cpp \
	&& g++ ${CMAKE_CXX_FLAGS} -c ../ga/main.cpp \
	&& g++ ${CMAKE_CXX_FLAGS} Cec22.o Cec22Impl.o main.o -o test.exe

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
