BUILDDIR=build
APP=test
OPT=O3

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
	&& g++ -std=c++20 -fmodules-ts -c -xc++-system-header iostream \
	&& g++ -std=c++20 -${OPT} -fmodules-ts -c -xc++-system-header vector \
	&& g++ -std=c++20 -${OPT} -fmodules-ts -c -xc++-system-header algorithm \
	&& g++ -std=c++20 -${OPT} -fmodules-ts -c -xc++-system-header functional \
	&& g++ -std=c++20 -${OPT} -fmodules-ts -c -xc++-system-header cmath \
	&& g++ -std=c++20 -${OPT} -fmodules-ts -c -xc++-system-header numeric \
	&& g++ -std=c++20 -${OPT} -fmodules-ts -c ../ga/Cec22.cpp ../ga/Cec22Impl.cpp \
	&& g++ -std=c++20 -${OPT} -fmodules-ts -c ../ga/main.cpp \
	&& g++ -std=c++20 -${OPT} Cec22.o Cec22Impl.o main.o -o test.exe

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
