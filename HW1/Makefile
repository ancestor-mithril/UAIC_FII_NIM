BUILDDIR=build
APP=main

cxx:
	export CXX=g++

builddir:
	if [[ -d "${BUILDDIR}" ]]; then echo "Directory ${BUILDDIR} already exists"; else echo "Creating ${BUILDDIR}"; mkdir ${BUILDDIR}; fi;

clean:
	rm -rf ${BUILDDIR}

release: builddir cxx 
	# cmake \ # TODO: Repair CmakeListss
	# 	-D CMAKE_BUILD_TYPE=Release \
	# 	-D CMAKE_CXX_COMPILER=g++ \
	# 	-G "MinGW Makefiles"  \
	# 	-S . \
	# 	-B ${BUILDDIR}
	# cmake --build ${BUILDDIR}
	cd ${BUILDDIR} \
	&& g++ -std=c++20 -fmodules-ts -xc++-system-header string \
	&& g++ -std=c++20 -fmodules-ts -xc++-system-header iostream \
	&& g++ -std=c++20 -fmodules-ts -c ../src/helloworld.cpp \
	&& g++ -std=c++20 -fmodules-ts -c ../src/main.cpp \
	&& g++ -std=c++20 helloworld.o main.o -o main 

debug: builddir cxx
	cmake \
		-D CMAKE_BUILD_TYPE=Debug \
		-D CMAKE_CXX_COMPILER=g++ \
		-G "MinGW Makefiles"  \
		-S . \
		-B ${BUILDDIR}
	cmake --build ${BUILDDIR}

run: release
	./${BUILDDIR}/${APP}.exe 

gdb: debug
	gdb ./${BUILDDIR}/${APP}.exe

all: clean release