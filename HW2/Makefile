BUILDDIR=build
APP=app
OPT=Ofast
CMAKE_CXX_FLAGS=-std=c++20 -${OPT} -Wall -Wextra -Wpedantic
GCC=g++

builddir:
	if [[ -d "${BUILDDIR}" ]]; then echo "Directory ${BUILDDIR} already exists"; else echo "Creating ${BUILDDIR}"; mkdir ${BUILDDIR}; fi;

clean:
	rm -rf ${BUILDDIR}

format:
	clang-format -i ./pso/cec22/Cec22.cpp
	clang-format -i ./pso/cec22/Cec22.h
	clang-format -i ./pso/functions/FunctionManager.h
	clang-format -i ./pso/functions/FunctionManager.cpp
	clang-format -i ./pso/main.cpp

release: builddir
	cd ${BUILDDIR} \
	&& ${GCC} ${CMAKE_CXX_FLAGS} -c ../pso/cec22/Cec22.cpp \
	&& ${GCC} ${CMAKE_CXX_FLAGS} -c ../pso/functions/FunctionManager.cpp \
	&& ${GCC} ${CMAKE_CXX_FLAGS} -c ../pso/main.cpp \
	&& ${GCC} ${CMAKE_CXX_FLAGS} Cec22.o main.o -o ${APP}.exe

run: release
	./${BUILDDIR}/${APP}.exe
