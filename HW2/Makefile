BUILDDIR=build
APP=app
RELEASE=-O3 -DNDEBUG
CMAKE_CXX_FLAGS=-std=c++20 -Wall -Wextra -Wpedantic -pthread
GCC=g++-11

builddir:
	if [[ -d "${BUILDDIR}" ]]; then echo "Directory ${BUILDDIR} already exists"; else echo "Creating ${BUILDDIR}"; mkdir ${BUILDDIR}; fi;

clean:
	rm -rf ${BUILDDIR}

format:
	clang-format -i ./pso/cec22/Cec22.cpp
	clang-format -i ./pso/cec22/Cec22.h
	clang-format -i ./pso/functions/FunctionManager.h
	clang-format -i ./pso/functions/FunctionManager.cpp
	clang-format -i ./pso/functions/CacheLayer.h
	clang-format -i ./pso/utils/Utils.h
	clang-format -i ./pso/utils/Utils.cpp
	clang-format -i ./pso/utils/Constants.h
	clang-format -i ./pso/pso/PSO.h
	clang-format -i ./pso/pso/PSO.cpp
	clang-format -i ./pso/main.cpp
	clang-format -i ./KDTree/KDTree.hpp
	clang-format -i ./KDTree/KDTree.cpp
	clang-format -i ./KDTree/tests/construction_time.cpp
	clang-format -i ./KDTree/tests/error_test.cpp
	clang-format -i ./KDTree/tests/toy_test.cpp

release: clean builddir
	cd ${BUILDDIR} \
	&& ${GCC} ${RELEASE} ${CMAKE_CXX_FLAGS} -c ../pso/cec22/Cec22.cpp \
	&& ${GCC} ${RELEASE} ${CMAKE_CXX_FLAGS} -c ../pso/functions/FunctionManager.cpp \
	&& ${GCC} ${RELEASE} ${CMAKE_CXX_FLAGS} -c ../KDTree/KDTree.cpp \
	&& ${GCC} ${RELEASE} ${CMAKE_CXX_FLAGS} -c ../pso/pso/PSO.cpp \
	&& ${GCC} ${RELEASE} ${CMAKE_CXX_FLAGS} -c ../pso/utils/Utils.cpp \
	&& ${GCC} ${RELEASE} ${CMAKE_CXX_FLAGS} -c ../pso/main.cpp \
	&& ${GCC} ${RELEASE} ${CMAKE_CXX_FLAGS} Cec22.o FunctionManager.o PSO.o main.o KDTree.o Utils.o -o ${APP}.exe

run: release
	./${BUILDDIR}/${APP}.exe

gdb: release # TODO: use debug
	gdb ./${BUILDDIR}/${APP}.exe
