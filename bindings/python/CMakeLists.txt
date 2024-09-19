cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-O3")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

project(unopy)

include_directories("${CMAKE_SOURCE_DIR}/../../uno")

file (GLOB SOURCE_FILES ${CMAKE_SOURCE_DIR}/../../uno/tools/Options.cpp)
file (GLOB HEADER_FILES
		${CMAKE_SOURCE_DIR}/../../uno/linear_algebra/SparseVector.hpp
		${CMAKE_SOURCE_DIR}/../../uno/tools/Options.hpp)
file (GLOB PYTHON_FILES SparseVector.cpp Options.cpp unopy.cpp)

find_package(pybind11 REQUIRED)

pybind11_add_module(unopy 
	${SOURCE_FILES}
	${HEADER_FILES}
	${PYTHON_FILES}
)