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

set(UNO_DIRECTORY ${CMAKE_SOURCE_DIR}/../../uno)
file (GLOB SOURCE_FILES ${UNO_DIRECTORY}/tools/Options.cpp)
file (GLOB HEADER_FILES
		${UNO_DIRECTORY}/linear_algebra/SparseVector.hpp
		${UNO_DIRECTORY}/linear_algebra/RectangularMatrix.hpp
		${UNO_DIRECTORY}/tools/Options.hpp)
file (GLOB PYTHON_FILES SparseVector.cpp RectangularMatrix.cpp Options.cpp unopy.cpp)

find_package(pybind11 REQUIRED)

pybind11_add_module(unopy 
	${SOURCE_FILES}
	${HEADER_FILES}
	${PYTHON_FILES}
)