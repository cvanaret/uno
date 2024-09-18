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

file (GLOB SOURCE_FILES "${CMAKE_SOURCE_DIR}/../../uno/*.cpp")
file (GLOB PYTHON_FILES "unopy.cpp")

find_package(pybind11 REQUIRED)

pybind11_add_module(unopy 
	${SOURCE_FILES}
	${PYTHON_FILES}
)

target_link_libraries(unopy PUBLIC)
