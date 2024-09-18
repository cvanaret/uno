# Copyright (c) 2018-2024 Charlie Vanaret
# Licensed under the MIT license. See LICENSE file in the project directory for details.

cmake_minimum_required(VERSION 3.7)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
	cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

# define the project name
project(Uno VERSION 1.0
			DESCRIPTION "Uno (Unifying Nonlinear Optimization)"
            LANGUAGES CXX C Fortran)

# set C++17 and enable other languages
set(CMAKE_CXX_STANDARD 17)
enable_language(CXX C Fortran)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wnon-virtual-dtor -pedantic -Wunused-value -Wconversion -Wmaybe-uninitialized")
set(CMAKE_CXX_FLAGS_DEBUG "-pg")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG") # disable asserts

# optional Gtest
option(WITH_GTEST "Enable GoogleTest" OFF)
message(STATUS "GoogleTest: WITH_GTEST=${WITH_GTEST}")

# directories
set(DIRECTORIES uno)

# source files
file(GLOB UNO_SOURCE_FILES
    uno/Uno.cpp
    uno/ingredients/globalization_mechanism/*.cpp
    uno/ingredients/globalization_strategy/*.cpp
    uno/ingredients/globalization_strategy/switching_methods/*.cpp
    uno/ingredients/globalization_strategy/switching_methods/filter_methods/*.cpp
    uno/ingredients/globalization_strategy/switching_methods/filter_methods/filters/*.cpp
    uno/ingredients/globalization_strategy/switching_methods/funnel_methods/*.cpp
    uno/ingredients/constraint_relaxation_strategy/*.cpp
    uno/ingredients/subproblem/*.cpp
    uno/ingredients/subproblem/inequality_constrained_methods/*.cpp
    uno/ingredients/subproblem/interior_point_methods/*.cpp
    uno/model/*.cpp
    uno/optimization/*.cpp
    uno/preprocessing/*.cpp
    uno/tools/*.cpp
)

#########################
# external dependencies #
#########################
set(LIBRARIES "")

# function that links an existing library to Uno
function(link_to_uno library_name library_path)
   # add the library
   set(LIBRARIES ${LIBRARIES} ${library_path} PARENT_SCOPE)
   # add a preprocessor definition
   string(TOUPPER ${library_name} library_name_upper)
   add_definitions("-D HAS_${library_name_upper}")
   # include the corresponding directory
   get_filename_component(directory ${library_path} DIRECTORY)
   set(DIRECTORIES ${DIRECTORIES} ${directory} PARENT_SCOPE)
   message(STATUS "Library ${library_name} was found.")
endfunction()

# detect optional libraries: ma57 metis bqpd
find_library(ma57 NAMES ma57 hsl)
if(${ma57} STREQUAL "ma57-NOTFOUND")
   message(WARNING "Optional library MA57 was not found.")
else()
   list(APPEND UNO_SOURCE_FILES uno/solvers/linear/MA57Solver.cpp)
   link_to_uno(ma57 ${ma57})
endif()

find_library(metis metis)
if(${metis} STREQUAL "metis-NOTFOUND")
   message(WARNING "Optional library METIS was not found.")
else()
   link_to_uno(metis ${metis})
endif()

find_library(bqpd bqpd)
if(${bqpd} STREQUAL "bqpd-NOTFOUND")
   message(WARNING "Optional library BQPD was not found.")
else()
   list(APPEND UNO_SOURCE_FILES uno/solvers/QP/BQPDSolver.cpp uno/solvers/QP/wdotd.f)
   link_to_uno(bqpd ${bqpd})
endif()

# detection of required libraries
find_library(dl dl REQUIRED)
list(APPEND LIBRARIES ${dl})

find_library(blas blas REQUIRED)
list(APPEND LIBRARIES ${blas})

find_library(lapack lapack REQUIRED)
list(APPEND LIBRARIES ${lapack})

find_package(OpenMP REQUIRED)
list(APPEND LIBRARIES OpenMP::OpenMP_CXX)

###############
# Uno library #
###############
add_library(uno STATIC ${UNO_SOURCE_FILES})
target_include_directories(uno PUBLIC ${DIRECTORIES})

# link the libraries
target_link_libraries(uno PUBLIC ${LIBRARIES})

# copy the option file
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/uno.options DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

######################
# optional AMPL main #
######################
find_library(amplsolver amplsolver)
if(${amplsolver} STREQUAL "amplsolver-NOTFOUND")
   message(WARNING "Optional library amplsolver (ASL) was not found.")
else()
   message(STATUS "Library amplsolver was found.")
   add_executable(uno_ampl bindings/AMPL/AMPLModel.cpp bindings/AMPL/uno_ampl.cpp)
   
   target_link_libraries(uno_ampl PUBLIC uno ${amplsolver})
   add_definitions("-D HAS_AMPLSOLVER")
   # include the corresponding directory
   get_filename_component(directory ${amplsolver} DIRECTORY)
   include_directories(${directory})
endif()

##################################
# optional GoogleTest unit tests #
##################################
if(WITH_GTEST)
    find_package(GTest CONFIG REQUIRED)
    if(NOT ${GTest}_DIR STREQUAL "${GTest}-NOTFOUND")
        file(GLOB TESTS_UNO_SOURCE_FILES
            unotest/*.cpp
        )
        add_executable(run_unotest ${TESTS_UNO_SOURCE_FILES})
        target_link_libraries(run_unotest PUBLIC GTest::gtest uno)
    endif()
endif()
