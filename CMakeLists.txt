# Copyright (c) 2018-2024 Charlie Vanaret
# Licensed under the MIT license. See LICENSE file in the project directory for details.

cmake_minimum_required(VERSION 3.7)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
	cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

# define the project name
project(Uno VERSION 1.3.0
        DESCRIPTION "Uno (Unifying Nonconvex Optimization)"
        LANGUAGES CXX C Fortran)
set(CMAKE_CXX_STANDARD 17)

include(FortranCInterface)
FortranCInterface_VERIFY(CXX)


FortranCInterface_HEADER(${CMAKE_BINARY_DIR}/include/fortran_interface.h
                            MACRO_NAMESPACE "FC_"
                            SYMBOL_NAMESPACE "FC_")

# directories
set(DIRECTORIES uno ${CMAKE_BINARY_DIR}/include)


set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wnon-virtual-dtor -pedantic -Wunused-value -Wconversion")
set(CMAKE_CXX_FLAGS_DEBUG "-pg")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG") # disable asserts
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
   SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wmaybe-uninitialized")
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_CURRENT_SOURCE_DIR}/cmake-library/finders)

# optional Gtest
option(WITH_GTEST "Enable GoogleTest" OFF)
message(STATUS "GoogleTest: WITH_GTEST=${WITH_GTEST}")


# source files
file(GLOB UNO_SOURCE_FILES
   uno/Uno.cpp
   uno/ingredients/constraint_relaxation_strategies/*.cpp
   uno/ingredients/globalization_mechanisms/*.cpp
   uno/ingredients/globalization_strategies/*.cpp
   uno/ingredients/globalization_strategies/switching_methods/*.cpp
   uno/ingredients/globalization_strategies/switching_methods/filter_methods/*.cpp
   uno/ingredients/globalization_strategies/switching_methods/filter_methods/filters/*.cpp
   uno/ingredients/globalization_strategies/switching_methods/funnel_methods/*.cpp
   uno/ingredients/hessian_models/*.cpp
   uno/ingredients/inequality_handling_methods/*.cpp
   uno/ingredients/inequality_handling_methods/inequality_constrained_methods/*.cpp
   uno/ingredients/inequality_handling_methods/interior_point_methods/*.cpp
   uno/ingredients/subproblem_solvers/*.cpp
   uno/ingredients/subproblems/*.cpp
   uno/model/*.cpp
   uno/optimization/*.cpp
   uno/options/*.cpp
   uno/preprocessing/*.cpp
   uno/reformulation/*.cpp
   uno/tools/*.cpp
)

# unit test source files
file(GLOB TESTS_UNO_SOURCE_FILES
   unotest/unit_tests/unotest.cpp
   unotest/unit_tests/CollectionAdapterTests.cpp
   unotest/unit_tests/ConcatenationTests.cpp
   unotest/unit_tests/COOSparseStorageTests.cpp
   unotest/unit_tests/CSCSparseStorageTests.cpp
   unotest/unit_tests/MatrixVectorProductTests.cpp
   unotest/unit_tests/RangeTests.cpp
   unotest/unit_tests/ScalarMultipleTests.cpp
   unotest/unit_tests/SparseVectorTests.cpp
   unotest/unit_tests/SumTests.cpp
   unotest/unit_tests/VectorTests.cpp
   unotest/unit_tests/VectorViewTests.cpp
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

# HSL or MA57
find_library(HSL hsl)
if(HSL)
   link_to_uno(hsl ${HSL})
else()
   message(WARNING "Optional library HSL was not found.")
   find_library(MA57 ma57)
   if(MA57)
      link_to_uno(ma57 ${MA57})
   else()
      message(WARNING "Optional library MA57 was not found.")
   endif()
   find_library(MA27 ma27)
   if(MA27)
      link_to_uno(ma27 ${MA27})
   else()
      message(WARNING "Optional library MA27 was not found.")
   endif()
endif()
if(HSL OR MA57)
   list(APPEND UNO_SOURCE_FILES uno/ingredients/subproblem_solvers/MA57/MA57Solver.cpp)
   list(APPEND TESTS_UNO_SOURCE_FILES unotest/functional_tests/MA57SolverTests.cpp)
   
   find_package(BLAS REQUIRED)
   list(APPEND LIBRARIES ${BLAS_LIBRARIES})
endif()
if(HSL OR MA27)
   list(APPEND UNO_SOURCE_FILES uno/ingredients/subproblem_solvers/MA27/MA27Solver.cpp)
   list(APPEND TESTS_UNO_SOURCE_FILES unotest/functional_tests/MA27SolverTests.cpp)

   find_package(BLAS REQUIRED)
   list(APPEND LIBRARIES ${BLAS_LIBRARIES})
endif()

# METIS
find_library(METIS metis)
if(NOT METIS)
   message(WARNING "Optional library METIS was not found.")
else()
   link_to_uno(metis ${METIS})
endif()

# BQPD
find_library(BQPD bqpd)
if(NOT BQPD)
   message(WARNING "Optional library BQPD was not found.")
else()
   list(APPEND UNO_SOURCE_FILES uno/ingredients/subproblem_solvers/BQPD/BQPDSolver.cpp)
   list(APPEND TESTS_UNO_SOURCE_FILES unotest/functional_tests/BQPDSolverTests.cpp)
   link_to_uno(bqpd ${BQPD})
endif()

# HiGHS
find_package(HIGHS)
if(NOT HIGHS)
   message(WARNING "Optional library HiGHS was not found.")
else()
   list(APPEND UNO_SOURCE_FILES uno/ingredients/subproblem_solvers/HiGHS/HiGHSSolver.cpp)
   list(APPEND TESTS_UNO_SOURCE_FILES unotest/functional_tests/HiGHSSolverTests.cpp)
   link_to_uno(highs ${HIGHS})
   list(APPEND LIBRARIES highs::highs)
endif()

# MUMPS
find_package(MUMPS)
if(NOT MUMPS_LIBRARY)
   message(WARNING "Optional library MUMPS was not found.")
else()
   list(APPEND UNO_SOURCE_FILES uno/ingredients/subproblem_solvers/MUMPS/MUMPSSolver.cpp)
   list(APPEND TESTS_UNO_SOURCE_FILES unotest/functional_tests/MUMPSSolverTests.cpp)
   list(APPEND LIBRARIES ${MUMPS_LIBRARY} ${MUMPS_COMMON_LIBRARY} ${MUMPS_PORD_LIBRARY})
   
   list(APPEND DIRECTORIES ${MUMPS_INCLUDE_DIR})

   if(NOT MUMPS_MPISEQ_LIBRARY)
      # parallel
      add_definitions("-D MUMPS_PARALLEL")
      find_package(MPI REQUIRED)
      list(APPEND LIBRARIES MPI::MPI_CXX MPI::MPI_Fortran)
      add_definitions("-D HAS_MPI")

      find_package(BLACS REQUIRED)
      list(APPEND LIBRARIES ${BLACS_LIBRARY})
      list(APPEND DIRECTORIES ${BLACS_INCLUDE_DIRS})

      find_package(ScaLAPACK REQUIRED)
      list(APPEND LIBRARIES ${ScaLAPACK_LIBRARY})
      list(APPEND DIRECTORIES ${ScaLAPACK_INCLUDE_DIRS})
   else()
      # link dummy parallel library (provided by MUMPS distribution)
      link_to_uno(MUMPS_MPISEQ_LIBRARY ${MUMPS_MPISEQ_LIBRARY})
   endif()
   
   find_package(METIS REQUIRED)
   list(APPEND LIBRARIES ${METIS_LIBRARY})
   list(APPEND DIRECTORIES ${METIS_INCLUDE_DIRS})
   
   find_package(BLAS REQUIRED)
   list(APPEND LIBRARIES ${BLAS_LIBRARIES})
   
   find_package(OpenMP REQUIRED)
   list(APPEND LIBRARIES OpenMP::OpenMP_CXX)
   
   add_definitions("-D HAS_MUMPS")
endif()

###############
# Uno library #
###############
add_library(uno STATIC ${UNO_SOURCE_FILES})
set_property(TARGET uno PROPERTY POSITION_INDEPENDENT_CODE ON)
target_include_directories(uno PUBLIC ${DIRECTORIES})
target_link_libraries(uno PUBLIC ${LIBRARIES}
   # link Fortran compiler even when static libraries are linked
   ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES} 
)

######################
# optional AMPL main #
######################
find_library(AMPLSOLVER amplsolver)
if(NOT AMPLSOLVER)
   message(WARNING "Optional library amplsolver (ASL) was not found.")
else()
   message(STATUS "Library amplsolver was found.")
   add_executable(uno_ampl bindings/AMPL/AMPLModel.cpp bindings/AMPL/AMPLUserCallbacks.cpp bindings/AMPL/uno_ampl.cpp)
   
   target_link_libraries(uno_ampl PUBLIC uno ${AMPLSOLVER} ${CMAKE_DL_LIBS})
   add_definitions("-D HAS_AMPLSOLVER")
   # include the corresponding directory
   get_filename_component(directory ${AMPLSOLVER} DIRECTORY)
   include_directories(${directory})
endif()

##################################
# optional GoogleTest unit tests #
##################################
if(WITH_GTEST)
   find_package(GTest CONFIG REQUIRED)
   if(NOT GTest_DIR)
      message(WARNING "Optional library GTest was not found.")
   else()
      add_executable(run_unotest ${TESTS_UNO_SOURCE_FILES})
      target_link_libraries(run_unotest PUBLIC GTest::gtest uno)
   endif()
endif()
