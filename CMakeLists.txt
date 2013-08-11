# Author: petter.strandmark@gmail.com (Petter Strandmark)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.0)
PROJECT(EASYIP C CXX)

OPTION(USE_MOSEK "Use Mosek solver." OFF)
OPTION(USE_CPLEX "Use Cplex solver." OFF)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/bin)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Default locations to search for on various platforms.
LIST(APPEND SEARCH_LIBS /usr/lib)
LIST(APPEND SEARCH_LIBS /usr/local/lib)
LIST(APPEND SEARCH_LIBS /usr/local/homebrew/lib) # Mac OS X
LIST(APPEND SEARCH_LIBS /opt/local/lib)

LIST(APPEND SEARCH_HEADERS /usr/include)
LIST(APPEND SEARCH_HEADERS /usr/local/include)
LIST(APPEND SEARCH_HEADERS /usr/local/homebrew/include) # Mac OS X
LIST(APPEND SEARCH_HEADERS /opt/local/include)

ENABLE_TESTING()

# Change the default build type from Debug to Release, while still
# supporting overriding the build type.
#
# The CACHE STRING logic here and elsewhere is needed to force CMake
# to pay attention to the value of these variables.
IF (NOT CMAKE_BUILD_TYPE)
  MESSAGE("-- No build type specified; defaulting to CMAKE_BUILD_TYPE=Release.")
  SET(CMAKE_BUILD_TYPE Release CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
    FORCE)
ENDIF (NOT CMAKE_BUILD_TYPE)

#
# gcc settings.
#
IF (CMAKE_COMPILER_IS_GNUCXX)
  MESSAGE("-- Enabling C++11 support for gcc.")

  # No warnings for C. The Meschach library contains really old code.
  # -fPIC for building a shared library.
  set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -w -fPIC")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-sign-compare -Wno-unused-parameter -std=c++0x -fPIC")
ENDIF (CMAKE_COMPILER_IS_GNUCXX)

#
# Clang settings
#
IF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  MESSAGE("-- Enabling C++11 support for clang.")

  # No warnings for C. The Meschach library contains really old code.
  # -fPIC for building a shared library.
  set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -w -fPIC")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-sign-compare -Wno-unused-parameter -std=c++11 -fPIC")
ENDIF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")

#
# MSVC settings
#
IF (MSVC)
  # Disable deprecation warning for standard functions.
  ADD_DEFINITIONS("/wd4996")
  # To get rid of annoying min and max macros if windows.h
  # is included.
  ADD_DEFINITIONS("-DNOMINMAX=1")
ENDIF (MSVC)


#
# Libraries
#
SET(COIN_LIB_DIR "" CACHE PATH "Coin (Cbc, Clp, ...) library directory")

FIND_LIBRARY(CBC_LIB NAMES Cbc libCbc PATHS ${COIN_LIB_DIR})

IF (NOT EXISTS ${CBC_LIB})
	MESSAGE(FATAL_ERROR
	        "Can't find Cbc. Please specify COIN_LIB_DIR")
ELSE (NOT EXISTS ${CBC_LIB})
	MESSAGE("-- Found Coin library (Cbc).")
ENDIF (NOT EXISTS ${CBC_LIB})

FIND_LIBRARY(CBCSOLVER_LIB NAMES CbcSolver libCbcSolver PATHS ${COIN_LIB_DIR})
FIND_LIBRARY(CGL_LIB NAMES Cgl libCgl PATHS ${COIN_LIB_DIR})
FIND_LIBRARY(CLP_LIB NAMES Clp libClp PATHS ${COIN_LIB_DIR})
FIND_LIBRARY(COINUTILS_LIB NAMES CoinUtils libCoinUtils PATHS ${COIN_LIB_DIR})
FIND_LIBRARY(OSI_LIB NAMES Osi libOsi PATHS ${COIN_LIB_DIR})
FIND_LIBRARY(OSI_CLP_LIB NAMES OsiClp libOsiClp PATHS ${COIN_LIB_DIR})

if (${USE_MOSEK})
	FIND_LIBRARY(OSI_MSK_LIB NAMES libOsiMsk OsiMsk PATHS ${COIN_LIB_DIR})
	FIND_PATH(OSI_MSK_INCLUDE NAMES coin/OsiMskSolverInterface.hpp PATHS ${COIN_INCLUDE_DIR} )

	SET(MOSEK_BIN_DIR "" CACHE PATH "Mosek binary directory.")
	FIND_LIBRARY(MOSEK_LIB NAMES mosek64_7_0 PATHS ${MOSEK_BIN_DIR})

	ADD_DEFINITIONS(-DHAS_MOSEK=1)
endif (${USE_MOSEK})

if (${USE_CPLEX})
	FIND_LIBRARY(OSI_CPX_LIB NAMES libOsiCpx OsiCpx PATHS ${COIN_LIB_DIR})
	FIND_PATH(OSI_CPX_INCLUDE NAMES coin/OsiCpxSolverInterface.hpp PATHS ${COIN_INCLUDE_DIR})

	SET(CPLEX_LIB_DIR "" CACHE PATH "Cplex library directory.")
	FIND_LIBRARY(CPLEX_LIB NAMES cplex1251 PATHS ${CPLEX_LIB_DIR})

	ADD_DEFINITIONS(-DHAS_CPLEX=1)
endif (${USE_CPLEX})

#
# Include directories
#
SET(COIN_INCLUDE_DIR "" CACHE PATH "Coin (Cbc, Clp, ...) include directory (where coin/ is)")

FIND_PATH(COIN_INCLUDE NAMES coin/CbcModel.hpp PATHS ${COIN_INCLUDE_DIR})
IF (NOT EXISTS ${COIN_INCLUDE})
	MESSAGE(FATAL_ERROR
	        "Can't find coin include directory. Please specify COIN_INCLUDE_DIR")
ELSE (NOT EXISTS ${COIN_INCLUDE})
	MESSAGE("-- Found Coin include directory.")
ENDIF (NOT EXISTS ${COIN_INCLUDE})

INCLUDE_DIRECTORIES(${COIN_INCLUDE})
INCLUDE_DIRECTORIES(${OSIMSK_INCLUDE})
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})
INCLUDE_DIRECTORIES(third-party/Catch)
FILE(GLOB EASYIP_HEADERS ${CMAKE_SOURCE_DIR}/*.h)

ADD_LIBRARY(easyip SHARED easy-ip.cpp easy-ip.h)

TARGET_LINK_LIBRARIES(easyip ${CBC_LIB} ${CBCSOLVER_LIB} ${CGL_LIB} ${CLP_LIB} ${COINUTILS_LIB} ${OSI_LIB} ${OSI_CLP_LIB})
if (${USE_MOSEK})
	TARGET_LINK_LIBRARIES(easyip ${MOSEK_LIB} ${OSI_MSK_LIB})
endif (${USE_MOSEK})
if (${USE_CPLEX})
	TARGET_LINK_LIBRARIES(easyip ${CPLEX_LIB} ${OSI_CPX_LIB})
endif (${USE_CPLEX})

INSTALL(TARGETS easyip
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

ADD_SUBDIRECTORY(examples)
ADD_SUBDIRECTORY(tests)
