CMAKE_MINIMUM_REQUIRED(VERSION 3.16)

IF (NOT DEFINED PATEFIELD_VERSION_INFO)
    SET(PATEFIELD_VERSION_MAJOR 0)
    SET(PATEFIELD_VERSION_MINOR 1)
    SET(PATEFIELD_VERSION_PATCH 0)
    SET(PATEFIELD_VERSION_INFO "${PATEFIELD_VERSION_MAJOR}.${PATEFIELD_VERSION_MINOR}.${PATEFIELD_VERSION_PATCH}")
ENDIF ()

PROJECT(
    patefield
    VERSION ${PATEFIELD_VERSION_INFO}
    DESCRIPTION "Modern C++ implementation of Patefield Algorithm."
    LANGUAGES CXX
)

################################################################################
#                                   SETTINGS                                   #
################################################################################
OPTION(PATEFIELD_DEV_MODE OFF)
OPTION(PATEFIELD_ENABLE_DEBUG OFF)
OPTION(PATEFIELD_ENABLE_OPENMP OFF)
OPTION(PATEFIELD_BUILDING_WHEELS OFF)
OPTION(PATEFIELD_ENABLE_EXT_TESTS OFF)

IF (PATEFIELD_ENABLE_EXT_TESTS OR PATEFIELD_DEV_MODE)
    ADD_DEFINITIONS(-DPATEFIELD_BUILD_CPP_TESTS=TRUE)
ENDIF ()

SET(CMAKE_CXX_STANDARD 14)
SET(CMAKE_POSITION_INDEPENDENT_CODE ON)

IF (PATEFIELD_ENABLE_OPENMP)
    FIND_PACKAGE(OpenMP REQUIRED)
ENDIF ()

################################################################################
#                                     FLAGS                                    #
################################################################################
INCLUDE(GNUInstallDirs)

SET(CMAKE_MODULE_PATH_SAVED ${CMAKE_MODULE_PATH})
LIST(INSERT CMAKE_MODULE_PATH 0 "${PROJECT_SOURCE_DIR}/cmake")

SET(DEFAULT_BUILD_TYPE "Release")
# Set build type to Release if not specified
INCLUDE(BuildType)

IF (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    # error C2039: '_snprintf': is not a member of 'std'
    ADD_DEFINITIONS(-DHAVE_SNPRINTF)
ENDIF ()

########################################################################################################
#                                          DYNAMIC LIBRARY                                             #
########################################################################################################
ADD_LIBRARY(patefield SHARED ${PROJECT_SOURCE_DIR}/src/patefield.cpp)

IF (OpenMP_CXX_FOUND)
    MESSAGE(STATUS "patefield: Building with OpenMP support")
    TARGET_COMPILE_DEFINITIONS(patefield PUBLIC PATEFIELD_HAS_OPENMP_SUPPORT=TRUE)
    TARGET_LINK_LIBRARIES(patefield PRIVATE OpenMP::OpenMP_CXX)
ENDIF ()

TARGET_INCLUDE_DIRECTORIES(patefield PUBLIC "${PROJECT_SOURCE_DIR}/include")
TARGET_INCLUDE_DIRECTORIES(patefield PRIVATE "${PROJECT_SOURCE_DIR}/external/pcg-cpp/include")
TARGET_COMPILE_DEFINITIONS(patefield
    INTERFACE
    PATEFIELD_VERSION_MAJOR=${PATEFIELD_VERSION_MAJOR}
    PATEFIELD_VERSION_MINOR=${PATEFIELD_VERSION_MINOR}
    PATEFIELD_VERSION_PATCH=${PATEFIELD_VERSION_PATCH}
    PATEFIELD_VERSION_NAME=${PATEFIELD_VERSION_INFO}
)

IF(PATEFIELD_DEV_MODE)
    MESSAGE(STATUS "patefield: Building in development mode")
    TARGET_COMPILE_OPTIONS(patefield PRIVATE -Wall -Wextra -Wunused-variable -Wunused-const-variable)
ENDIF()

TARGET_COMPILE_FEATURES(patefield INTERFACE cxx_std_14)
TARGET_COMPILE_OPTIONS(patefield PRIVATE "$<$<CONFIG:RELEASE>:${PATEFIELD_RELEASE_OPTIONS}>")
SET_PROPERTY(TARGET patefield PROPERTY CXX_STANDARD_REQUIRED ON)
SET_PROPERTY(TARGET patefield PROPERTY CXX_STANDARD 14)
SET_PROPERTY(TARGET patefield PROPERTY POSITION_INDEPENDENT_CODE ON)

INSTALL(
    TARGETS patefield
    ARCHIVE DESTINATION ${PROJECT_SOURCE_DIR}
    LIBRARY DESTINATION ${PROJECT_SOURCE_DIR}
    RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}
)

ADD_LIBRARY(patefield::patefield ALIAS patefield)

########################################################################################################
#                                          HEADER TARGET                                               #
########################################################################################################
ADD_LIBRARY(patefield_headers INTERFACE)
ADD_LIBRARY(patefield::headers ALIAS patefield_headers)

TARGET_COMPILE_FEATURES(patefield_headers INTERFACE cxx_std_14)
TARGET_COMPILE_DEFINITIONS(patefield_headers
    INTERFACE
    PATEFIELD_VERSION_MAJOR=${PATEFIELD_VERSION_MAJOR}
    PATEFIELD_VERSION_MINOR=${PATEFIELD_VERSION_MINOR}
    PATEFIELD_VERSION_PATCH=${PATEFIELD_VERSION_PATCH}
    PATEFIELD_VERSION_NAME=${PATEFIELD_VERSION_INFO}
)

# Set required C++14 flags
SET_PROPERTY(TARGET patefield_headers PROPERTY CXX_STANDARD_REQUIRED ON)
SET_PROPERTY(TARGET patefield_headers PROPERTY CXX_STANDARD 14)
SET_PROPERTY(TARGET patefield_headers PROPERTY POSITION_INDEPENDENT_CODE ON)

# -- Clear cache --
UNSET(PATEFIELD_DEV_MODE CACHE)
UNSET(PATEFIELD_ENABLE_DEBUG CACHE)
UNSET(PATEFIELD_ENABLE_OPENMP CACHE)
UNSET(PATEFIELD_BUILDING_WHEELS CACHE)
UNSET(PATEFIELD_ENABLE_EXT_TESTS CACHE)
