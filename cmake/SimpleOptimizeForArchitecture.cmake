SET(PATEFIELD_ARCHITECTURE_FLAGS "")

INCLUDE(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG(-march=native COMPILER_SUPPORTS_MARCH_NATIVE)
IF(COMPILER_SUPPORTS_MARCH_NATIVE)
    SET(PATEFIELD_ARCHITECTURE_FLAGS ${PATEFIELD_ARCHITECTURE_FLAGS} -march=native)
ENDIF()
CHECK_CXX_COMPILER_FLAG(-mtune=native COMPILER_SUPPORTS_MTUNE_NATIVE)
IF(COMPILER_SUPPORTS_MTUNE_NATIVE)
    SET(PATEFIELD_ARCHITECTURE_FLAGS ${PATEFIELD_ARCHITECTURE_FLAGS} -mtune=native)
ENDIF()
CHECK_CXX_COMPILER_FLAG(-ftree-vectorize COMPILER_SUPPORTS_FTREE)
IF(COMPILER_SUPPORTS_FTREE)
    SET(PATEFIELD_ARCHITECTURE_FLAGS ${PATEFIELD_ARCHITECTURE_FLAGS} -ftree-vectorize)
ENDIF()
CHECK_CXX_COMPILER_FLAG(-msse3 COMPILER_SUPPORTS_SSE3)
IF(COMPILER_SUPPORTS_SSE3)
    SET(PATEFIELD_ARCHITECTURE_FLAGS ${PATEFIELD_ARCHITECTURE_FLAGS} -msse3)
ENDIF()
CHECK_CXX_COMPILER_FLAG(-mavx COMPILER_SUPPORTS_MAVX)
IF(COMPILER_SUPPORTS_MAVX)
    SET(PATEFIELD_ARCHITECTURE_FLAGS ${PATEFIELD_ARCHITECTURE_FLAGS} -mavx)
ENDIF()
CHECK_CXX_COMPILER_FLAG(-mavx2 COMPILER_SUPPORTS_MAVX2)
IF(COMPILER_SUPPORTS_MAVX2)
    SET(PATEFIELD_ARCHITECTURE_FLAGS ${PATEFIELD_ARCHITECTURE_FLAGS} -mavx2)
ENDIF()
