CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(minkowski)

INCLUDE_DIRECTORIES(../libnormaliz)
INCLUDE_DIRECTORIES(..)

#Eigen library
INCLUDE_DIRECTORIES(../external/eigen)
INCLUDE_DIRECTORIES(../external)

include_directories("/usr/local/include/flint/flintxx")
include_directories("/usr/local/include/flint/flint")
include_directories("/usr/local/include/flint")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -lflint -lmpfr -lgmp")


INCLUDE_DIRECTORIES(${GMP_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${GMPXX_INCLUDE_DIR})

#find_package (FLINT)
#if (FLINT_FOUND)
#    include_directories(${FLINT_INCLUDE_DIRS})
    #target_link_libraries (helloworld ${BZIP2_LIBRARIES})
#endif (FLINT_FOUND)




ADD_EXECUTABLE (minkowski minkowski.cpp ../external/rref.hpp)

IF (NMZ_BUILD_STATIC)
    TARGET_LINK_LIBRARIES (minkowski normaliz-static)
    TARGET_LINK_LIBRARIES (minkowski ${SCIP_LIBRARY} ${COCOA_LIBRARY} ${GMPXX_STATIC_LIBRARIES} ${GMP_STATIC_LIBRARIES})
    #Static start and end
    SET_TARGET_PROPERTIES(minkowski PROPERTIES LINK_SEARCH_START_STATIC 1)
    SET_TARGET_PROPERTIES(minkowski PROPERTIES LINK_SEARCH_END_STATIC 1)
    # Don't try to modify the RPATH
    SET_TARGET_PROPERTIES(minkowski PROPERTIES INSTALL_RPATH "")
ELSE (NMZ_BUILD_STATIC)
    TARGET_LINK_LIBRARIES (minkowski normaliz-shared)
    TARGET_LINK_LIBRARIES (minkowski ${SCIP_LIBRARY} ${COCOA_LIBRARY} ${GMPXX_LIBRARIES} ${GMP_LIBRARIES})
ENDIF (NMZ_BUILD_STATIC)