cmake_minimum_required(VERSION 2.8)
project(envtools)

if(NOT MSVC)
set(CMAKE_CXX_FLAGS "--std=c++11  -fopenmp -Wall -Wpedantic")
endif(NOT MSVC)

set(IMAGESRC
envremap.cpp
envremap.h
)

add_library(envlib ${IMAGESRC})

option(LINK_LIB_TIFF "Link with libtiff" OFF)
option(LINK_LIB_OPENEXR "Link with OpenEXR" OFF)

#set(TIFF_LIBRARY_LOCATION  "/usr/local/lib" CACHE PATH "FFTW Library Directory")
#set(TIFF_INCLUDE_PATH "/usr/local/include" CACHE PATH "FFTW Include Directory")

if(LINK_LIB_TIFF)
find_package(TIFF REQUIRED)
IF(TIFF_FOUND)
add_definitions(-DTIFF_SUPPORT)
set(DEPENDENCY_INCLUDE ${DEPENDENCY_INCLUDE} ${TIFF_INCLUDE_DIR})
set(DEPENDENCY_LIB ${DEPENDENCY_LIB} ${TIFF_LIBRARIES})
endif(TIFF_FOUND)
endif()


if(LINK_LIB_OPENEXR)

find_package(OpenEXR REQUIRED)

if(OPENEXR_FOUND)
add_definitions(-DEXR_SUPPORT)
set(DEPENDENCY_INCLUDE ${DEPENDENCY_INCLUDE} ${OPENEXR_INCLUDE_PATHS})
set(DEPENDENCY_LIB ${DEPENDENCY_LIB} ${OPENEXR_LIBRARIES})
endif()

endif()


add_executable(envtool apps/envtool.cpp envremap.h image/inout.h image/inout.cpp)
target_link_libraries(envtool envlib ${DEPENDENCY_LIB})

add_executable(envconv apps/envconv.cpp envremap.h image/inout.h image/inout.cpp)
target_link_libraries(envconv envlib ${DEPENDENCY_LIB})

add_executable(sherror apps/sherror.cpp)
add_executable(shvisualize apps/shimage.cpp image/inout.cpp image/inout.h)
target_link_libraries(shvisualize envlib ${DEPENDENCY_LIB})
add_executable(shtrans apps/shtrans.cpp image/inout.cpp image/inout.h)
target_link_libraries(shtrans envlib ${DEPENDENCY_LIB})

include_directories(${DEPENDENCY_INCLUDE})