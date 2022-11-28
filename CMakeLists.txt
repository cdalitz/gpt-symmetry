# set to the available number of cores
# when your compiler does not support OpenMP, it is automatacilly set to one
set(NUM_THREADS 4)

# finding some libraries on OSX and Windows only works with newer cmake version
if (APPLE AND MSVC)
  cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
else ()
  cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
endif ()

# for compilation with optimization
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

project(gptsymmetry)

if (MSVC)
  # with MS Visual C++ we must explicity switch on proper exception handling
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc")
  # Under windows CMP0074 must be set to new so <PackageName>_ROOT
  # environment variables can be used for find package
  cmake_policy(SET CMP0074 NEW)
endif (MSVC)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
#set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# libpng dependency
find_package(PNG REQUIRED)
include_directories(${PNG_INCLUDE_DIR})

# Check if more than one thread should be used
if(NOT ${NUM_THREADS} MATCHES 1)
  # check for OpenMP
  find_package(OpenMP)
  if (OPENMP_FOUND)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  else()
    # OpenMP is not supportet
    # reset NUM_THREADS to one
    set(NUM_THREADS 1)
  endif()
endif()

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNUM_THREADS=${NUM_THREADS}")

# all source files
set(SRC src/main.cpp src/image.cpp src/image_processing.cpp src/symmetry_transform.cpp)

# default target (created with "make")
add_executable (gptsymmetry ${SRC})
# target_compile_options(symmetry PUBLIC -Wall -Wextra -pedantic)

if (APPLE)
  target_link_libraries(gptsymmetry ${PNG_LIBRARY} OpenMP::OpenMP_CXX)
else ()
  target_link_libraries(gptsymmetry ${PNG_LIBRARY})
endif()

