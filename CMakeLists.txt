cmake_minimum_required(VERSION 3.25)

project(RNAnue VERSION 0.2.0)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS -fopenmp)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# Configure header file to pass the version number to the source code
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/Config.hpp.in ${CMAKE_CURRENT_SOURCE_DIR}/include/Config.hpp)

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

### SeqAn3
list(APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/submodules/seqan3/build_system")
find_package(seqan3 3.3.0 REQUIRED)

### OpenMP
find_package(OpenMP) 

### Boost
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF)
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(Boost_DEBUG ON)
endif()

add_subdirectory(submodules/boost EXCLUDE_FROM_ALL)

### ViennaRNA
if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(MACOSX_UNIVERSAL_BINARY --enable-universal-binary)
endif()

set(VIENNA_BUILD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/submodules/ViennaRNA/ViennaRNA-2.6.4/build)
file(MAKE_DIRECTORY ${VIENNA_BUILD_DIR})

    execute_process(
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/submodules/ViennaRNA/ViennaRNA-2.6.4/configure
        --prefix=${VIENNA_BUILD_DIR}
        --enable-lto
        --without-perl
        --without-python
        --without-doc
        --without-kinfold
        --without-forester
        --without-kinwalker
        --without-cla
        --without-cla-pdf
        --without-doc-html
        --without-doc-pdf
        --without-tutorial
        --without-tutorial-html
        --without-tutorial-pdf
        CXX=$ENV{CXX}
        CC=$ENV{CC}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/submodules/ViennaRNA/ViennaRNA-2.6.4
)

execute_process(
    COMMAND make
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/submodules/ViennaRNA/ViennaRNA-2.6.4
)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

file(GLOB SOURCES "./src/*.cpp")
include_directories(
./include
${VIENNA_BUILD_DIR}/include
) 

add_executable (RNAnue ${SOURCES})

target_link_libraries (RNAnue seqan3::seqan3)
target_link_libraries (RNAnue Boost::program_options Boost::filesystem Boost::property_tree Boost::foreach Boost::dynamic_bitset)
target_link_libraries(RNAnue ${VIENNA_BUILD_DIR}/lib/libRNA.a)

if(OpenMP_CXX_FOUND)
    target_link_libraries (RNAnue OpenMP::OpenMP_CXX)
endif()

target_compile_options(RNAnue PRIVATE -Wno-interference-size)