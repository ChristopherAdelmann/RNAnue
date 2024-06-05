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
include (FindPkgConfig)
pkg_check_modules (ViennaRNA RNAlib2)

if (NOT ViennaRNA_FOUND)
    include (ExternalProject)
    ExternalProject_Add (
            "ViennaRNA"
            PREFIX "${CMAKE_CURRENT_BINARY_DIR}/ext/ViennaRNA"
            LOG_DOWNLOAD 1
            DOWNLOAD_NO_PROGRESS 1
            TIMEOUT 240
            URL https://github.com/ViennaRNA/ViennaRNA/releases/download/v2.4.17/ViennaRNA-2.4.17.tar.gz
            URL_HASH MD5=530b5ce7caf0e7cf58535cad22265215
            LOG_CONFIGURE 1
            USES_TERMINAL_CONFIGURE 1
            CONFIGURE_COMMAND <SOURCE_DIR>/configure CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER}
                --prefix=<INSTALL_DIR> --without-swig --without-kinfold --disable-simd
                --without-forester --without-rnalocmin --without-doc --without-tutorial
            LOG_BUILD 1
            USES_TERMINAL_BUILD 1
            BUILD_IN_SOURCE 1
            BUILD_COMMAND make -j
            LOG_INSTALL 1
    )
    ExternalProject_Get_Property (ViennaRNA INSTALL_DIR)
endif()

# set(VIENNA_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/submodules/ViennaRNA/build)
# file(MAKE_DIRECTORY ${VIENNA_BUILD_DIR})

# add_custom_target(
#     BuildViennaRNA ALL
#     COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/submodules/ViennaRNA/configure
#         --prefix=${VIENNA_BUILD_DIR}
#         --enable-lto
#         --without-perl
#         --without-python
#         --without-doc
#         --without-kinfold
#         --without-forester
#         --without-kinwalker
#         --without-cla
#         --without-cla-pdf
#         --without-doc-html
#         --without-doc-pdf
#         --without-tutorial
#         --without-tutorial-html
#         --without-tutorial-pdf
#         --enable-universal-binary
#         CXX=$ENV{CXX}
#         CC=$ENV{CC}
#     COMMAND make
#     COMMAND make install
#     WORKING_DIRECTORY ${VIENNA_BUILD_DIR}
# )

# add_library(ViennaRNA STATIC IMPORTED GLOBAL)
# set_property(TARGET ViennaRNA PROPERTY
#     IMPORTED_LOCATION ${VIENNA_BUILD_DIR}/src/ViennaRNA/libRNA.a
# )

# add_dependencies(ViennaRNA BuildViennaRNA)

### samtools
execute_process(
    COMMAND autoheader
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/submodules/samtools
)

execute_process(
    COMMAND autoconf -Wno-syntax
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/submodules/samtools
)

execute_process(
    COMMAND ./configure --prefix=${CMAKE_CURRENT_SOURCE_DIR}/submodules/samtools/build
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/submodules/samtools
)

execute_process(
    COMMAND make
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/submodules/samtools
)

execute_process(
    COMMAND make install
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/submodules/samtools
)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

file(GLOB SOURCES "./src/*.cpp")
include_directories(include) 

add_executable (RNAnue ${SOURCES})

target_link_libraries (RNAnue PRIVATE seqan3::seqan3)
target_link_libraries (RNAnue PRIVATE Boost::program_options Boost::filesystem Boost::property_tree Boost::foreach Boost::dynamic_bitset)

if (ViennaRNA_FOUND)
    target_include_directories (RNAnue SYSTEM PUBLIC ${ViennaRNA_INCLUDE_DIRS})
    target_link_libraries (RNAnue PUBLIC ${ViennaRNA_LDFLAGS})
else ()
    add_dependencies (RNAnue ViennaRNA)
    target_include_directories (RNAnue SYSTEM PUBLIC ${INSTALL_DIR}/include)
    target_link_libraries (RNAnue PUBLIC ${INSTALL_DIR}/lib/libRNA.a)
endif ()

if(OpenMP_CXX_FOUND)
    target_link_libraries (RNAnue PRIVATE OpenMP::OpenMP_CXX)
endif()

target_compile_options(RNAnue PRIVATE -Wno-interference-size)