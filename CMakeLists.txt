cmake_minimum_required(VERSION 3.25)

project(RNAnue VERSION 1.0.0)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-fopenmp")
set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# Configure header file to pass the version number to the source code
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/Config.hpp.in ${CMAKE_CURRENT_SOURCE_DIR}/include/Config.hpp)

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")


### SeqAn3
list(APPEND CMAKE_PREFIX_PATH ${CMAKE_CURRENT_SOURCE_DIR}/submodules/seqan3/build_system)
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
set(VIENNA_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/submodules/ViennaRNA/build)
file(MAKE_DIRECTORY ${VIENNA_BUILD_DIR})

add_custom_target(
    ViennaRNA ALL
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/submodules/ViennaRNA/configure
        CXX=$ENV{CXX}
        CC=$ENV{CC}
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
        --enable-universal-binary
    COMMAND make CXX=$ENV{CXX} CC=$ENV{CC}
    COMMAND make install
    WORKING_DIRECTORY ${VIENNA_BUILD_DIR}
)

### htslib
include(htslib)

### segemehl
# set(SEGEMEHL_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/submodules/segemehl/build)
# set(SEGEMEHL_BIN_DIR ${CMAKE_CURRENT_BINARY_DIR}/submodules/segemehl/bin)
# file(MAKE_DIRECTORY ${SEGEMEHL_BUILD_DIR})
# file(MAKE_DIRECTORY ${SEGEMEHL_BIN_DIR})

# set(ENV{PKG_CONFIG_PATH} ${CMAKE_CURRENT_BINARY_DIR}/submodules/htslib-install/lib/pkgconfig:ENV{PKG_CONFIG_PATH})

# add_custom_target(
#     segemehl ALL
#     COMMAND make all TARGETDIR=${SEGEMEHL_BIN_DIR} BUILDDIR=${SEGEMEHL_BUILD_DIR} CXX=$ENV{CXX} CC=$ENV{CC} PKG_CONFIG_PATH=${CMAKE_CURRENT_BINARY_DIR}/submodules/htslib-install/lib/pkgconfig
#     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/submodules/segemehl
# )

# add_dependencies(segemehl htslib)

add_subdirectory(submodules/segemehl)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

file(GLOB SOURCES CONFIGURE_DEPENDS "./src/*.cpp")
list(APPEND SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/submodules/bam_sort/bam_sort.c
                    ${CMAKE_CURRENT_SOURCE_DIR}/submodules/bam_sort/sam_opts.c
                    ${CMAKE_CURRENT_SOURCE_DIR}/submodules/bam_sort/bedidx.c
                    ${CMAKE_CURRENT_SOURCE_DIR}/submodules/bam_sort/bam.c
                    ${CMAKE_CURRENT_SOURCE_DIR}/submodules/bam_sort/sam_utils.c)

include_directories(include ${CMAKE_CURRENT_SOURCE_DIR}/submodules/bam_sort) 

add_executable (RNAnue ${SOURCES} ${SMATOOLS_SOURCES})

if(NOT HTSlib_FOUND)
    add_dependencies(RNAnue htslib)
endif()

target_link_libraries(RNAnue PRIVATE ${HTSlib_LIBRARIES})
target_include_directories(RNAnue PRIVATE ${HTSlib_INCLUDE_DIRS} ${HTSlib_INCLUDE_DIRS}/htslib)

add_dependencies(RNAnue ViennaRNA)
target_include_directories (RNAnue PRIVATE ${VIENNA_BUILD_DIR}/include)
target_link_libraries (RNAnue PRIVATE ${VIENNA_BUILD_DIR}/lib/libRNA.a)

target_link_libraries (RNAnue PRIVATE seqan3::seqan3 Boost::program_options Boost::filesystem Boost::property_tree Boost::foreach Boost::dynamic_bitset segemehl)

if(OpenMP_CXX_FOUND)
    target_link_libraries (RNAnue PRIVATE OpenMP::OpenMP_CXX)
endif()

install(TARGETS RNAnue RUNTIME DESTINATION bin)

set(CPACK_PACKAGE_NAME RNAnue)
set(CPACK_PACKAGE_VENDOR "RNA Biology & Bioinformatics; Institute of Biomedical Genetics; University of Stuttgart")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "RNAnue: efficient data analysis for RNA-RNA interactomics")
set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})
set(CPACK_VERBATIM_VARIABLES TRUE)
set(CPACK_PACKAGE_INSTALL_DIRECTORY ${CPACK_PACKAGE_NAME})

set(CPACK_RESOURCE_FILE_LICENSE ${CMAKE_CURRENT_SOURCE_DIR}/LICENSE)
#set(CPACK_RESOURCE_FILE_README ${CMAKE_CURRENT_SOURCE_DIR}/README.md)

if(APPLE)
    message(STATUS "Building for MacOS")
    set(CPACK_GENERATOR TGZ productbuild)
elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    message(STATUS "Building for Linux")
    set(CPACK_GENERATOR TGZ RPM)
else()
    message(STATUS "Building for Other OS")
    set(CPACK_GENERATOR TGZ)
endif()

include(CPack)

cpack_add_component(RNAnue
    DISPLAY_NAME "RNAnue"
    DESCRIPTION "RNAnue: efficient data analysis for RNA-RNA interactomics"
    REQUIRED
)