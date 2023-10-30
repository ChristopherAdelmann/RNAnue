cmake_minimum_required (VERSION 3.20)
project (RNAnue VERSION 0.1.0)
enable_testing()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

option(BUILD_STAR "Install STAR mapping tool. If OFF STAR should be in $PATH. Default is ON." ON)

execute_process(COMMAND git submodule init submodules/boost submodules/seqan3 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

if (BUILD_STAR)
    execute_process(COMMAND git submodule init submodules/STAR WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif()

execute_process(COMMAND git submodule update --recursive WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

message("Detected system: ${CMAKE_SYSTEM_NAME}")

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin" AND BUILD_STAR)
    execute_process(COMMAND STAR --version WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/submodules/STAR/MacOSX_x86_64/STAR RESULT_VARIABLE STAR_OUTPUT)
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Linux" AND BUILD_STAR)
    execute_process(COMMAND STAR --version WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/submodules/STAR/bin/Linux_x86_64/STAR RESULT_VARIABLE STAR_OUTPUT)
else()
    execute_process(COMMAND STAR --version RESULT_VARIABLE STAR_OUTPUT)
endif()

if(STAR_OUTPUT EQUAL "1")
    message(FATAL_ERROR "Could not find compatible STAR executable. Please install STAR seperatly and set BUILD_STAR to OFF.")
else()
    message("Compatible STAR version was found.")
endif()

message("Configuring Config.h")
message("Input file: ${CMAKE_CURRENT_SOURCE_DIR}/include/Config.h.in")
message("Output file: ${CMAKE_CURRENT_BINARY_DIR}/include/Config.h")

# configure header file to pass the version number to  the source code
#configure_file(./include/Config.h.in ./include/Config.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/Config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/include/Config.h)

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

### SeqAn3
list (APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/submodules/seqan3/build_system")
find_package (seqan3 3.0.2 REQUIRED)

### OpenMP
find_package(OpenMP)

### boost
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF)
set(Boost_DEBUG ON)

add_subdirectory(submodules/boost EXCLUDE_FROM_ALL)

file(GLOB SOURCES "./source/*.cpp")
include_directories(./include)

add_executable (RNAnue ${SOURCES} ${STAR_SOURCES})

target_compile_features(RNAnue PRIVATE cxx_std_20)
target_link_libraries (RNAnue seqan3::seqan3)
target_link_libraries (RNAnue Boost::program_options Boost::filesystem Boost::property_tree Boost::foreach Boost::dynamic_bitset)

if(OpenMP_CXX_FOUND)
    target_link_libraries (RNAnue OpenMP::OpenMP_CXX)
endif()