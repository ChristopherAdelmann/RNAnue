cmake_minimum_required (VERSION 3.25)
project (RNAnue VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

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

add_executable (RNAnue ${SOURCES})

target_compile_features(RNAnue PRIVATE cxx_std_20)
target_link_libraries (RNAnue seqan3::seqan3)
target_link_libraries (RNAnue Boost::program_options Boost::filesystem Boost::property_tree Boost::foreach Boost::dynamic_bitset)

if(OpenMP_CXX_FOUND)
    target_link_libraries (RNAnue OpenMP::OpenMP_CXX)
endif()