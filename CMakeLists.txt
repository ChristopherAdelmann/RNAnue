cmake_minimum_required (VERSION 3.20)
project (RNAnue VERSION 0.1.1)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

message("Configuring Config.h")
message("Input file: ${CMAKE_CURRENT_SOURCE_DIR}/include/Config.h.in")
message("Output file: ${CMAKE_CURRENT_BINARY_DIR}/include/Config.h")

# configure header file to pass the version number to the source code
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/Config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/include/Config.h)

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

### SeqAn3
list (APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/submodules/seqan3/build_system")
find_package (seqan3 3.3.0 REQUIRED)

### OpenMP
find_package(OpenMP) 

### boost
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF)
set(Boost_DEBUG ON)

add_subdirectory(submodules/boost EXCLUDE_FROM_ALL)

##ViennaRNA
include(ExternalProject)
set(VIENNA_BUILD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/submodules/ViennaRNA/build)
file(MAKE_DIRECTORY ${VIENNA_BUILD_DIR})

ExternalProject_Add(ViennaRNA
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/submodules/ViennaRNA
    CONFIGURE_COMMAND 
    COMMAND "cd ${CMAKE_CURRENT_SOURCE_DIR}/submodules/ViennaRNA"
    COMMAND "autoreconf -i"
    COMMAND "./configure --prefix=${VIENNA_BUILD_DIR} --enable-lto --without-perl --without-python --without-doc --without-kinfold --without-forester --without-kinwalker --without-rnalocmin --without-cla --without-cla-pdf --without-doc-html --without-doc-pdf --without-svm --without-swig --without-tutorial --without-tutorial-html --without-tutorial-pdf"          
    BUILD_COMMAND "make"
    INSTALL_COMMAND "make install"
    TEST_COMMAND ""
)

#Compiler flags
#set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} "-lRNA -fno-lto")




file(GLOB SOURCES "./source/*.cpp")
include_directories(
./include
${VIENNA_BUILD_DIR}/include
) 

add_executable (RNAnue ${SOURCES}) #${VIENNA_SOURCES}

target_compile_features(RNAnue PRIVATE cxx_std_20)
target_link_libraries (RNAnue seqan3::seqan3)
target_link_libraries (RNAnue Boost::program_options Boost::filesystem Boost::property_tree Boost::foreach Boost::dynamic_bitset)
target_link_libraries(RNAnue ${VIENNA_BUILD_DIR}/build/lib/libRNA.a)

if(OpenMP_CXX_FOUND)
    target_link_libraries (RNAnue OpenMP::OpenMP_CXX)
endif()