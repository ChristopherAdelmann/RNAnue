# ----------------------------------------------------------------------------
# Testing
# ----------------------------------------------------------------------------
message(STATUS "Building tests")

include(FetchContent)
FetchContent_Declare(
    googletest
    DOWNLOAD_EXTRACT_TIMESTAMP true
    URL https://github.com/google/googletest/archive/1d17ea141d2c11b8917d2c7d029f1c4e2b9769b2.zip
)

FetchContent_MakeAvailable(googletest)

enable_testing()

add_executable(tests    ${CMAKE_CURRENT_SOURCE_DIR}/InteractionClusterTests.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/FeatureAnnotatorTests.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/Analyze.cpp 
                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/DataTypes.cpp 
                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/FeatureAnnotator.cpp 
                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/Utility.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/FeatureParser.cpp)

target_link_libraries(tests     PRIVATE 
                                GTest::gtest_main
                                seqan3::seqan3
                                Boost::program_options 
                                Boost::filesystem 
                                Boost::property_tree 
                                Boost::foreach 
                                Boost::dynamic_bitset
                                Boost::math)


include(GoogleTest)
gtest_discover_tests(tests)