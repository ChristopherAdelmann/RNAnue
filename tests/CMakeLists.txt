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
                        ${CMAKE_CURRENT_SOURCE_DIR}/DetectTests.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/EvaluatedSplitRecordsTests.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/Analyze.cpp 
                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/DataTypes.cpp 
                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/annotation/FeatureAnnotator.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/annotation/FeatureParser.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/annotation/FeatureWriter.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/Utility.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/Detect.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/EvaluatedSplitRecords.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/CooptimalPairwiseAligner.cpp)

                        
add_dependencies(tests ViennaRNA)
target_include_directories (tests PRIVATE ${VIENNA_BUILD_DIR}/include)
target_link_libraries (tests PRIVATE ${VIENNA_BUILD_DIR}/lib/libRNA.a)

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