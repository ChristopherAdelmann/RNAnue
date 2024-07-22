include(ExternalProject)

message(STATUS "ViennaRNA build dir: ${VIENNA_BUILD_DIR}")

find_library(VIENNA_RNA_LIBRARY NAMES RNA)

if(NOT VIENNA_RNA_LIBRARY)
  message(FATAL_ERROR "ViennaRNA library not found")
else()
  message(STATUS "ViennaRNA library found: ${VIENNA_RNA_LIBRARY}")
endif()