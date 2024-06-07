include(ExternalProject)

set(zlib_PREFIX ${CMAKE_BINARY_DIR}/submodules/zlib-prefix)
set(zlib_INSTALL ${CMAKE_BINARY_DIR}/submodules/zlib-install)

ExternalProject_Add(
  zlib
  PREFIX ${zlib_PREFIX}
  GIT_REPOSITORY "https://github.com/madler/zlib.git"
  GIT_TAG "v1.3.1"
  UPDATE_COMMAND ""
  # BUILD_IN_SOURCE 1 CONFIGURE_COMMAND ${zlib_PREFIX}/src/zlib/configure --prefix=${zlib_INSTALL}
  # --static
  INSTALL_DIR ${zlib_INSTALL}
  CMAKE_ARGS -DCMAKE_C_COMPILER=$ENV{CC}-DCMAKE_CXX_COMPILER=$ENV{CXX}
             -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DCMAKE_INSTALL_PREFIX=${zlib_INSTALL}
             -DCMAKE_MACOSX_RPATH=0
)

include_directories(${zlib_INSTALL}/include)
set(zlib_LIBRARIES ${zlib_INSTALL}/lib/libz.a)