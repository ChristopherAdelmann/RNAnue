cmake_minimum_required(VERSION 3.10)
project(segemehl)

find_library(HTSLIB htslib)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -pedantic -g -O3 -DSORTEDUNMAPPED -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DDBGNFO -DSHOWALIGN -DDBGLEVEL=0 -DPROGNFO")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} $(pkg-config --cflags htslib) -lm -lpthread -lz -lncurses -lform -lmenu")

include_directories(include)

file(GLOB sources "libs/*.c")

add_library(segemehl STATIC ${sources})
target_include_directories(segemehl PUBLIC include)