cmake_minimum_required (VERSION 3.14) # CMake minimum version required. It is always good to ensure!
project (ParallelAssembleBenchmark) # name of your project

# ----- Find PZ. By default searches one directory up -----
find_package(NeoPZ REQUIRED HINTS ${CMAKE_SOURCE_DIR}/../neopz_install/ ${CMAKE_SOURCE_DIR}/neopz_install/)

file(GLOB_RECURSE PZ_SRC
     "../neopz/*.h"
     "../neopz/*.cpp"
)
add_custom_target( GlobalFiles SOURCES ${PZ_SRC})

add_executable(pbench main.cpp) # creates an executable called myExecutable using the file main.cpp
target_link_libraries(pbench NeoPZ::pz) # linking the executable myExecutable with the static lib fakeiset
