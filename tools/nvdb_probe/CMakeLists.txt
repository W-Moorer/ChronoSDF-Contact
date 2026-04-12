cmake_minimum_required(VERSION 3.20)

project(nvdb_probe LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(OpenVDB CONFIG REQUIRED)

add_executable(nvdb_probe
    src/main.cpp
)

target_link_libraries(nvdb_probe PRIVATE OpenVDB::openvdb)

if(MSVC)
    target_compile_options(nvdb_probe PRIVATE /W4 /permissive-)
else()
    target_compile_options(nvdb_probe PRIVATE -Wall -Wextra -Wpedantic)
endif()
