cmake_minimum_required(VERSION 3.20)

project(mesh_to_vdb_nvdb LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(OpenVDB CONFIG REQUIRED)

add_executable(mesh_to_vdb_nvdb
    src/main.cpp
    src/MeshIO.cpp
    ../../chrono/src/chrono_thirdparty/libstl/stlfile.c
)

target_include_directories(mesh_to_vdb_nvdb
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/../../chrono/src
)

target_compile_definitions(mesh_to_vdb_nvdb PRIVATE NANOVDB_USE_OPENVDB)

target_link_libraries(mesh_to_vdb_nvdb PRIVATE OpenVDB::openvdb)

if(MSVC)
    target_compile_options(mesh_to_vdb_nvdb PRIVATE /W4 /permissive- /bigobj)
else()
    target_compile_options(mesh_to_vdb_nvdb PRIVATE -Wall -Wextra -Wpedantic)
endif()
