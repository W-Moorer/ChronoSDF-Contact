cmake_minimum_required(VERSION 3.18)
cmake_policy(SET CMP0091 NEW)

project(sdf_contact_benchmarks)

find_package(Chrono CONFIG)

if(NOT Chrono_FOUND)
  message("Could not find Chrono core package")
  return()
endif()

set(BENCHMARK_SOURCES
    src/main.cpp
)

add_executable(sdf_contact_benchmarks ${BENCHMARK_SOURCES})
target_compile_features(sdf_contact_benchmarks PRIVATE cxx_std_17)

target_compile_definitions(sdf_contact_benchmarks PRIVATE "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\"")
target_compile_definitions(sdf_contact_benchmarks PRIVATE "BENCHMARK_TOOL_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}\"")

if(MSVC)
    set_target_properties(sdf_contact_benchmarks PROPERTIES MSVC_RUNTIME_LIBRARY ${CHRONO_MSVC_RUNTIME_LIBRARY})
endif()

target_link_libraries(sdf_contact_benchmarks PRIVATE ${CHRONO_TARGETS})

add_CHRONO_DLLS_copy_command()

if(TARGET COPY_CHRONO_DLLS)
    add_dependencies(sdf_contact_benchmarks COPY_CHRONO_DLLS)
endif()

if(WIN32 AND TARGET Chrono::Chrono_core)
    add_custom_target(sync_benchmark_chrono_core
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                $<TARGET_FILE:Chrono::Chrono_core>
                ${PROJECT_BINARY_DIR}/$<CONFIG>
        VERBATIM
    )
    add_dependencies(sdf_contact_benchmarks sync_benchmark_chrono_core)
endif()
