cmake_minimum_required(VERSION 3.16)
project(cg_cvgen LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CG_CVGEN_CORE_SOURCES
    src/Config.cpp
    src/FileUtils.cpp
    src/TextUtils.cpp
)

add_executable(cg_cvgen
    ${CG_CVGEN_CORE_SOURCES}
    src/Runner.cpp
    src/main.cpp
)

target_include_directories(cg_cvgen PRIVATE include)

if (MSVC)
    target_compile_options(cg_cvgen PRIVATE /W4)
else()
    target_compile_options(cg_cvgen PRIVATE -Wall -Wextra -Wpedantic)
endif()

# Dependency-free unit tests for the pure-logic modules (TextUtils, Config).
# Does not link Runner.cpp, since Runner's only job is to shell out to the
# Python backend - there is nothing to unit test there without actually
# running it. Needs nothing beyond the standard library: no Python, torch,
# or mlcolvar required.
enable_testing()
add_executable(cg_cvgen_tests
    ${CG_CVGEN_CORE_SOURCES}
    tests/test_main.cpp
    tests/test_text_utils.cpp
    tests/test_config.cpp
)
target_include_directories(cg_cvgen_tests PRIVATE include tests)
add_test(NAME cg_cvgen_tests COMMAND cg_cvgen_tests)
