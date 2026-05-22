cmake_minimum_required(VERSION 3.16)
project(cg_cvgen LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(cg_cvgen
    src/main.cpp
    src/Config.cpp
    src/FileUtils.cpp
    src/Runner.cpp
    src/TextUtils.cpp
)

target_include_directories(cg_cvgen PRIVATE include)

if (MSVC)
    target_compile_options(cg_cvgen PRIVATE /W4)
else()
    target_compile_options(cg_cvgen PRIVATE -Wall -Wextra -Wpedantic)
endif()
