cmake_minimum_required(VERSION 3.16)
project(cg_md LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(cg_md
    src/Config.cpp
    src/Dssp.cpp
    src/FileUtils.cpp
    src/GroUtils.cpp
    src/IndexBuilder.cpp
    src/MdpWriter.cpp
    src/PlumedWriter.cpp
    src/SamplingMonitor.cpp
    src/Shell.cpp
    src/StringUtils.cpp
    src/TopologyEditor.cpp
    src/Workflow.cpp
    src/WorkflowAdaptive.cpp
    src/WorkflowPrepare.cpp
    src/WorkflowSimulation.cpp
    src/WorkflowMetadynamics.cpp
    src/main.cpp
)

target_include_directories(cg_md PRIVATE include)

target_compile_options(cg_md PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall -Wextra -Wpedantic>
)
