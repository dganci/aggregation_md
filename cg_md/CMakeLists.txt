cmake_minimum_required(VERSION 3.16)
project(cg_md LANGUAGES CXX)

# Default to an optimized build. build/CMakeCache.txt was configured with an
# empty CMAKE_BUILD_TYPE, which means no optimisation flags at all: cg_md
# itself mostly orchestrates external tools, but the per-chunk analysis
# (TimeSeriesStats, ColvarTable, SamplingMonitor, ClusterAnalysis) runs on
# every adaptive chunk of every run, and -O0 is not what you want there.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CG_MD_CORE_SOURCES
    src/ClusterAnalysis.cpp
    src/ColvarTable.cpp
    src/ColvarTruncate.cpp
    src/Config.cpp
    src/ConfigHelp.cpp
    src/ConfigIdentity.cpp
    src/ConfigOptions.cpp
    src/ConfigValidate.cpp
    src/FileUtils.cpp
    src/GroUtils.cpp
    src/GromacsDriver.cpp
    src/IndexBuilder.cpp
    src/MdpProduction.cpp
    src/MdpWriter.cpp
    src/PdbUtils.cpp
    src/PlumedDescriptors.cpp
    src/PlumedLog.cpp
    src/PlumedWriter.cpp
    src/PlumedWriterMetad.cpp
    src/Reporter.cpp
    src/RunDiagnostics.cpp
    src/RunLedger.cpp
    src/CvReadiness.cpp
    src/SamplingEstimators.cpp
    src/SamplingMonitor.cpp
    src/Shell.cpp
    src/StringUtils.cpp
    src/TimeSeriesStats.cpp
    src/TopologyEditor.cpp
)

add_executable(cg_md
    ${CG_MD_CORE_SOURCES}
    src/AdaptiveSampler.cpp
    src/AdaptiveSamplerFinalize.cpp
    src/MetadynamicsCvParams.cpp
    src/MetadynamicsBatch.cpp
    src/MetadynamicsFinalize.cpp
    src/MetadynamicsGuards.cpp
    src/MetadynamicsInternal.cpp
    src/MetadynamicsRunner.cpp
    src/Preparer.cpp
    src/PreparerChecks.cpp
    src/PreparerCoarseGrain.cpp
    src/PreparerPackmol.cpp
    src/RelaxRunner.cpp
    src/RelaxRunnerSelect.cpp
    src/SimulationRunner.cpp
    src/Workflow.cpp
    src/WorkflowChecks.cpp
    src/main.cpp
)

target_include_directories(cg_md PRIVATE include)

target_compile_options(cg_md PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall -Wextra -Wpedantic>
)

# Dependency-free unit tests for the pure-logic modules (StringUtils,
# GroUtils, TopologyEditor, SamplingMonitor, IndexBuilder parsing). These do
# not spawn gmx/martinize2/packmol/plumed, so `ctest` (or running the
# executable directly) works without any MD software installed.
enable_testing()
add_executable(cg_md_tests
    ${CG_MD_CORE_SOURCES}
    tests/test_main.cpp
    tests/test_string_utils.cpp
    tests/test_gro_utils.cpp
    tests/test_topology_editor.cpp
    tests/test_sampling_monitor.cpp
    tests/test_sampling_estimators.cpp
    tests/test_cv_readiness.cpp
    tests/test_index_builder.cpp
    tests/test_cluster_analysis.cpp
    tests/test_time_series_stats.cpp
    tests/test_run_diagnostics.cpp
    tests/test_colvar_table.cpp
    tests/test_plumed_writer.cpp
    tests/test_plumed_metad.cpp
    tests/test_plumed_log.cpp
    tests/test_mdp_writer.cpp
    tests/test_mdp_reference.cpp
    tests/test_run_ledger.cpp
    tests/test_config_identity.cpp
    tests/test_pdb_utils.cpp
    tests/test_file_utils.cpp
)
target_include_directories(cg_md_tests PRIVATE include tests)
add_test(NAME cg_md_tests COMMAND cg_md_tests)
