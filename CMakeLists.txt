cmake_minimum_required(VERSION 3.15)
project(LunarMasconsOrbit C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -funroll-loops")
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-march=native" HAS_MARCH_NATIVE)
if(HAS_MARCH_NATIVE)
  add_compile_options(-march=native)
endif()

add_executable(lunarmc src/main.cpp)
# Optional native HDF5 for compressed scientific output; the executable still
# builds without it for CSV-only simulations, but checkpoints require HDF5.
find_package(HDF5 COMPONENTS C QUIET)
if(HDF5_FOUND)
  target_include_directories(lunarmc PRIVATE ${HDF5_INCLUDE_DIRS})
  target_link_libraries(lunarmc PRIVATE ${HDF5_LIBRARIES})
  target_compile_definitions(lunarmc PRIVATE MASCONS_HAS_HDF5=1)
  message(STATUS "HDF5 output and rolling HDF5 checkpoints enabled")
else()
  message(STATUS "HDF5 not found: CSV-only output remains available; restart requires HDF5")
endif()

# Apple clang ships without OpenMP; point CMake at Homebrew libomp if present.
if(APPLE AND NOT OpenMP_CXX_FOUND)
  execute_process(COMMAND brew --prefix libomp
                  OUTPUT_VARIABLE LIBOMP OUTPUT_STRIP_TRAILING_WHITESPACE
                  ERROR_QUIET RESULT_VARIABLE BREW_RC)
  if(BREW_RC EQUAL 0 AND EXISTS "${LIBOMP}/lib/libomp.dylib")
    set(OpenMP_CXX_FLAGS "-Xclang -fopenmp -I${LIBOMP}/include")
    set(OpenMP_CXX_LIB_NAMES omp)
    set(OpenMP_omp_LIBRARY "${LIBOMP}/lib/libomp.dylib")
  endif()
endif()

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  target_link_libraries(lunarmc PRIVATE OpenMP::OpenMP_CXX)
else()
  message(WARNING "OpenMP not found - Monte Carlo will run single-threaded")
endif()

include(CTest)
if(BUILD_TESTING)
  add_executable(mascons_core_tests tests/core_tests.cpp)
  target_include_directories(mascons_core_tests PRIVATE src)
  add_test(NAME core_physics_and_integrators COMMAND mascons_core_tests)
endif()
if(BUILD_TESTING)
  add_executable(gravity_probe tests/gravity_probe.cpp)
  target_include_directories(gravity_probe PRIVATE src)
  find_package(Python3 COMPONENTS Interpreter QUIET)
  if(Python3_Interpreter_FOUND)
    execute_process(COMMAND "${Python3_EXECUTABLE}" -c "import scipy" RESULT_VARIABLE SCIPY_RC OUTPUT_QUIET ERROR_QUIET)
    if(SCIPY_RC EQUAL 0)
      add_test(NAME independent_scipy_gravity
        COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/tests/validate_scipy.py"
                "$<TARGET_FILE:gravity_probe>" "${CMAKE_SOURCE_DIR}/data/gravity.bin")
    endif()
  endif()
endif()
if(BUILD_TESTING AND Python3_Interpreter_FOUND)
  add_test(NAME lunar_impact_regression
    COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/tests/test_impact_regression.py"
            "$<TARGET_FILE:lunarmc>" "${CMAKE_SOURCE_DIR}/data")
endif()

if(BUILD_TESTING)
  add_executable(third_body_tests tests/test_third_bodies.cpp)
  target_include_directories(third_body_tests PRIVATE src)
  add_test(NAME third_body_tides_and_ephemeris COMMAND third_body_tests)
endif()

if(BUILD_TESTING AND Python3_Interpreter_FOUND)
  add_test(NAME lunar_third_body_end_to_end
    COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/tests/test_third_body_regression.py"
            "$<TARGET_FILE:lunarmc>" "${CMAKE_SOURCE_DIR}/data")
endif()

if(BUILD_TESTING AND Python3_Interpreter_FOUND)
  add_test(NAME kozai_lidov_master_switch
    COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/tests/test_kl_switch.py"
            "$<TARGET_FILE:lunarmc>" "${CMAKE_SOURCE_DIR}/data")
endif()

if(BUILD_TESTING AND HDF5_FOUND AND Python3_Interpreter_FOUND)
  execute_process(COMMAND "${Python3_EXECUTABLE}" -c "import h5py" RESULT_VARIABLE H5PY_RC OUTPUT_QUIET ERROR_QUIET)
  if(H5PY_RC EQUAL 0)
    add_test(NAME checkpoint_restart_and_interruption
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/tests/test_checkpoint_restart.py"
              "$<TARGET_FILE:lunarmc>" "${CMAKE_SOURCE_DIR}/data")
    add_test(NAME native_hdf5_output_and_restart
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/tests/test_hdf5_output.py"
              "$<TARGET_FILE:lunarmc>" "${CMAKE_SOURCE_DIR}/data")
  endif()
endif()

if(BUILD_TESTING AND Python3_Interpreter_FOUND)
  add_test(NAME ensemble_stop_all_impacted
    COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/tests/test_stop_all_impacted.py"
            "$<TARGET_FILE:lunarmc>" "${CMAKE_SOURCE_DIR}/data")
endif()

if(BUILD_TESTING AND Python3_Interpreter_FOUND)
  add_test(NAME maximum_available_harmonic_degree
    COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/tests/test_max_degree.py"
            "$<TARGET_FILE:lunarmc>" "${CMAKE_SOURCE_DIR}/data")
endif()

if(BUILD_TESTING AND Python3_Interpreter_FOUND)
  add_test(NAME adaptive_terrain_event_regression
    COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/tests/test_adaptive_event.py"
            "$<TARGET_FILE:lunarmc>" "${CMAKE_SOURCE_DIR}/data")
endif()

if(BUILD_TESTING AND HDF5_FOUND AND Python3_Interpreter_FOUND AND H5PY_RC EQUAL 0)
  add_test(NAME adaptive_high_altitude_checkpoint_restart
    COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/tests/test_adaptive_restart.py"
      "$<TARGET_FILE:lunarmc>" "${CMAKE_SOURCE_DIR}/data")
endif()

if(BUILD_TESTING AND HDF5_FOUND AND Python3_Interpreter_FOUND AND H5PY_RC EQUAL 0)
  add_test(NAME live_progress_and_resume
    COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/tests/test_live_progress.py"
            "$<TARGET_FILE:lunarmc>" "${CMAKE_SOURCE_DIR}/data")
endif()

if(BUILD_TESTING)
  add_executable(progress_schedule_tests tests/test_progress_schedule.cpp)
  target_include_directories(progress_schedule_tests PRIVATE src)
  add_test(NAME ten_minute_poll_thirty_minute_heartbeat COMMAND progress_schedule_tests)
endif()
