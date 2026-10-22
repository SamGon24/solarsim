# =============================================================================
#  SolarSim - CMake build
#
#    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
#    cmake --build build -j
#    ./build/solarsim            (Windows: build\Release\solarsim.exe with MSVC)
#
#  Options:
#    -DSOLARSIM_BUILD_TESTS=ON/OFF   build the unit-test executable (default ON)
#    -DSOLARSIM_NATIVE=ON            optimise for the build machine (-march=native)
#    -DSOLARSIM_WARNINGS_AS_ERRORS=ON
# =============================================================================
cmake_minimum_required(VERSION 3.10)
project(SolarSim VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel)
endif()

option(SOLARSIM_BUILD_TESTS "Build the unit tests" ON)
option(SOLARSIM_NATIVE "Optimise for the host CPU (-march=native)" OFF)
option(SOLARSIM_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)

find_package(Threads REQUIRED)

# ---- library with everything except main() --------------------------------------
file(GLOB SOLARSIM_LIB_SOURCES CONFIGURE_DEPENDS
  ${CMAKE_CURRENT_SOURCE_DIR}/src/core/*.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/astro/*.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/physics/*.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/analysis/*.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/render/*.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/io/*.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/platform/*.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/app/*.cpp)

add_library(solarsim_core STATIC ${SOLARSIM_LIB_SOURCES})
target_include_directories(solarsim_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(solarsim_core PUBLIC Threads::Threads)

if(MSVC)
  target_compile_options(solarsim_core PUBLIC /W4 /utf-8 /permissive- /EHsc /bigobj)
  target_compile_definitions(solarsim_core PUBLIC _CRT_SECURE_NO_WARNINGS NOMINMAX WIN32_LEAN_AND_MEAN _USE_MATH_DEFINES)
  if(SOLARSIM_WARNINGS_AS_ERRORS)
    target_compile_options(solarsim_core PUBLIC /WX)
  endif()
else()
  target_compile_options(solarsim_core PUBLIC -Wall -Wextra -Wpedantic)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(solarsim_core PUBLIC -Wshadow=local)
  endif()
  if(SOLARSIM_WARNINGS_AS_ERRORS)
    target_compile_options(solarsim_core PUBLIC -Werror)
  endif()
  if(SOLARSIM_NATIVE)
    target_compile_options(solarsim_core PUBLIC -march=native)
  endif()
endif()

# MinGW: C99-conformant printf (%zu, %lld) and a static runtime so the .exe
# runs without extra DLLs.
if(MINGW)
  target_compile_definitions(solarsim_core PUBLIC __USE_MINGW_ANSI_STDIO=1)
  target_link_libraries(solarsim_core PUBLIC -static -static-libgcc -static-libstdc++)
endif()

# ---- executable -----------------------------------------------------------------------
add_executable(solarsim src/main.cpp)
target_link_libraries(solarsim PRIVATE solarsim_core)

# ---- tests --------------------------------------------------------------------------------
if(SOLARSIM_BUILD_TESTS)
  enable_testing()
  file(GLOB SOLARSIM_TEST_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp)
  add_executable(solarsim_tests ${SOLARSIM_TEST_SOURCES})
  target_link_libraries(solarsim_tests PRIVATE solarsim_core)
  add_test(NAME unit_tests COMMAND solarsim_tests)
  add_test(NAME selftest COMMAND solarsim selftest)
endif()

install(TARGETS solarsim RUNTIME DESTINATION bin)
install(DIRECTORY data/ DESTINATION share/solarsim)
