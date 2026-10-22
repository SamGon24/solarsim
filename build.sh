#!/usr/bin/env sh
# =============================================================================
#  SolarSim - one-shot build script for Linux / macOS / BSD / MSYS2
#
#    ./build.sh            Release build
#    ./build.sh Debug      Debug build (-O0 -g)
#    ./build.sh test       Release build + run the unit tests
#
#  Uses CMake when available, otherwise GNU Make, otherwise a single direct
#  compiler invocation - so it works with nothing but a C++17 compiler.
# =============================================================================
set -e
cd "$(dirname "$0")"

MODE="${1:-Release}"
RUN_TESTS=0
if [ "$MODE" = "test" ]; then
  MODE=Release
  RUN_TESTS=1
fi

# Pick a compiler: $CXX, then c++, g++, clang++
if [ -z "$CXX" ]; then
  for c in c++ g++ clang++; do
    if command -v "$c" >/dev/null 2>&1; then CXX="$c"; break; fi
  done
fi
if [ -z "$CXX" ]; then
  echo "error: no C++ compiler found. Install g++ (>= 8) or clang++ (>= 7)." >&2
  exit 1
fi

JOBS=2
if command -v nproc >/dev/null 2>&1; then JOBS=$(nproc); elif command -v sysctl >/dev/null 2>&1; then JOBS=$(sysctl -n hw.ncpu 2>/dev/null || echo 2); fi

if command -v cmake >/dev/null 2>&1; then
  echo "==> CMake build ($MODE, $JOBS jobs, compiler $CXX)"
  CXX="$CXX" cmake -S . -B build -DCMAKE_BUILD_TYPE="$MODE"
  cmake --build build -j "$JOBS"
  BIN=build/solarsim
  [ "$RUN_TESTS" = 1 ] && ./build/solarsim_tests
elif command -v make >/dev/null 2>&1; then
  echo "==> Make build ($JOBS jobs, compiler $CXX)"
  if [ "$MODE" = "Debug" ]; then
    make -j "$JOBS" CXX="$CXX" debug
    BIN=./solarsim_debug
  else
    make -j "$JOBS" CXX="$CXX"
    BIN=./solarsim
  fi
  [ "$RUN_TESTS" = 1 ] && make CXX="$CXX" tests
else
  echo "==> Direct compiler build (no CMake/Make found)"
  OPT="-O2"
  [ "$MODE" = "Debug" ] && OPT="-O0 -g"
  # shellcheck disable=SC2086
  "$CXX" -std=c++17 $OPT -Isrc src/core/*.cpp src/astro/*.cpp src/physics/*.cpp src/analysis/*.cpp \
      src/render/*.cpp src/io/*.cpp src/platform/*.cpp src/app/*.cpp src/main.cpp -o solarsim -pthread
  BIN=./solarsim
fi

echo
echo "Build complete: $BIN"
echo "Try:  $BIN selftest      $BIN validate      $BIN view      $BIN help"
