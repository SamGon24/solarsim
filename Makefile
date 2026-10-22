# =============================================================================
#  SolarSim - plain GNU Make build (no CMake required)
#
#    make              -> ./solarsim           (release, -O2)
#    make debug        -> ./solarsim_debug     (-O0 -g, assertions)
#    make tests        -> ./solarsim_tests and runs them
#    make clean
#
#  Override the compiler with e.g.  make CXX=clang++
# =============================================================================
CXX      ?= g++
CXXFLAGS ?= -O2
STD      := -std=c++17
WARN     := -Wall -Wextra -Wpedantic
INC      := -Isrc
LDLIBS   := -pthread

# MinGW / MSYS2 on Windows: C99 printf and a self-contained executable.
ifeq ($(OS),Windows_NT)
  CXXFLAGS += -D__USE_MINGW_ANSI_STDIO=1
  LDLIBS   += -static -static-libgcc -static-libstdc++
endif

SRC_DIRS := src/core src/astro src/physics src/analysis src/render src/io src/platform src/app
LIB_SRC  := $(foreach d,$(SRC_DIRS),$(wildcard $(d)/*.cpp))
TEST_SRC := $(wildcard tests/*.cpp)

BUILD    := build-make
LIB_OBJ  := $(patsubst %.cpp,$(BUILD)/release/%.o,$(LIB_SRC))
DBG_OBJ  := $(patsubst %.cpp,$(BUILD)/debug/%.o,$(LIB_SRC))
TEST_OBJ := $(patsubst %.cpp,$(BUILD)/release/%.o,$(TEST_SRC))

.PHONY: all debug tests clean run

all: solarsim

solarsim: $(LIB_OBJ) $(BUILD)/release/src/main.o
	$(CXX) $(STD) $(CXXFLAGS) -o $@ $^ $(LDLIBS)

solarsim_debug: $(DBG_OBJ) $(BUILD)/debug/src/main.o
	$(CXX) $(STD) -O0 -g -o $@ $^ $(LDLIBS)

debug: solarsim_debug

solarsim_tests: $(LIB_OBJ) $(TEST_OBJ)
	$(CXX) $(STD) $(CXXFLAGS) -o $@ $^ $(LDLIBS)

tests: solarsim_tests
	./solarsim_tests

$(BUILD)/release/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(STD) $(CXXFLAGS) $(WARN) $(INC) -MMD -MP -c $< -o $@

$(BUILD)/debug/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(STD) -O0 -g $(WARN) $(INC) -MMD -MP -c $< -o $@

run: solarsim
	./solarsim

clean:
	rm -rf $(BUILD) solarsim solarsim_debug solarsim_tests

-include $(LIB_OBJ:.o=.d) $(DBG_OBJ:.o=.d) $(TEST_OBJ:.o=.d)
