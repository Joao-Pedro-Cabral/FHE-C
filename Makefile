# ----------------------------------------------------------------------
# Toolchain
# ----------------------------------------------------------------------

CC := gcc

# ARM cross compiler for Ubuntu 22.04
ARM_CC ?= aarch64-linux-gnu-gcc


# ----------------------------------------------------------------------
# Architecture
# ----------------------------------------------------------------------

ARCH ?= x86

ifeq ($(ARCH),x86)

CC := gcc

ARCH_CFLAGS := \
					-O3 -fopenmp \
					-mavx2 -mfma -march=native

BUILD_DIR := build

else ifeq ($(ARCH),arm)

CC := $(ARM_CC)

ARCH_CFLAGS := \
					-O3 -fopenmp \
					-march=armv8.2-a+simd+crypto

BUILD_DIR := build-arm

else

$(error Unsupported ARCH='$(ARCH)'. Use ARCH=x86 or ARCH=arm)

endif

# ----------------------------------------------------------------------
# Optimization mode
# ----------------------------------------------------------------------

OPT ?= 0

ifeq ($(OPT),1)

		ifeq ($(ARCH),x86)
				OPT_CFLAGS := -D__AVX2_MODE__

		else ifeq ($(ARCH),arm)
				OPT_CFLAGS := -D__NEON_MODE__

		endif

else

		OPT_CFLAGS :=

endif

# ----------------------------------------------------------------------
# Configuration parameters
# ----------------------------------------------------------------------

N2			?= 16
L			 ?= 4
SET_N	 ?= 0
N			 ?= 512
N_TESTS ?= 1000

PARAM_CFLAGS := \
					-DN2=$(N2) \
					-DL=$(L) \
					-DN=$(N) \
					-DSET_N=$(SET_N) \
					-DN_TESTS=$(N_TESTS)

# ----------------------------------------------------------------------
# Common compiler flags
# ----------------------------------------------------------------------

COMMON_CFLAGS := \
					-Wall -Wextra \
					-Wshadow -Wconversion \
					-Wdouble-promotion \
					-Wformat=2 \
					-Isrc

CFLAGS := $(COMMON_CFLAGS) $(ARCH_CFLAGS) $(OPT_CFLAGS) $(PARAM_CFLAGS)

LDLIBS := -lm

# ----------------------------------------------------------------------
# Directories
# ----------------------------------------------------------------------

SRC_DIR := src
TEST_DIR := test


# ----------------------------------------------------------------------
# Sources
# ----------------------------------------------------------------------

SRC_C := $(wildcard $(SRC_DIR)/*.c)

SRC_S := $(wildcard $(SRC_DIR)/*.S)


# ----------------------------------------------------------------------
# Test selection
# ----------------------------------------------------------------------

TEST ?= ckks_test

TARGET := $(BUILD_DIR)/$(TEST)

SRC_OBJ := $(patsubst %.c,$(BUILD_DIR)/%.o,$(SRC_C)) \
					 $(patsubst %.S,$(BUILD_DIR)/%.o,$(SRC_S))

TEST_OBJ := $(BUILD_DIR)/$(TEST).o


# ----------------------------------------------------------------------
# Phony targets
# ----------------------------------------------------------------------

.PHONY: all build run debug valgrind vtune cloud_scp cloud_ssh clean FORCE \
				openfhe openfhe_build openfhe_run openfhe_clean

all: run


# ----------------------------------------------------------------------
# Build
# ----------------------------------------------------------------------

build: $(TARGET)

$(TARGET): $(SRC_OBJ) $(TEST_OBJ)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $^ -o $@ $(LDLIBS)


$(BUILD_DIR)/%.o: %.c FORCE
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -c $< -o $@


$(BUILD_DIR)/%.o: %.S FORCE
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -c $< -o $@


$(BUILD_DIR)/%.o: $(TEST_DIR)/%.c FORCE
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -c $< -o $@

# ----------------------------------------------------------------------
# OpenMP runtime configuration
# ----------------------------------------------------------------------

OMP_NUM_THREADS ?= 4

export OMP_NUM_THREADS

# ----------------------------------------------------------------------
# Run
# ----------------------------------------------------------------------

run: CFLAGS += -D__RUN_MODE__
run: $(TARGET)
	./$(TARGET)


# ----------------------------------------------------------------------
# Debug
# ----------------------------------------------------------------------

debug: CFLAGS += -D__DEBUG_MODE__ -g
debug: $(TARGET)
	./$(TARGET)


# ----------------------------------------------------------------------
# Valgrind
# ----------------------------------------------------------------------

VALGRIND_OPTIONS := \
	-s \
	--leak-check=full \
	--track-origins=yes \
	--show-leak-kinds=all

valgrind: CFLAGS += -D__DEBUG_MODE__ -g
valgrind: $(TARGET)
	valgrind $(VALGRIND_OPTIONS) ./$(TARGET)


# ----------------------------------------------------------------------
# Intel VTune
# ----------------------------------------------------------------------

vtune: CFLAGS += -D__VTUNE_MODE__ -g
vtune: $(TARGET)
	@result_dir=$$(mktemp -d /tmp/vtune-XXXXXX); \
	echo "VTune result directory: $$result_dir"; \
	vtune \
		-collect hotspots \
		-result-dir "$$result_dir" \
		-- \
		"$(abspath $(TARGET))"


# ----------------------------------------------------------------------
# gem5
# ----------------------------------------------------------------------

GEM5 ?= $(GEM5_PATH)/build/ARM/gem5.opt
GEM5_CONFIG ?= gem5/gem5_core.py

gem5: ARCH=arm
gem5: CFLAGS += -D__ARM_NEON__
gem5: $(TARGET)
	$(GEM5) -d $(BUILD_DIR) $(GEM5_CONFIG) $(abspath $(TARGET))

# ----------------------------------------------------------------------
# Cloud
# ----------------------------------------------------------------------

cloud_ssh:
	gcloud compute ssh axion-arm-3 --zone=northamerica-northeast2-b --ssh-key-file=~/.ssh/gcloud_key

cloud_scp:
	gcloud compute scp --recurse Makefile	axion-arm-3:~/FHE-C/ --zone=northamerica-northeast2-b --ssh-key-file=~/.ssh/gcloud_key
	gcloud compute scp --recurse test	axion-arm-3:~/FHE-C/ --zone=northamerica-northeast2-b --ssh-key-file=~/.ssh/gcloud_key
	gcloud compute scp --recurse src	axion-arm-3:~/FHE-C/ --zone=northamerica-northeast2-b --ssh-key-file=~/.ssh/gcloud_key


# ----------------------------------------------------------------------
# OpenFHE Toolchain
# ----------------------------------------------------------------------

CXX := g++

ARM_CXX ?= aarch64-linux-gnu-g++

OPENFHE_TEST := openfhe_perf

OPENFHE_SRC := test/$(OPENFHE_TEST).cpp

OPENFHE_BUILD_DIR := build-openfhe

ifeq ($(ARCH),x86)
		OPENFHE_CXX := g++
else ifeq ($(ARCH),arm)
		OPENFHE_CXX := $(ARM_CXX)
endif

OPENFHE_CXXFLAGS := \
		-O3 \
		-std=c++17 \
		-DCKKS_N2=$(N2) \
		-DCKKS_L=$(L) \
		-DCKKS_N=$(N) \
		-DCKKS_N_TESTS=$(N_TESTS)

# ----------------------------------------------------------------------
# OpenFHE Libraries
# ----------------------------------------------------------------------

OPENFHE_PREFIX ?= /usr/local

OPENFHE_INCLUDE := $(OPENFHE_PREFIX)/include/openfhe

OPENFHE_LIBDIR := $(OPENFHE_PREFIX)/lib

OPENFHE_LDLIBS := \
		-lOPENFHEpke \
		-lOPENFHEcore \
		-lOPENFHEbinfhe \
		-Wl,-rpath,$(OPENFHE_LIBDIR)

OPENFHE_CXXFLAGS += \
		-I$(OPENFHE_INCLUDE) \
		-I$(OPENFHE_INCLUDE)/core \
		-I$(OPENFHE_INCLUDE)/pke \
		-I$(OPENFHE_INCLUDE)/binfhe

OPENFHE_LDFLAGS := \
		-L$(OPENFHE_LIBDIR)

# ----------------------------------------------------------------------
# OpenFHE CKKS Benchmark
# ----------------------------------------------------------------------

OPENFHE_TARGET := $(OPENFHE_BUILD_DIR)/$(OPENFHE_TEST)

openfhe: openfhe_run

openfhe_build: $(OPENFHE_TARGET)

$(OPENFHE_TARGET): $(OPENFHE_SRC) FORCE
	@mkdir -p $(OPENFHE_BUILD_DIR)
	$(OPENFHE_CXX) \
			$(OPENFHE_CXXFLAGS) \
			$< \
			-o $@ \
			$(OPENFHE_LDFLAGS) \
			$(OPENFHE_LDLIBS)

openfhe_run: $(OPENFHE_TARGET)
	./$(OPENFHE_TARGET)

openfhe_clean:
	rm -rf $(OPENFHE_BUILD_DIR)

# ----------------------------------------------------------------------
# Clean
# ----------------------------------------------------------------------

clean:
	rm -rf build build-arm build-openfhe

# FORCE causes objects to be recompiled for each mode.
FORCE:
