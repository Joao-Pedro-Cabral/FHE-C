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
          -O3 \
          -mavx2 -mfma -march=native

BUILD_DIR := build

else ifeq ($(ARCH),arm)

CC := $(ARM_CC)

ARCH_CFLAGS := \
          -O3 \
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
# Common compiler flags
# ----------------------------------------------------------------------

COMMON_CFLAGS := \
          -Wall -Wextra \
          -Wshadow -Wconversion \
          -Wdouble-promotion \
          -Wformat=2 \
          -Isrc

CFLAGS := $(COMMON_CFLAGS) $(ARCH_CFLAGS) $(OPT_CFLAGS)

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

.PHONY: all build run debug valgrind vtune cloud_scp cloud_ssh clean FORCE

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
	gcloud compute ssh axion-arm --zone=us-central1-c

cloud_scp:
	gcloud compute scp --recurse Makefile  axion-arm:~/FHE-C/ --zone=us-central1-c
	gcloud compute scp --recurse test  axion-arm:~/FHE-C/ --zone=us-central1-c
	gcloud compute scp --recurse src  axion-arm:~/FHE-C/ --zone=us-central1-c

# ----------------------------------------------------------------------
# Clean
# ----------------------------------------------------------------------

clean:
	rm -rf build build-arm

# FORCE causes objects to be recompiled for each mode.
FORCE:
