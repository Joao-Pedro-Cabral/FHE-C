# ----------------------------------------------------------------------
# Toolchain
# ----------------------------------------------------------------------

CC := gcc

# ARM cross compiler for Ubuntu 22.04
ARM_CC ?= aarch64-linux-gnu-gcc


# ----------------------------------------------------------------------
# Common compiler flags
# ----------------------------------------------------------------------

COMMON_CFLAGS := \
          -Wall -Wextra \
          -Wshadow -Wconversion \
          -Wundef -Wdouble-promotion \
          -Wformat=2 \
          -Isrc

LDLIBS := -lm

# ----------------------------------------------------------------------
# Optional AES
# ----------------------------------------------------------------------

ifeq ($(AES),1)
COMMON_CFLAGS += -D__AES_MODE__
endif

# ----------------------------------------------------------------------
# Architectures
# ----------------------------------------------------------------------

# Native x86-64
NATIVE_CFLAGS := \
          -O3 \
          -mavx2 -mfma -march=native \
          -Isrc/blake3 

# ARMv8.2-A / Cortex-A55
#
# Cortex-A55 supports Advanced SIMD (NEON) and FP.
#
ARM_CFLAGS := \
          -static \
          -O3 \
          -mcpu=cortex-a55 \
          -march=armv8.2-a+simd+crypto

CFLAGS := $(COMMON_CFLAGS) $(NATIVE_CFLAGS)

ARM_CFLAGS := $(COMMON_CFLAGS) $(ARM_CFLAGS)

# ----------------------------------------------------------------------
# Directories
# ----------------------------------------------------------------------

SRC_DIR := src
BLAKE3_DIR := $(SRC_DIR)/blake3
TEST_DIR := test
BUILD_DIR := build
ARM_BUILD_DIR := build-arm

# ----------------------------------------------------------------------
# Sources
# ----------------------------------------------------------------------

SRC_C := $(wildcard $(SRC_DIR)/*.c) \
         $(wildcard $(BLAKE3_DIR)/*.c)

SRC_S := $(wildcard $(SRC_DIR)/*.S) \
         $(wildcard $(BLAKE3_DIR)/*.S)

ARM_SRC_C := $(wildcard $(SRC_DIR)/*.c)

ARM_SRC_S := $(wildcard $(SRC_DIR)/*.S)

# ----------------------------------------------------------------------
# Test selection
# ----------------------------------------------------------------------

TEST ?= ckks_test

TARGET := $(BUILD_DIR)/$(TEST)

# Native objects
SRC_OBJ := $(patsubst %.c,$(BUILD_DIR)/%.o,$(SRC_C)) \
           $(patsubst %.S,$(BUILD_DIR)/%.o,$(SRC_S))

TEST_OBJ := $(BUILD_DIR)/$(TEST).o

# ARM objects
ARM_SRC_OBJ := $(patsubst %.c,$(ARM_BUILD_DIR)/%.o,$(ARM_SRC_C)) \
               $(patsubst %.S,$(ARM_BUILD_DIR)/%.o,$(ARM_SRC_S))

ARM_TEST_OBJ := $(ARM_BUILD_DIR)/$(TEST).o

ARM_TARGET := $(ARM_BUILD_DIR)/$(TEST)

# ----------------------------------------------------------------------
# Phony targets
# ----------------------------------------------------------------------

.PHONY: all build run debug valgrind vtune gem5 clean FORCE

all: run

# ----------------------------------------------------------------------
# Normal native build
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
# Normal run
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
# ARM / gem5 build
# ----------------------------------------------------------------------

gem5-build: ARM_CFLAGS += -D__ARM_NEON__
gem5-build: $(ARM_TARGET)


$(ARM_TARGET): $(ARM_SRC_OBJ) $(ARM_TEST_OBJ)
	@mkdir -p $(dir $@)
	$(ARM_CC) $(ARM_CFLAGS) $^ -o $@ $(LDLIBS)


$(ARM_BUILD_DIR)/%.o: %.c FORCE
	@mkdir -p $(dir $@)
	$(ARM_CC) $(COMMON_CFLAGS) $(ARM_CFLAGS) -c $< -o $@


$(ARM_BUILD_DIR)/%.o: %.S FORCE
	@mkdir -p $(dir $@)
	$(ARM_CC) $(COMMON_CFLAGS) $(ARM_CFLAGS) -c $< -o $@


$(ARM_BUILD_DIR)/%.o: $(TEST_DIR)/%.c FORCE
	@mkdir -p $(dir $@)
	$(ARM_CC) $(COMMON_CFLAGS) $(ARM_CFLAGS) -c $< -o $@


# ----------------------------------------------------------------------
# gem5 execution
# ----------------------------------------------------------------------

GEM5 ?= $(GEM5_PATH)/build/ARM/gem5.opt
GEM5_CONFIG ?= gem5/gem5_core.py

gem5: gem5-build
	$(GEM5) -d $(ARM_BUILD_DIR) $(GEM5_CONFIG) $(abspath $(ARM_TARGET))


# ----------------------------------------------------------------------
# Clean
# ----------------------------------------------------------------------

clean:
	rm -rf $(BUILD_DIR) $(ARM_BUILD_DIR)

# FORCE causes objects to be recompiled for each mode.
FORCE:
