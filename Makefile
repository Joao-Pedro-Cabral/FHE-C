CC := gcc

CFLAGS := -O3 -mavx2 -mfma -march=native \
          -Wall -Wextra \
          -Wshadow -Wconversion \
          -Wundef -Wdouble-promotion \
          -Wformat=2 \
          -Isrc -Isrc/blake3

ifeq ($(AES),1)
CFLAGS += -D__AES_AVX2__
endif

LDLIBS := -lm

SRC_DIR := src
BLAKE3_DIR := $(SRC_DIR)/blake3
TEST_DIR := test
BUILD_DIR := build

SRC_C := $(wildcard $(SRC_DIR)/*.c) \
         $(wildcard $(BLAKE3_DIR)/*.c)

SRC_S := $(wildcard $(SRC_DIR)/*.S) \
         $(wildcard $(BLAKE3_DIR)/*.S)

TEST ?= ckks_test

TARGET := $(BUILD_DIR)/$(TEST)

SRC_OBJ := $(patsubst %.c,$(BUILD_DIR)/%.o,$(SRC_C)) \
           $(patsubst %.S,$(BUILD_DIR)/%.o,$(SRC_S))

TEST_OBJ := $(BUILD_DIR)/$(TEST).o

.PHONY: all build run debug valgrind vtune clean FORCE

# Default: compile and run in normal mode
all: run

# Compile only
build: $(TARGET)

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

VALGRIND_OPTIONS := -s --leak-check=full --track-origins=yes --show-leak-kinds=all

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
# Build
# ----------------------------------------------------------------------

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


# FORCE causes the objects to be recompiled for each mode.
FORCE:


clean:
	rm -rf $(BUILD_DIR)
