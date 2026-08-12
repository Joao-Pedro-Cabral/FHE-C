CC := gcc

CFLAGS := -O3 -mavx2 -mfma -march=native \
          -Wall -Wextra \
          -Wshadow -Wconversion \
          -Wundef -Wdouble-promotion \
          -Wformat=2 \
          -Isrc -Isrc/blake3 \
          $(CFLAGS_EXTRA)
LDLIBS = -lm -lcrypto

VALGRIND_OPTIONS := -s --leak-check=full --track-origins=yes --show-leak-kinds=all

SRC_DIR := src
BLAKE3_DIR := $(SRC_DIR)/blake3
TEST_DIR := test
BUILD_DIR := build

SRC_C := $(wildcard $(SRC_DIR)/*.c) \
         $(wildcard $(BLAKE3_DIR)/*.c)

SRC_S := $(wildcard $(SRC_DIR)/*.S) \
         $(wildcard $(BLAKE3_DIR)/*.S)

SRC := $(SRC_C) $(SRC_S)

TEST ?=

ifeq ($(TEST),)
$(error Usage: make TEST=<test_name>)
endif

TARGET := $(BUILD_DIR)/$(TEST)

SRC_OBJ := $(patsubst %.c,$(BUILD_DIR)/%.o,$(SRC_C)) \
           $(patsubst %.S,$(BUILD_DIR)/%.o,$(SRC_S))
TEST_OBJ := $(BUILD_DIR)/$(TEST).o

.PHONY: all run valgrind clean

all: $(TARGET)

$(TARGET): $(SRC_OBJ) $(TEST_OBJ)
	@mkdir -p $(BUILD_DIR)
	$(CC) $(CFLAGS) $^ -o $@ $(LDLIBS)

$(BUILD_DIR)/%.o: %.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: %.S
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: $(TEST_DIR)/%.c
	@mkdir -p $(BUILD_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

run: $(TARGET)
	./$(TARGET)

valgrind: CFLAGS += -g
valgrind: $(TARGET)
	valgrind $(VALGRIND_OPTIONS) ./$(TARGET)

clean:
	rm -rf $(BUILD_DIR)
