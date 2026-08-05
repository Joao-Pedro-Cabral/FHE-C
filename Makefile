CC := gcc

CFLAGS := -O3 \
          -Wall -Wextra \
          -Wshadow -Wconversion \
          -Wundef -Wdouble-promotion \
          -Wformat=2 \
          -Isrc
LDLIBS = -lm

VALGRIND_OPTIONS := -s --leak-check=full --track-origins=yes --show-leak-kinds=all

SRC_DIR := src
TEST_DIR := test
BUILD_DIR := build

SRC := $(wildcard $(SRC_DIR)/*.c)

TEST ?=

ifeq ($(TEST),)
$(error Usage: make TEST=<test_name>)
endif

TARGET := $(BUILD_DIR)/$(TEST)

SRC_OBJ := $(patsubst $(SRC_DIR)/%.c,$(BUILD_DIR)/%.o,$(SRC))
TEST_OBJ := $(BUILD_DIR)/$(TEST).o

.PHONY: all run valgrind clean

all: $(TARGET)

$(TARGET): $(SRC_OBJ) $(TEST_OBJ)
	@mkdir -p $(BUILD_DIR)
	$(CC) $(CFLAGS) $^ -o $@ $(LDLIBS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(BUILD_DIR)
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
