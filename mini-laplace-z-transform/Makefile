CC = gcc
CFLAGS = -Wall -Wextra -std=c11 -O2 -g
INCLUDES = -Iinclude
LDFLAGS = -lm

SRC_DIR = src
TEST_DIR = tests
EXAMPLE_DIR = examples
BUILD_DIR = build

SRCS = $(SRC_DIR)/laplace_z_transform_core.c $(SRC_DIR)/laplace_transform.c $(SRC_DIR)/z_transform.c $(SRC_DIR)/partial_fraction.c $(SRC_DIR)/transfer_function.c $(SRC_DIR)/stability_analysis.c $(SRC_DIR)/inverse_transform.c

OBJS = $(SRCS:$(SRC_DIR)/%.c=$(BUILD_DIR)/%.o)

TESTS = $(TEST_DIR)/test_core.c $(TEST_DIR)/test_laplace.c $(TEST_DIR)/test_z.c $(TEST_DIR)/test_stability.c

TEST_BINS = $(TESTS:$(TEST_DIR)/%.c=$(BUILD_DIR)/%)

EXAMPLES = $(EXAMPLE_DIR)/ex1_first_order.c $(EXAMPLE_DIR)/ex2_second_order.c $(EXAMPLE_DIR)/ex3_bilinear_filter.c

EXAMPLE_BINS = $(EXAMPLES:$(EXAMPLE_DIR)/%.c=$(BUILD_DIR)/%)

.PHONY: all clean test examples

all: test examples

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

$(BUILD_DIR)/test_%: $(TEST_DIR)/test_%.c $(OBJS) | $(BUILD_DIR)
	$(CC) $(CFLAGS) $(INCLUDES) $< $(OBJS) -o $@ $(LDFLAGS)

$(BUILD_DIR)/ex%: $(EXAMPLE_DIR)/ex%.c $(OBJS) | $(BUILD_DIR)
	$(CC) $(CFLAGS) $(INCLUDES) $< $(OBJS) -o $@ $(LDFLAGS)

test: $(TEST_BINS)
	@echo "=== Running Tests ==="
	@for t in $(TEST_BINS); do echo "--- $$t ---"; ./$$t || exit 1; done
	@echo "=== All tests passed ==="

examples: $(EXAMPLE_BINS)
	@echo "=== Running Examples ==="
	@for e in $(EXAMPLE_BINS); do echo "--- $$e ---"; ./$$e; done

clean:
	rm -rf $(BUILD_DIR)
