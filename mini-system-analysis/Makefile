# Makefile for mini-system-analysis
CC = gcc
CFLAGS = -std=c11 -Wall -Wextra -O1
LDFLAGS = -lm
INCLUDE = -Iinclude

SRC_DIR = src
TEST_DIR = tests
EXAMPLE_DIR = examples
DEMO_DIR = demos
BENCH_DIR = benches
BUILD_DIR = build

TESTS = test_system test_stability test_kalman
EXAMPLES = ex1_step_response ex2_stability ex3_state_space
DEMOS = demo_bode_analysis
BENCHES = bench_convolution

.PHONY: all test demo bench clean

all: $(addprefix $(BUILD_DIR)/, $(EXAMPLES) $(DEMOS))
	@echo "All examples and demos built."

test: $(addprefix $(BUILD_DIR)/, $(TESTS))
	@echo "=== Running tests ==="
	@$(BUILD_DIR)/test_system || (echo "test_system FAILED" && exit 1)
	@$(BUILD_DIR)/test_stability || (echo "test_stability FAILED" && exit 1)
	@cd $(BUILD_DIR) && ./test_kalman || ./test_kalman.exe || (echo "test_kalman FAILED" && exit 1)
	@echo "=== All tests passed ==="

demo: $(addprefix $(BUILD_DIR)/, $(DEMOS))
	@echo "=== Running demos ==="
	@$(BUILD_DIR)/demo_bode_analysis

bench: $(addprefix $(BUILD_DIR)/, $(BENCHES))
	@echo "=== Running benchmarks ==="
	@$(BUILD_DIR)/bench_convolution

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/test_system: $(SRC_DIR)/system_defs.c $(SRC_DIR)/convolution.c $(TEST_DIR)/test_system.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) $(INCLUDE) $(SRC_DIR)/system_defs.c $(SRC_DIR)/convolution.c $(TEST_DIR)/test_system.c -o $@ $(LDFLAGS)

$(BUILD_DIR)/test_stability: $(SRC_DIR)/system_defs.c $(SRC_DIR)/transfer_function.c $(SRC_DIR)/state_space.c $(SRC_DIR)/stability.c $(SRC_DIR)/frequency_response.c $(TEST_DIR)/test_stability.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) $(INCLUDE) $(SRC_DIR)/system_defs.c $(SRC_DIR)/transfer_function.c $(SRC_DIR)/state_space.c $(SRC_DIR)/stability.c $(SRC_DIR)/frequency_response.c $(TEST_DIR)/test_stability.c -o $@ $(LDFLAGS)

$(BUILD_DIR)/ex1_step_response: $(SRC_DIR)/system_defs.c $(SRC_DIR)/convolution.c $(SRC_DIR)/transfer_function.c $(SRC_DIR)/frequency_response.c $(EXAMPLE_DIR)/ex1_step_response.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) $(INCLUDE) $(SRC_DIR)/system_defs.c $(SRC_DIR)/convolution.c $(SRC_DIR)/transfer_function.c $(SRC_DIR)/frequency_response.c $(EXAMPLE_DIR)/ex1_step_response.c -o $@ $(LDFLAGS)

$(BUILD_DIR)/ex2_stability: $(SRC_DIR)/system_defs.c $(SRC_DIR)/transfer_function.c $(SRC_DIR)/stability.c $(SRC_DIR)/state_space.c $(SRC_DIR)/frequency_response.c $(EXAMPLE_DIR)/ex2_stability.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) $(INCLUDE) $(SRC_DIR)/system_defs.c $(SRC_DIR)/transfer_function.c $(SRC_DIR)/stability.c $(SRC_DIR)/state_space.c $(SRC_DIR)/frequency_response.c $(EXAMPLE_DIR)/ex2_stability.c -o $@ $(LDFLAGS)

$(BUILD_DIR)/ex3_state_space: $(SRC_DIR)/system_defs.c $(SRC_DIR)/state_space.c $(SRC_DIR)/stability.c $(SRC_DIR)/transfer_function.c $(SRC_DIR)/frequency_response.c $(EXAMPLE_DIR)/ex3_state_space.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) $(INCLUDE) $(SRC_DIR)/system_defs.c $(SRC_DIR)/state_space.c $(SRC_DIR)/stability.c $(SRC_DIR)/transfer_function.c $(SRC_DIR)/frequency_response.c $(EXAMPLE_DIR)/ex3_state_space.c -o $@ $(LDFLAGS)

$(BUILD_DIR)/demo_bode_analysis: $(SRC_DIR)/system_defs.c $(SRC_DIR)/transfer_function.c $(SRC_DIR)/frequency_response.c $(DEMO_DIR)/demo_bode_analysis.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) $(INCLUDE) $(SRC_DIR)/system_defs.c $(SRC_DIR)/transfer_function.c $(SRC_DIR)/frequency_response.c $(DEMO_DIR)/demo_bode_analysis.c -o $@ $(LDFLAGS)

$(BUILD_DIR)/bench_convolution: $(SRC_DIR)/system_defs.c $(SRC_DIR)/convolution.c $(SRC_DIR)/transfer_function.c $(SRC_DIR)/state_space.c $(SRC_DIR)/stability.c $(SRC_DIR)/frequency_response.c $(SRC_DIR)/kalman_filter.c $(BENCH_DIR)/bench_convolution.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) $(INCLUDE) $(SRC_DIR)/system_defs.c $(SRC_DIR)/convolution.c $(SRC_DIR)/transfer_function.c $(SRC_DIR)/state_space.c $(SRC_DIR)/stability.c $(SRC_DIR)/frequency_response.c $(SRC_DIR)/kalman_filter.c $(BENCH_DIR)/bench_convolution.c -o $@ $(LDFLAGS)

$(BUILD_DIR)/test_kalman: $(SRC_DIR)/system_defs.c $(SRC_DIR)/kalman_filter.c $(SRC_DIR)/state_space.c $(TEST_DIR)/test_kalman.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) $(INCLUDE) $(SRC_DIR)/system_defs.c $(SRC_DIR)/kalman_filter.c $(SRC_DIR)/state_space.c $(TEST_DIR)/test_kalman.c -o $@ $(LDFLAGS)

clean:
	rm -rf $(BUILD_DIR)
