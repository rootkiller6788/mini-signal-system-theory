# Makefile — mini-filter-theory
# Filter Theory: analog/digital filter design, analysis, and implementation
#
# Targets:
#   make          — build library and tests
#   make test     — run test suite
#   make examples — build all examples
#   make all      — build everything
#   make clean    — remove build artifacts

CC       = gcc
CFLAGS   = -std=c11 -Wall -Wextra -O2 -g
INCLUDES = -Iinclude
LDFLAGS  = -lm

SRC_DIR  = src
TST_DIR  = tests
EX_DIR   = examples
DM_DIR   = demos
BN_DIR   = benches
BLD_DIR  = build

# Source files
SRCS = $(SRC_DIR)/filter_analog.c \
       $(SRC_DIR)/filter_transfer.c \
       $(SRC_DIR)/filter_digital.c \
       $(SRC_DIR)/filter_design.c \
       $(SRC_DIR)/filter_analysis.c \
       $(SRC_DIR)/filter_math.c

OBJS = $(SRCS:$(SRC_DIR)/%.c=$(BLD_DIR)/%.o)

# Library
LIB = $(BLD_DIR)/libfilter.a

# Tests
TEST_SRC = $(TST_DIR)/test_filter.c
TEST_BIN = $(BLD_DIR)/test_filter

# Examples
EX1_SRC = $(EX_DIR)/ex1_butterworth_lpf.c
EX2_SRC = $(EX_DIR)/ex2_fir_design.c
EX3_SRC = $(EX_DIR)/ex3_audio_crossover.c
EX4_SRC = $(EX_DIR)/ex4_bessel_pulse.c
EX1_BIN = $(BLD_DIR)/ex1_butterworth_lpf
EX2_BIN = $(BLD_DIR)/ex2_fir_design
EX3_BIN = $(BLD_DIR)/ex3_audio_crossover
EX4_BIN = $(BLD_DIR)/ex4_bessel_pulse

# Demo and bench
DEMO_SRC = $(DM_DIR)/demo_response.c
DEMO_BIN = $(BLD_DIR)/demo_response
BENCH_SRC = $(BN_DIR)/bench_filter.c
BENCH_BIN = $(BLD_DIR)/bench_filter

.PHONY: all test examples demo bench clean

all: $(LIB) $(TEST_BIN) $(EX1_BIN) $(EX2_BIN) $(EX3_BIN) $(EX4_BIN)

# Build directory
$(BLD_DIR):
	mkdir -p $(BLD_DIR)

# Library
$(LIB): $(BLD_DIR) $(OBJS)
	ar rcs $@ $(OBJS)

# Compile source files
$(BLD_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# Test binary
$(TEST_BIN): $(LIB) $(TEST_SRC)
	$(CC) $(CFLAGS) $(INCLUDES) $(TEST_SRC) -L$(BLD_DIR) -lfilter $(LDFLAGS) -o $@

# Example binaries
$(EX1_BIN): $(LIB) $(EX1_SRC)
	$(CC) $(CFLAGS) $(INCLUDES) $(EX1_SRC) -L$(BLD_DIR) -lfilter $(LDFLAGS) -o $@

$(EX2_BIN): $(LIB) $(EX2_SRC)
	$(CC) $(CFLAGS) $(INCLUDES) $(EX2_SRC) -L$(BLD_DIR) -lfilter $(LDFLAGS) -o $@

$(EX3_BIN): $(LIB) $(EX3_SRC)
	$(CC) $(CFLAGS) $(INCLUDES) $(EX3_SRC) -L$(BLD_DIR) -lfilter $(LDFLAGS) -o $@

$(EX4_BIN): $(LIB) $(EX4_SRC)
	$(CC) $(CFLAGS) $(INCLUDES) $(EX4_SRC) -L$(BLD_DIR) -lfilter $(LDFLAGS) -o $@

# Demo binary
$(DEMO_BIN): $(LIB) $(DEMO_SRC)
	$(CC) $(CFLAGS) $(INCLUDES) $(DEMO_SRC) -L$(BLD_DIR) -lfilter $(LDFLAGS) -o $@

# Bench binary
$(BENCH_BIN): $(LIB) $(BENCH_SRC)
	$(CC) $(CFLAGS) $(INCLUDES) $(BENCH_SRC) -L$(BLD_DIR) -lfilter $(LDFLAGS) -o $@

# Run test
test: $(TEST_BIN)
	./$(TEST_BIN)

# Build examples
examples: $(EX1_BIN) $(EX2_BIN) $(EX3_BIN) $(EX4_BIN)

# Build demo
demo: $(DEMO_BIN)

# Build bench
bench: $(BENCH_BIN)

# Clean
clean:
	rm -rf $(BLD_DIR)
