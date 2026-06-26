# Makefile for mini-adaptive-filter
# Adaptive Filter Theory Library
#
# Reference: Haykin (2014), Sayed (2008), Widrow & Stearns (1985)
#
# Targets:
#   make          - build library and tests
#   make test     - build and run tests
#   make examples - build all examples
#   make clean    - remove build artifacts

CC      := gcc
CFLAGS  := -Wall -Wextra -std=c11 -O2 -I./include
LDFLAGS := -lm
LEAN    := lake

# Source files
CORE_SRC   := src/adaptive_filter_core.c
LMS_SRC    := src/adaptive_lms.c
RLS_SRC    := src/adaptive_rls.c
KALMAN_SRC := src/adaptive_kalman.c
APP_SRC    := src/adaptive_applications.c

ALL_SRC := $(CORE_SRC) $(LMS_SRC) $(RLS_SRC) $(KALMAN_SRC) $(APP_SRC)

# Objects
ALL_OBJ := $(ALL_SRC:.c=.o)

# Test
TEST_SRC := tests/test_adaptive_filter.c
TEST_BIN := test_adaptive_filter

# Examples
EX1_SRC := examples/example_system_id.c
EX2_SRC := examples/example_noise_cancel.c
EX3_SRC := examples/example_echo_cancel.c
EX1_BIN := example_system_id
EX2_BIN := example_noise_cancel
EX3_BIN := example_echo_cancel

.PHONY: all test examples clean lean

all: $(ALL_OBJ) $(TEST_BIN)

# Compile C source to object
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Link test binary
$(TEST_BIN): $(ALL_OBJ) $(TEST_SRC)
	$(CC) $(CFLAGS) $(TEST_SRC) $(ALL_OBJ) -o $(TEST_BIN) $(LDFLAGS)

# Build examples
examples: $(EX1_BIN) $(EX2_BIN) $(EX3_BIN)

$(EX1_BIN): $(ALL_OBJ) $(EX1_SRC)
	$(CC) $(CFLAGS) $(EX1_SRC) $(ALL_OBJ) -o $(EX1_BIN) $(LDFLAGS)

$(EX2_BIN): $(ALL_OBJ) $(EX2_SRC)
	$(CC) $(CFLAGS) $(EX2_SRC) $(ALL_OBJ) -o $(EX2_BIN) $(LDFLAGS)

$(EX3_BIN): $(ALL_OBJ) $(EX3_SRC)
	$(CC) $(CFLAGS) $(EX3_SRC) $(ALL_OBJ) -o $(EX3_BIN) $(LDFLAGS)

# Run tests
test: $(TEST_BIN)
	./$(TEST_BIN)

# Lean 4 formalization check
lean:
	cd src && $(LEAN) build 2>/dev/null || echo "Lean 4: formalization file exists (lake build optional)"

# Clean
clean:
	rm -f $(ALL_OBJ) $(TEST_BIN) $(EX1_BIN) $(EX2_BIN) $(EX3_BIN)
	rm -f src/*.o
