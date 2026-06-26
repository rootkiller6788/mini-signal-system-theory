# Makefile — mini-fourier-analysis
# Builds library, tests, examples, demos, and benchmarks.

CC       = gcc
CFLAGS   = -std=c11 -Wall -Wextra -O2 -g -Iinclude
LDFLAGS  = -lm

# Source files for the library
LIB_SRCS = src/fourier_core.c \
           src/fourier_fft.c \
           src/fourier_window.c \
           src/fourier_spectrum.c \
           src/fourier_convolution.c \
           src/fourier_applications.c

LIB_OBJS = $(LIB_SRCS:.c=.o)

# Test
TEST_SRC   = tests/test_fourier.c
TEST_BIN   = tests/test_fourier

# Examples
EX_SRC1 = examples/example_spectrum_analyzer.c
EX_SRC2 = examples/example_radar_pulse.c
EX_SRC3 = examples/example_vibration_analysis.c
EX_BIN1 = examples/example_spectrum_analyzer
EX_BIN2 = examples/example_radar_pulse
EX_BIN3 = examples/example_vibration_analysis

.PHONY: all clean test examples

all: test examples

# ─── Library objects ───

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# ─── Test ───

$(TEST_BIN): $(TEST_SRC) $(LIB_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

test: $(TEST_BIN)
	./$(TEST_BIN)

# ─── Examples ───

$(EX_BIN1): $(EX_SRC1) $(LIB_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(EX_BIN2): $(EX_SRC2) $(LIB_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(EX_BIN3): $(EX_SRC3) $(LIB_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

examples: $(EX_BIN1) $(EX_BIN2) $(EX_BIN3)

# ─── Run examples ───

run-example1: $(EX_BIN1)
	./$(EX_BIN1)

run-example2: $(EX_BIN2)
	./$(EX_BIN2)

run-example3: $(EX_BIN3)
	./$(EX_BIN3)

# ─── Clean ───

clean:
	rm -f src/*.o
	rm -f $(TEST_BIN)
	rm -f $(EX_BIN1) $(EX_BIN2) $(EX_BIN3)
