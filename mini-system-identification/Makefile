# Makefile for mini-system-identification
# System Identification module — L1-L9 knowledge coverage
#
# Reference: Ljung (1999), Söderström & Stoica (1989),
#            Van Overschee & De Moor (1996)

CC       = gcc
CFLAGS   = -std=c11 -Wall -Wextra -O1 -g
INC      = -Iinclude
LDFLAGS  = -lm

SRCDIR   = src
TESTDIR  = tests
EXDIR    = examples

# Source files
SRCS     = $(SRCDIR)/sysid_types.c \
           $(SRCDIR)/sysid_models.c \
           $(SRCDIR)/sysid_estimation.c \
           $(SRCDIR)/sysid_validation.c \
           $(SRCDIR)/sysid_signals.c \
           $(SRCDIR)/sysid_subspace.c \
           $(SRCDIR)/sysid_adaptive.c \
           $(SRCDIR)/sysid_nonlinear.c \
           $(SRCDIR)/sysid_frequency.c

OBJS     = $(SRCS:.c=.o)

# Headers
HDRS     = include/sysid_types.h \
           include/sysid_models.h \
           include/sysid_estimation.h \
           include/sysid_validation.h \
           include/sysid_signals.h \
           include/sysid_subspace.h \
           include/sysid_adaptive.h \
           include/sysid_frequency.h \
           include/sysid_nonlinear.h

# Test binary
TEST_BIN = $(TESTDIR)/test_sysid
TEST_SRC = $(TESTDIR)/test_sysid.c

# Example binaries
EX_BINS  = $(EXDIR)/example_dc_motor \
           $(EXDIR)/example_rc_circuit \
           $(EXDIR)/example_mass_spring

.PHONY: all clean test examples check

all: libsysid.a test examples demos benches

# Static library
libsysid.a: $(OBJS)
	ar rcs $@ $(OBJS)

# Object compilation
$(SRCDIR)/%.o: $(SRCDIR)/%.c $(HDRS)
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

# Test binary
$(TEST_BIN): $(TEST_SRC) $(OBJS)
	$(CC) $(CFLAGS) $(INC) $< $(OBJS) $(LDFLAGS) -o $@

test: $(TEST_BIN)
	@echo "Running system identification tests..."
	$(TEST_BIN)

# Examples
$(EXDIR)/example_dc_motor: $(EXDIR)/example_dc_motor.c $(OBJS)
	$(CC) $(CFLAGS) $(INC) $< $(OBJS) $(LDFLAGS) -o $@

$(EXDIR)/example_rc_circuit: $(EXDIR)/example_rc_circuit.c $(OBJS)
	$(CC) $(CFLAGS) $(INC) $< $(OBJS) $(LDFLAGS) -o $@

$(EXDIR)/example_mass_spring: $(EXDIR)/example_mass_spring.c $(OBJS)
	$(CC) $(CFLAGS) $(INC) $< $(OBJS) $(LDFLAGS) -o $@

examples: $(EX_BINS)

# Demos
DEMO_BIN = demos/demo_frequency_response

$(DEMO_BIN): demos/demo_frequency_response.c $(OBJS)
	$(CC) $(CFLAGS) $(INC) $< $(OBJS) $(LDFLAGS) -o $@

demos: $(DEMO_BIN)

# Benchmarks
BENCH_BIN = benches/bench_estimation

$(BENCH_BIN): benches/bench_estimation.c $(OBJS)
	$(CC) $(CFLAGS) $(INC) $< $(OBJS) $(LDFLAGS) -o $@

benches: $(BENCH_BIN)

# Quick check: compile only
check:
	$(CC) $(CFLAGS) $(INC) -c $(SRCS) $(LDFLAGS)

# Line count
count:
	@echo "=== Line Count (include + src) ==="
	@find include -name "*.h" -exec cat {} + | wc -l
	@find src -name "*.c" -exec cat {} + | wc -l
	@echo "=== Total ==="
	@find include src -name "*.c" -o -name "*.h" | xargs cat | wc -l

# Safety scans (per SKILL.md §10)
safety:
	@echo "=== Filler Detection ==="
	@grep -rn "_fn[0-9]" --include="*.c" --include="*.h" src include 2>/dev/null || echo "  _fnN: CLEAN"
	@grep -rn "_aux[0-9]" --include="*.c" --include="*.h" src include 2>/dev/null || echo "  _auxN: CLEAN"
	@grep -rn "_ext[0-9][^a-z]" --include="*.c" --include="*.h" src include 2>/dev/null || echo "  _extN: CLEAN"
	@echo "=== Stub Detection ==="
	@grep -rn "TODO\|FIXME\|stub\|placeholder" --include="*.c" --include="*.h" src include tests examples 2>/dev/null || echo "  CLEAN"
	@echo "=== Empty File Detection ==="
	@find src tests examples -name "*.c" -size -200c 2>/dev/null || echo "  CLEAN"

clean:
	rm -f $(OBJS) libsysid.a $(TEST_BIN) $(EX_BINS) $(DEMO_BIN) $(BENCH_BIN)
