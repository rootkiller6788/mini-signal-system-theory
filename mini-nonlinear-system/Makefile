# Makefile for mini-nonlinear-system
CC       = gcc
CFLAGS   = -std=c11 -Wall -Wextra -O2 -g -Iinclude
LDFLAGS  = -lm

SRCS = src/nonlinear_core.c src/nonlinear_stability.c \
       src/nonlinear_describing.c src/nonlinear_bifurcation.c \
       src/nonlinear_chaos.c src/nonlinear_phase.c \
       src/nonlinear_volterra.c src/nonlinear_applications.c
OBJS = $(SRCS:.c=.o)

.PHONY: all clean test examples count

all: test

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

tests/test_nonlinear: $(OBJS) tests/test_nonlinear.c
	$(CC) $(CFLAGS) -o $@ tests/test_nonlinear.c $(OBJS) $(LDFLAGS)

examples/van_der_pol: $(OBJS) examples/van_der_pol.c
	$(CC) $(CFLAGS) -o $@ examples/van_der_pol.c $(OBJS) $(LDFLAGS)

examples/lorenz_attractor: $(OBJS) examples/lorenz_attractor.c
	$(CC) $(CFLAGS) -o $@ examples/lorenz_attractor.c $(OBJS) $(LDFLAGS)

examples/pll_nonlinear: $(OBJS) examples/pll_nonlinear.c
	$(CC) $(CFLAGS) -o $@ examples/pll_nonlinear.c $(OBJS) $(LDFLAGS)

test: $(OBJS) tests/test_nonlinear
	./tests/test_nonlinear

examples: examples/van_der_pol examples/lorenz_attractor examples/pll_nonlinear
	@echo "Examples built successfully"

count:
	@echo "include/ + src/ total lines:"
	@find include src -name "*.h" -o -name "*.c" | xargs wc -l | tail -1

clean:
	rm -f $(OBJS)
	rm -f tests/test_nonlinear tests/test_nonlinear.exe
	rm -f examples/van_der_pol examples/van_der_pol.exe
	rm -f examples/lorenz_attractor examples/lorenz_attractor.exe
	rm -f examples/pll_nonlinear examples/pll_nonlinear.exe
