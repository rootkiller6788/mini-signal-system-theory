CC = gcc
CFLAGS = -Wall -Wextra -std=c11 -O2 -g -Iinclude
LDFLAGS = -lm

SRCS = src/signal_basis.c src/signal_ops.c src/signal_correlation.c src/signal_decomposition.c src/signal_energy.c src/signal_applications.c src/signal_random.c
OBJS = src/signal_basis.o src/signal_ops.o src/signal_correlation.o src/signal_decomposition.o src/signal_energy.o src/signal_applications.o src/signal_random.o

TESTS = tests/test_signal_basis.c tests/test_signal_ops.c tests/test_signal_correlation.c tests/test_signal_decomposition.c tests/test_signal_energy.c
EXAMPLES = examples/example_audio_tone.c examples/example_biomedical.c examples/example_preamble_detection.c examples/example_radar_pulse.c

.PHONY: all clean test examples bench

all: test examples

%.o: %.c include/*.h
	$(CC) $(CFLAGS) -c $< -o $@

tests/%.test: tests/%.c $(OBJS)
	$(CC) $(CFLAGS) $< $(OBJS) -o $@ $(LDFLAGS)

examples/%.exe: examples/%.c $(OBJS)
	$(CC) $(CFLAGS) $< $(OBJS) -o $@ $(LDFLAGS)

test: $(OBJS) $(TESTS:.c=.test)
	@echo "=== Running Tests ==="
	@for t in $(TESTS:.c=.test); do echo "--- $$t ---"; ./$$t || exit 1; done
	@echo "=== All Tests Passed ==="

examples: $(OBJS) $(EXAMPLES:.c=.exe)
	@echo "=== Running Examples ==="
	@for e in $(EXAMPLES:.c=.exe); do echo "--- $$e ---"; ./$$e; done
	@echo "=== All Examples Complete ==="

bench: $(OBJS) benches/signal_bench.c
	$(CC) $(CFLAGS) benches/signal_bench.c $(OBJS) -o benches/signal_bench.bench $(LDFLAGS)
	./benches/signal_bench.bench

clean:
	rm -f src/*.o tests/*.test examples/*.exe benches/*.bench
