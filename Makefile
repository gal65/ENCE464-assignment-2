CFLAGS = -O3 -Ofast -lpthread -Wall -march=native
TARGET = main

%: %.c
	gcc $(CFLAGS) $< -o $@

all: $(TARGET)

# Make with profiling hooks
prof: CFLAGS += -pg
prof: $(TARGET)

# Make with debug symbols
dbg: CFLAGS += -g
dbg: $(TARGET)

cachegrind:
	valgrind --tool=cachegrind ./$(TARGET) -n 101 -i 100 -t 4

callgrind:
	valgrind --tool=callgrind --cache-sim=yes --dump-instr=yes --collect-jumps=yes --callgrind-out-file=./callgrind.out.latest ./$(TARGET) -n 101 -i 100 -t 4
	kcachegrind ./callgrind.out.latest

clean:
	rm -f $(TARGET)

.PHONY: all prof dbg clean cachegrind callgrind
