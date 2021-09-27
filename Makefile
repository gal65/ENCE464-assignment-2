CFLAGS = -O3 -Ofast -lpthread -Wall -march=native --std=gnu99
CFLAGS_NO_OPT = -Og -lpthread -Wall -march=native --std=gnu99

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

# Make with no optimization and debug symbols
dbg_no_opt: CFLAGS = $(CFLAGS_NO_OPT) -g
dbg_no_opt: $(TARGET)

cachegrind:
	valgrind --tool=cachegrind ./$(TARGET) -n 101 -i 100 -t 4

callgrind:
	valgrind --tool=callgrind \
		--cache-sim=yes \
		--dump-instr=yes \
		--collect-jumps=yes \
		--callgrind-out-file=./callgrind.out.latest \
		./$(TARGET) -n 101 -i 100 -t 4
	kcachegrind ./callgrind.out.latest

clean:
	rm -f $(TARGET)

.PHONY: all prof dbg clean cachegrind callgrind
