MAKEFLAGS += -j4

CFLAGS = -O3 -Ofast -lpthread -Wall -Wno-maybe-uninitialized -march=native --std=gnu99
CFLAGS_NO_OPT = -Og -lpthread -Wall -Wno-maybe-uninitialized --std=gnu99

CFLAGS += $(PFLAGS)
CFLAGS_NO_OPT += $(PFLAGS)

TARGET = main
A_TARGET = main_a
B_TARGET = main_b

all: $(TARGET)

clean:
	rm -f $(TARGET) $(A_TARGET) $(B_TARGET)

%: %.c
	gcc $(CFLAGS) $< -o $@

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
		./$(A_TARGET) -n 101 -i 100 -t 4
	kcachegrind ./callgrind.out.latest

checkbuild_a:
	gcc $(CFLAGS) -g -DA  $(TARGET).c -o $(A_TARGET)
	valgrind --tool=cachegrind ./$(A_TARGET) -n 101 -i 100 -t 4
	valgrind --tool=callgrind \
		--cache-sim=yes \
		--dump-instr=yes \
		--collect-jumps=yes \
		--callgrind-out-file=./callgrind.out.a_latest \
		./$(A_TARGET) -n 101 -i 100 -t 4
	kcachegrind ./callgrind.out.a_latest

checkbuild_b:
	gcc $(CFLAGS) -g -DB $(TARGET).c -o $(B_TARGET)
	valgrind --tool=cachegrind ./$(B_TARGET) -n 101 -i 100 -t 4
	valgrind --tool=callgrind \
		--cache-sim=yes \
		--dump-instr=yes \
		--collect-jumps=yes \
		--callgrind-out-file=./callgrind.out.b_latest \
		./$(B_TARGET) -n 101 -i 100 -t 4
	kcachegrind ./callgrind.out.b_latest

abtest: checkbuild_a checkbuild_b

.PHONY: all prof dbg clean cachegrind callgrind checkbuild_a checkbuild_b abtest
