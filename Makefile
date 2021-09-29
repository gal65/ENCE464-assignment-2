MAKEFLAGS += -j4

CFLAGS = -O3 -Ofast -fwhole-program -flto -pthread -Wall -Wno-maybe-uninitialized -march=native --std=gnu99
CFLAGS_NO_OPT = -Og -pthread -Wall -Wno-maybe-uninitialized --std=gnu99

CFLAGS += $(PFLAGS)
CFLAGS_NO_OPT += $(PFLAGS)

TARGET = main
A_TARGET = main_a
B_TARGET = main_b

# arguments passed to main when running test w/ cachegrind:
TEST_ARGS = -n 151 -i 100 -t $(shell nproc --all)

super_optimized: CFLAGS+= -fprofile-use
super_optimized: $(TARGET)

all: $(TARGET)

clean:
	rm -f $(TARGET) $(A_TARGET) $(B_TARGET)

full_clean: clean
	rm -f {cache,call}grind.out.* gmon.out

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

# Generate executable that will output PGO file
pgo_gen: CFLAGS += -fprofile-generate -fprofile-correction
pgo_gen: $(TARGET)

# Make executable with a PGO file
pgo_use: CFLAGS += -g -fprofile-use
pgo_use: $(TARGET)

pgo_full:
	rm -f $(TARGET) $(TARGET).gcda
	gcc $(CFLAGS) -fprofile-generate -fprofile-correction $(TARGET).c -o $(TARGET)
	./$(TARGET) $(TEST_ARGS)
	rm -f $(TARGET)
	gcc $(CFLAGS) -g -fprofile-use $(TARGET).c -o $(TARGET)

cachegrind:
	valgrind --tool=cachegrind ./$(TARGET) $(TEST_ARGS)

callgrind:
	valgrind --tool=callgrind \
		--cache-sim=yes \
		--dump-instr=yes \
		--collect-jumps=yes \
		--callgrind-out-file=./callgrind.out.latest \
		./$(TARGET) $(TEST_ARGS)
	kcachegrind ./callgrind.out.latest

build_a:
	gcc $(CFLAGS) -g -DA  $(TARGET).c -o $(A_TARGET)

build_b: 
	gcc $(CFLAGS) -g -DB $(TARGET).c -o $(B_TARGET)

checkbuild_a:
	gcc $(CFLAGS) -g -DA  $(TARGET).c -o $(A_TARGET)
	valgrind --tool=callgrind \
		--cache-sim=yes \
		--dump-instr=yes \
		--collect-jumps=yes \
		--callgrind-out-file=./callgrind.out.a_latest \
		./$(A_TARGET) $(TEST_ARGS)
	kcachegrind ./callgrind.out.a_latest

checkbuild_b:
	gcc $(CFLAGS) -g -DB $(TARGET).c -o $(B_TARGET)
	valgrind --tool=callgrind \
		--cache-sim=yes \
		--dump-instr=yes \
		--collect-jumps=yes \
		--callgrind-out-file=./callgrind.out.b_latest \
		./$(B_TARGET) $(TEST_ARGS)
	kcachegrind ./callgrind.out.b_latest

abtest: checkbuild_a checkbuild_b

abcheck: build_a build_b
	./test.sh ./main_a
	./test.sh ./main_b

.PHONY: \
super_optimized \
all \
prof \
dbg \
clean \
cachegind_clean \
cachegrind \
callgrind \
checkbuild_a \
checkbuild_b \
abtest \
abcheck \
pgo_use \
pgo_gen \
pgo_full
