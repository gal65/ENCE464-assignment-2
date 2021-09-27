CFLAGS = -O3 -Ofast -lpthread -Wall -march=native
TARGET = main

%: %.c
	gcc $(CFLAGS) $< -o $@

all: $(TARGET)

prof: CFLAGS += -pg
prof: $(TARGET)

cache: CFLAGS += -g
cache: $(TARGET)

clean:
	rm -f $(TARGET)

.PHONY: all prof cache clean
