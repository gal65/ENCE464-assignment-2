CFLAGS = -O3 -lpthread -Wall
TARGET = main

%: %.c
	gcc $(CFLAGS) $< -o $@

all: $(TARGET)

prof: CFLAGS += -pg
prof: $(TARGET)

clean:
	rm -f $(TARGET)

.PHONY: all prof clean
