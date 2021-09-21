%: %.c
	gcc -pthread $< -o $@

all: main

.PHONY: all
