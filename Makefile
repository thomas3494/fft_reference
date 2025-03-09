CC     ?= gcc
CFLAGS ?= -Ofast -march=native -mtune=native
LFLAGS  = -lfftw3 -lm

all: $(patsubst src/%.c, bin/%, $(wildcard src/*.c))

bin/%: src/%.c
	$(CC) $(CFLAGS) $< -o $@ $(LFLAGS)

clean:
	$(RM) bin/*
