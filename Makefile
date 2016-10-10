CFLAGS = -Wall -std=c11 -ggdb3
LDFLAGS = -lm

all: swrend

swrend: swrend.c
	gcc $(CFLAGS) $(LDFLAGS) $^ -o $@
