CC = gcc
CFLAGS = -Wall -O0 -ggdb

all: rand-cov-sim

rand-cov-sim: main.o
	$(CC) $(CFLAGS) main.o -o rand-cov-sim

main.o: main.c
	$(CC) $(CFLAGS) main.c

