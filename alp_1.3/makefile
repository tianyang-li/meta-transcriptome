CC = g++
CFLAGS = -Wall -ggdb

all: alp13
	
alp13: src/*
	cd src && $(CC) $(CFLAGS) -c *.cpp
	$(CC) $(CFLAGS) src/*.o -o alp13
