# Makefile for propagator

CC=g++
CFLAGS=-ggdb -O3 -Wall
# CFLAGS=-O3 -Wall
LDFLAGS=-lm

parser:	output_parser.cpp
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f *.o

