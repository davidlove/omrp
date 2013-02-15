# Makefile for newsvendor

CC=g++
CFLAGS=-ggdb -O0
# CFLAGS=-ggdb -O0 -Wall
# CFLAGS=-O3 -Wall
LDFLAGS=-lm
SOURCES=main.cpp rand01_twister.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=../decomp

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

main.o:	main.cpp main.h
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f *.o

