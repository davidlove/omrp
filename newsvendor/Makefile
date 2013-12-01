# Makefile for newsvendor

CC=g++
CFLAGS=-ggdb -O3 -Wall
# CFLAGS=-O3 -Wall
LDFLAGS=-lm
SOURCES=multiplerep.cpp circ_sort.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=newsvendor

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

multiplerep.o:	  multiplerep.cpp multiplerep.h
	$(CC) $(CFLAGS) -c $< -o $@

circ_sort.o: circ_sort.cpp circ_sort.h
	$(CC) $(CFLAGS) -c $< -o $@

# .cpp.o: phys.h numc.h
# 	$(CC) $(CFLAGS) -c $< -o $@

parser:	output_parser.cpp
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f *.o

