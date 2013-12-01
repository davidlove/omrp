# Makefile for newsvendor

CC=g++
# CFLAGS=-ggdb #-O3
CFLAGS=-O3 -Wall
LDFLAGS=-lm 
SOURCES=skp_omrp.cpp sip_solver.cpp mersenne.cpp stoc1.cpp stoc2.cpp userintf.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=sip

$(EXECUTABLE): $(OBJECTS) defs.h
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

skp_omrp.o:  skp_omrp.cpp
	$(CC) $(CFLAGS) -c $< -o $@

sip_solver.o: sip_solver.cpp sip_solver.h defs.h
	$(CC) $(CFLAGS) -c $< -o $@

mersenne.o: mersenne.cpp randomc.h
	$(CC) $(CFLAGS) -c $< -o $@

stoc1.o: stoc1.cpp stocc.h
	$(CC) $(CFLAGS) -c $< -o $@

stoc2.o: stoc2.cpp stocc.h
	$(CC) $(CFLAGS) -c $< -o $@

stoc3.o: stoc3.cpp stocc.h
	$(CC) $(CFLAGS) -c $< -o $@

userintf.o: userintf.cpp
	$(CC) $(CFLAGS) -c $< -o $@

parser:	output_parser.cpp
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f *.o

