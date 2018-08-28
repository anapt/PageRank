SHELL := /bin/bash


# ============================================
# COMMANDS

CC = gcc
RM = rm -f
CFLAGS=-lm -O3 -Wall

# ==========================================
# TARGETS


EXECUTABLES = serial gauss parallel

all: $(EXECUTABLES)

serial: serial.c
	$(CC) $< help_methods.c -o $@ $(CFLAGS)

gauss: gauss.c
	$(CC) $< help_methods.c -o $@ $(CFLAGS)

parallel: parallel.c
	$(CC) $< help_methods.c -o $@ $(CFLAGS) -fopenmp

clean:
	$(RM) *.o *~ $(EXECUTABLES)
