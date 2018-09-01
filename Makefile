SHELL := /bin/bash


# ============================================
# COMMANDS

CC = gcc
RM = rm -f
CFLAGS=-lm -O3 -Wall

# ==========================================
# TARGETS

EXECUTABLES = classic serial parallel compare

all: $(EXECUTABLES)

classic: src/classic.c
	$(CC) $< helpers/help_methods.c -o $@ $(CFLAGS)

serial: src/serial.c
	$(CC) $< -o $@ $(CFLAGS)

parallel: src/parallel.c
	$(CC) $< -o $@ $(CFLAGS) -fopenmp

compare: helpers/compare.c
	$(CC) $< -o $@ $(CFLAGS) 

clean:
	$(RM) *.o *~ $(EXECUTABLES)
