SHELL := /bin/bash


# ============================================
# COMMANDS

CC = gcc
RM = rm -f
CFLAGS=-lm -O3

# ==========================================
# TARGETS


EXECUTABLES = serial

all: $(EXECUTABLES)

serial: serial.c
	$(CC) $< help_methods.c -o $@ $(CFLAGS)

clean:
	$(RM) *.o *~ $(EXECUTABLES)
