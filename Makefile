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
	$(CC) $< -o $@ $(CFLAGS)

clean:
	$(RM) *.o *~ $(EXECUTABLES)
