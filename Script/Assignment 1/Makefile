# Makefile for Assignment 1 Numerical Programming

CC     = gcc
CFLAGS = -Wall -std=c99
EXE    = flow
OBJ    = main.o tasks.o list.o
# add any new files here ^

$(EXE): $(OBJ)
	$(CC) $(CFLAGS) -o $(EXE) $(OBJ)

main.o: tasks.h list.h
tasks.o: tasks.h list.h
list.o: tasks.h list.h

clean:
	rm -f $(OBJ) $(EXE)