# Makefile for Assignment 1 Numerical Programming

CC     = gcc
CFLAGS = -Wall -std=c99
EXE    = exec
OBJ    = main.o tasks.o
# add any new files here ^

$(EXE): $(OBJ)
	$(CC) $(CFLAGS) -o $(EXE) $(OBJ)

main.o: tasks.h
tasks.o: tasks.h
list.o: tasks.h

clean:
	rm -f $(OBJ) $(EXE)