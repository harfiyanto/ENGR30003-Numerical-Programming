/***************************************************************************
 *
 *   File        : tasks.h
 *   Student Id  : 772503
 *   Name        : Harfiyanto Dharma Santoso
 *
 ***************************************************************************/

#ifndef TASKS_H
#define TASKS_H

#include "list.h"

// read the file into a list of points
void read_file(const char* flow_file, List* list);

// read a point
Point* read_point(FILE* file);

// print points for task 1
void print_point(FILE* file, Point* p);

// function to compare 2 cells based on S (for qsort)
int compare_cells(const void * a, const void *b);

// Task 1
void maxveldiff(List* list);

// Task 2
void coarsegrid(List* list, int res);

// Task 3
void velstat(List* list);

// Task 4
void wakevis(List* list);

#endif
