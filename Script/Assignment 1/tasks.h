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

void print_point(FILE* file, Point* p);

void maxveldiff(List* list);

void coarsegrid(List* list, int res);

int grid_index(float x, float y, int res);

void velstat(List* list);

void wakevis(List* list);

#endif
