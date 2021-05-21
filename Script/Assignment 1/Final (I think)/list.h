/* 	Module for creating and managing linked list of points
	created for ENGR30003 Numerical Programming for Engineers 2017
	by Harfiyanto Dharma Santoso <harfiyantos@student.unimelb.edu.au>
	Student No: 772503
	File : list.h
*/
#ifndef LIST_H
#define LIST_H

typedef struct list List;
typedef struct point Point;
typedef struct cell Cell;

// structure for a point
struct point {
    float x, y, u, v;
    Point* next;
};

// cell structure to store points in a grid
struct cell {
	float x, y, u, v, s;
	int k;
};

// a list points to its first and last nodes, and stores its size (num. nodes)
struct list {
	Point* head;
	Point* tail;
	int size;
};

// create a new list of points and return its address
List* new_list();

// create a new point and return its address
Point* new_point();

// delete a list and free the memory used
void free_list(List *list);

// free the memory of a node
void free_point(Point* point);

// add a point to the start of the list
void list_add_start(List* list, Point* point);

// add a point to the end of the list 
void list_add_end(List* list, Point* point);


#endif