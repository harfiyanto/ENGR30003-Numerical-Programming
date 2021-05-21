/* 	Module for creating and managing linked list of points
	created for ENGR30003 Numerical Programming for Engineers 2017
	by Harfiyanto Dharma Santoso <harfiyantos@student.unimelb.edu.au>
	Student No: 772503
*/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "list.h"

// FUNCTION DEFINITIONS

// create a new list of points
List*
new_list() {
	List* list = malloc(sizeof *list);
	assert(list);
	
	list->head = NULL;
	list->tail = NULL;
	list->size = 0;

	return list;
}

// delete a list and free the memory used
void
free_list(List *list) {
	assert(list);
	// free each point
	Point* point = list->head;
	Point* next;
	while (point) {
		next = point->next;
		free_point(point);
		point = next;
	}
	// free the list 
	free(list);
}

// create a new point and return its address
Point*
new_point() {
	Point* point = malloc(sizeof *point);
	assert(point);
	point->next = NULL;
	return point;
}

// free the memory of a node
void
free_point(Point* point) {
	free(point);
}

// add a point to the start of the list
void
list_add_start(List* list, Point* point) {
	assert(list);
	assert(point);
	point->next = list->head;
	list->head = point;
	if (list->size == 0) {
		list->tail = point;
	}
	list->size+=1;
}

// add a point to the end of the list 
void
list_add_end(List* list, Point* point) {
	assert(list);
	assert(point);
	if (list->size == 0) {
		list->tail = point;
		list->head = point;
	} else {
		list->tail->next = point;
		list->tail = point;
	}
	list->size+=1;
}


