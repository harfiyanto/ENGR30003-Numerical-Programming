/***************************************************************************
 *
 *   File        : main.c
 *   Student Id  : 772503
 *   Name        : Harfiyanto Dharma Santoso
 *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>

#include "tasks.h"
#include "list.h"


int main(int argc, char *argv[]) {
	struct timeval start;
	struct timeval stop;
	float time_task_1, time_task_2, time_task_3, time_task_4;
	/* TODO: Parse Command Line Arguments */
	char* flow_file = NULL;
	int resolution = 0;

	if (argc != 3) {
		printf("Wrong input. Usage is: %s filename resolution\n", argv[0]);
		exit(EXIT_FAILURE);
	} else {
		flow_file = argv[1];
		resolution = atoi(argv[2]);
	}

	/* TODO: Add timing for each task and output running time in ms */
    

	/* Task 1: Find the maximum velocity difference */
	/* Plus read the file */
	List* list = new_list();
	list->head = NULL;
	list->tail = NULL;
	list->size = 0;
    gettimeofday(&start, NULL);
	read_file(flow_file, list);
	maxveldiff(list);
	gettimeofday(&stop, NULL);
	time_task_1 = (stop.tv_sec - start.tv_sec) * 1000.0;
	time_task_1 += (stop.tv_usec - start.tv_usec) /1000.0;
	
	
	/* Task 2: Coarser Grid */
	gettimeofday(&start, NULL);
	coarsegrid(list, resolution);
	gettimeofday(&stop, NULL);
	time_task_2 = (stop.tv_sec - start.tv_sec) * 1000.0;
	time_task_2 += (stop.tv_usec - start.tv_usec) /1000.0;

	/* Task 3: Statistics */
	gettimeofday(&start, NULL);
	velstat(list);
	gettimeofday(&stop, NULL);
	time_task_3 = (stop.tv_sec - start.tv_sec) * 1000.0;
	time_task_3 += (stop.tv_usec - start.tv_usec) /1000.0;

	/* Task 4: Wake height and visualisation */
	gettimeofday(&start, NULL);
	wakevis(list);
	gettimeofday(&stop, NULL);
	time_task_4 = (stop.tv_sec - start.tv_sec) * 1000.0;
	time_task_4 += (stop.tv_usec - start.tv_usec) /1000.0;

	printf("TASK 1: %.2f milliseconds\n", time_task_1);
	printf("TASK 2: %.2f milliseconds\n", time_task_2);
	printf("TASK 3: %.2f milliseconds\n", time_task_3);
	printf("TASK 4: %.2f milliseconds\n", time_task_4);

	free_list(list);
    
	return (EXIT_SUCCESS);
}
