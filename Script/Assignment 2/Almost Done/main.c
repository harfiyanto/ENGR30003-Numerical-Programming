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

int main(int argc, char *argv[]) {
	
	/* TODO: Parse Command Line Arguments
	DONOT explicitly set arguments to filenames */
	if (argc != 6) {
		printf("Wrong input.");
		exit(EXIT_FAILURE);
	} 
	
	char* q2_file = argv[1];
	char* q4_file = argv[2];
	char* q5_file = argv[3];
	double xo = atof(argv[4]);
	char* q6_file = argv[5];

	/* TODO: Add timing for each task and output running time in ms */
	/* Question 2 */
	shockwave(q2_file);
	/* Question 4 */
	linalgbsys(q4_file);
	
	/* Question 5 */
	interp(q5_file,xo);
	
	/* Question 6 */
	heateqn(q6_file);
    
	return (EXIT_SUCCESS);
}
