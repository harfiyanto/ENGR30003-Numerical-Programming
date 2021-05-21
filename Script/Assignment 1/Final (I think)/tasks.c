/***************************************************************************
 *
 *   File        : tasks.c
 *   Student Id  : 772503
 *   Name        : Harfiyanto Dharma Santoso
 *
 ***************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>

#include "tasks.h"
#include "list.h"

#define BOUND 0.050000
#define THRESHOLD 0.500000
#define MIN_X_TASK_1 20.000000
#define XMIN 10.000000
#define XMAX 70.000000
#define YMIN -20.000000
#define YMAX 20.000000
#define PERCENT 100.0
#define XPOINTS 12
#define MAX_HEADER_SIZE 100

#define FNAME_TASK_1 "task1.csv"
#define FNAME_TASK_2 "task2.csv"
#define FNAME_TASK_3 "task3.csv"
#define FNAME_TASK_4_1 "task4_1.csv"

#define HEADER_TASK_1 "x,y,u,v\n"
#define HEADER_TASK_2 "x,y,u,v,S\n"
#define HEADER_TASK_3 "threshold,points,percentage\n"
#define HEADER_TASK_4_1 "x,y_h\n"


// helper function to read a file into a linked list structure
void
read_file(const char* flow_file, List* list) {
    FILE* file = fopen(flow_file, "r");
    if (file == 0) {
        printf("File can't be opened\n");
    } else {
        char* header = malloc((sizeof *header) * MAX_HEADER_SIZE);
        assert(header);

        // scan the header
        fscanf(file, "%s", header);

        // read points and insert them into the list
        while(1){
            Point* p = read_point(file);
            if (!p) {
                break;
            }
            list_add_end(list, p);
        }
        fclose(file);

        free(header);
    }
}

// function to read a point and return its address
Point*
read_point(FILE* file) {
    Point* p = new_point();
    assert(p);
    p->next = NULL;

    // scan a line into point structure and return a pointer to it if success
    if (fscanf(file, "%f,%f,%f,%f", &(p->x), &(p->y), &(p->u), &(p->v)) == 4) {
        return p;
    } else {
        free(p);
        return NULL;
    }
}

// ================ TASK 1 =================
void
maxveldiff(List* list)
{
    float max_u, max_v, min_u, min_v;
    Point *current = list->head;
    Point *max_u_p, *min_u_p, *max_v_p, *min_v_p;

    // initialise the values
    max_u_p = current;
    min_u_p = current;
    max_v_p = current;
    min_v_p = current;
    min_u = current->u;
    min_v = current->v;
    max_u = current->u;
    max_v = current->v;

    // update the maximums and the corresponding point. if two points have
    // the same interest value, choose the one with smaller y
    while(current) {
        if (current->x > MIN_X_TASK_1) {
            // update max_u (if necessary)
            if (current->u > max_u) {
                max_u = current->u;
                max_u_p = current;
            } else if ((current->u == max_u)&& (current->x < max_u_p->x)) {
                max_u_p = current;
            }
            // update min_u (if necessary)
            if (current->u < min_u) {
                min_u = current->u;
                min_u_p = current;
            } else if ((current->u == min_u)&& (current->x < min_u_p->x)) {
                min_u_p = current;
            }
            // update max_v (if necessary)
            if (current->v > max_v) {
                max_v = current->v;
                max_v_p = current;
            } else if ((current->v == max_v)&& (current->x < max_v_p->x)) {
                max_v_p = current;
            }
            // update min_v (if necessary)
            if (current->v < min_v) {
                min_v = current->v;
                min_v_p = current;
            } else if ((current->v == min_v)&& (current->x < min_v_p->x)) {
                min_v_p = current;
            }
        }
        current = current->next;
    }

    // creating the file
    FILE *file;
    file = fopen(FNAME_TASK_1, "w");
    assert(file);
    fprintf(file,HEADER_TASK_1);
    print_point(file, max_u_p);
    print_point(file, min_u_p);
    print_point(file, max_v_p);
    print_point(file, min_v_p);
    fclose(file);
}

// helper function to print points for maxveldiff
void
print_point(FILE* file, Point* p) {
    fprintf(file, "%.6f,%.6f,%.6f,%.6f\n",p->x, p->y, p->u, p->v);    
}

// ================ TASK 2 =================
void
coarsegrid(List* list, int res)
{
    assert(list);
    assert(res > 0);

    int i;
    float x, y, u, v, s;
    int grid_x, grid_y;
    double x_mult = (XMAX - XMIN) / res;
    double y_mult = (YMAX - YMIN) / res;
    Point* current = list->head;
    Cell* grid = malloc((sizeof *grid) * res * res);    // array of grids
    assert(grid);

    // initialise each grid's values
    for (i = 0; i < res * res; i++) {
        grid[i].x = 0;
        grid[i].y = 0;
        grid[i].u = 0;
        grid[i].v = 0;
        grid[i].k = 0;
        grid[i].s = 0;
    }

    // go through each point and insert them to their corresponding grid(s)
    while (current) {
        grid_x = 0;
        grid_y = 0;

        // find the grid it belongs to 
        while (((grid_x + 1) * x_mult + XMIN < current->x) && (grid_x < res)) {
            grid_x++;
        }
        while (((grid_y + 1) * y_mult + YMIN < current->y) && (grid_y < res)) {
            grid_y++;
        }

        // insert the point into the grid
        i = grid_x + res * grid_y;
        grid[i].x += current->x;
        grid[i].y += current->y;
        grid[i].u += current->u;
        grid[i].v += current->v;
        grid[i].k += 1;

        // the point also belongs to the grid beside it 
        // except if the grid is off the boundary
        if (((grid_x + 1) * x_mult + XMIN) == current->x) {
            if (grid_x <= (res - 2)) {
                i = (grid_x + 1) + res * grid_y;
                grid[i].x += current->x;
                grid[i].y += current->y;
                grid[i].u += current->u;
                grid[i].v += current->v;
                grid[i].k += 1;
            }

        }

        // the point also belongs to the grid below it
        if (((grid_y + 1) * y_mult + YMIN) == current->y) {
            if (grid_y <= (res - 2)) {
                i = grid_x + res * (grid_y + 1);
                grid[i].x += current->x;
                grid[i].y += current->y;
                grid[i].u += current->u;
                grid[i].v += current->v;
                grid[i].k += 1;
            }

        }

        // the point also belongs to the grid diagonal to it
        if ((((grid_x + 1) * x_mult + XMIN) == current->x) && (((grid_y + 1) * y_mult + YMIN) == current->y)) {
            if (grid_x <= (res - 2) && grid_y <= (res-2)) {    
                i = (grid_x + 1) + res * (grid_y + 1);
                grid[i].x += current->x;
                grid[i].y += current->y;
                grid[i].u += current->u;
                grid[i].v += current->v;
                grid[i].k += 1;
            }
        }
        current = current->next;
    }

    // create a file and write the data
    FILE* file;
    file = fopen(FNAME_TASK_2, "w");
    assert(file);
    fprintf(file, HEADER_TASK_2);

    // calculating the S
    for (i = 0; i < res * res; i++) {
        x = grid[i].x / grid[i].k;
        y = grid[i].y / grid[i].k;
        u = grid[i].u / grid[i].k;
        v = grid[i].v / grid[i].k;
        grid[i].s = 100 * sqrt(u*u + v*v) / sqrt(x*x + y*y);
    }

    // sort the array using qsort
    qsort(grid, res * res, sizeof(*grid), compare_cells);

    // print the average values and S
    for (i = 0; i < res * res; i++) {
        x = grid[i].x / grid[i].k;
        y = grid[i].y / grid[i].k;
        u = grid[i].u / grid[i].k;
        v = grid[i].v / grid[i].k;
        s = grid[i].s;
        fprintf(file, "%.6f,%.6f,%.6f,%.6f,%.6f\n",x,y,u,v,s);
    }

    free(grid);
    fclose(file);
}

// helper function to compare two cells based on S (for qsort)
int
compare_cells(const void *a, const void *b) {
    Cell c1 = *((Cell*) a);
    Cell c2 = *((Cell*) b);

    if (c1.s > c2.s) return -1;
    if (c1.s < c2.s) return 1;
    return 0;
}

// ================ TASK 3 =================
void
velstat(List *list)
{
    double threshold = THRESHOLD;
    int count = 0;
    Point* current = list->head;

    // creating the file
    FILE* file;
    file = fopen(FNAME_TASK_3, "w");
    assert(file);
    fprintf(file, HEADER_TASK_3);

    // increment the threshold and print the percentage of points within it
    while(count != list->size) {
        count = 0;
        current = list->head;
        while(current) {
            if (current->u < threshold) {
                count++;
            }
            current = current->next;
        }
        fprintf(file, "%.6f,%d,%.6f\n", threshold, count, PERCENT * count / list->size);
        threshold+=0.10000;
    }
    fclose(file);
}

// ================ TASK 4 =================
void
wakevis(List* list)
{
    assert(list);
    int i,j,c;
    int n = XPOINTS;                     // Location in x for wake visualization
    float* yheight;
    yheight = (float*) calloc(n,sizeof(float));
    assert(yheight);

    float* x = malloc((sizeof *x) * XPOINTS);           // array to keep track of the x locations
    float* max_u = malloc((sizeof *max_u) * XPOINTS);   // array to keep track of the max u for each x loc
    Point** p = malloc((sizeof *p) * XPOINTS);          // array to keep the points with closest x(and max u)
    Point* current = list->head;
    int* initialised = malloc((sizeof *initialised) * XPOINTS);   // array to keep track of uninitialized points

    assert(x);
    assert(max_u);
    assert(p);
    assert(initialised);

    // set the locations in x and initialise p and initialised array
    for(c = 0; c < n; c++) {
        x[c] = XMIN + c * 5;
        p[c] = current;
        initialised[c] = 0;
    }

    // find the maximum value of u at each of the x location
    // initialise the value if necessary
    while(current) {
        for (c = 0; c < n; c++) {
            // if the x value is within the boundary
            if ((current->x >= (x[c] - BOUND)) && (current->x <= (x[c] + BOUND))) {
                // if not initialised
                if (initialised[c] == 0) {
                    // set u to the current point's value  
                    max_u[c] = current->u;
                    initialised[c] = 1;
                    p[c] = current;
                // if initialised
                } else {
                    // find the closest x-value
                    if (fabs(current->x - x[c]) < fabs(p[c]->x - x[c])) {
                        max_u[c] = current->u;
                        p[c] = current;
                    // then the max u if the current point have the same x
                    } else if (fabs(current->x - x[c]) == fabs(p[c]->x - x[c])) {
                        if (current->u > max_u[c]) {
                            max_u[c] = current->u;
                            p[c] = current;
                        }
                    }
                }
            }
        }
        current = current->next;
    }

    // creating and printing the data into the file
    FILE* file = fopen(FNAME_TASK_4_1, "w");
    assert(file);
    fprintf(file,HEADER_TASK_4_1);

    for (c = 0; c < n; c++) {
        fprintf(file,"%.6f,%.6f\n", p[c]->x, fabs(p[c]->y));
    }
    fclose(file);

    // calculate and set the spacing at each location
    for (c = 0; c < n; c++) {
        yheight[c] = ceil(10 * fabs(p[c]->y));
    }

    // housekeeping....
    free(x);
    free(max_u);
    free(p);
    free(initialised);

    /* Task 4: Part 2, nothing is to be changed here
       Remember to output the spacing into the array yheight
       for this to work. You also need to initialize i,j and 
       yheight so the skeleton as it stands will not compile */
     
    FILE *ft42;
    ft42 = fopen("task4_2.txt","w");
    assert(ft42);
    for (j = 11; j>=0; j--){
	for (i=0;i<yheight[j]-yheight[0]+4;i++){
 	    fprintf(ft42, " ");
	}
    	fprintf(ft42, "*\n");
    }
    for (i=0;i<5; i++){
    	fprintf(ft42, "III\n");
    }
    for(j = 0; j<12; j++ ){
    	for (i=0;i<yheight[j]-yheight[0]+4;i++){
    	    fprintf(ft42, " ");
    	}
    	fprintf(ft42, "*\n");
    }
    fclose(ft42);
    
    /* Cleanup */
    free(yheight);

}
