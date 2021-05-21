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

#define BOUND 0.05
#define XMIN 10.000000
#define XMAX 70.000000
#define YMIN -20.000000
#define YMAX 20.000000

#define FNAME_TASK_1 "task1.csv"
#define FNAME_TASK_2 "task2.csv"
#define FNAME_TASK_3 "task3.csv"
#define FNAME_TASK_4_1 "task4_1.csv"

#define HEADER_TASK_1 "x,y,u,v\n"
#define HEADER_TASK_2 "x,y,u,v,S\n"
#define HEADER_TASK_3 "threshold,points,percentage\n"
#define HEADER_TASK_4_1 "x,y_h\n"

void
read_file(const char* flow_file, List* list) {
    FILE* file = fopen(flow_file, "r");
    if (file == 0) {
        printf("File can't be opened\n");
    } else {
        char header[100];
        fscanf(file, "%s", header);

        while(1){
            Point* p = read_point(file);
            if (!p) {
                break;
            }
            list_add_end(list, p);
        }
        fclose(file);
    }
}

// function to read a point and return its address
Point*
read_point(FILE* file) {
    Point* p = new_point();
    assert(p);
    p->next = NULL;

    if (fscanf(file, "%f,%f,%f,%f", &(p->x), &(p->y), &(p->u), &(p->v)) == 4) {
        return p;
    } else {
        free(p);
        return NULL;
    }
}


void
print_point(FILE* file, Point* p) {
    fprintf(file, "%.6f,%.6f,%.6f,%.6f\n",p->x, p->y, p->u, p->v);
}

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
        if (current->x > 20) {
            // update max_u (if necessary)
            if (current->u > max_u) {
                max_u = current->u;
                max_u_p = current;
            } else if ((current->u == max_u)&& (current->y < max_u_p->y)) {
                max_u_p = current;
            }
            // update min_u (if necessary)
            if (current->u < min_u) {
                min_u = current->u;
                min_u_p = current;
            } else if ((current->u == min_u)&& (current->y < min_u_p->y)) {
                min_u_p = current;
            }
            // update max_v (if necessary)
            if (current->v > max_v) {
                max_v = current->v;
                max_v_p = current;
            } else if ((current->v == max_v)&& (current->y < max_v_p->y)) {
                max_v_p = current;
            }
            // update min_v (if necessary)
            if (current->v < min_v) {
                min_v = current->v;
                min_v_p = current;
            } else if ((current->v == min_v)&& (current->y < min_v_p->y)) {
                min_v_p = current;
            }
        }
        current = current->next;
    }

    // creating the file
    FILE *file;
    file = fopen(FNAME_TASK_1, "w");
    fprintf(file,HEADER_TASK_1);
    print_point(file, max_u_p);
    print_point(file, min_u_p);
    print_point(file, max_v_p);
    print_point(file, min_v_p);
    fclose(file);
}

void
coarsegrid(List* list, int res)
{
    assert(list);
    assert(res > 0);
    int i;
    float x, y, u, v, s;
    //float avg_x, avg_y, avg_u, avg_v, s;
    Point* current = list->head;

    Cell* grid = malloc((sizeof *grid) * res * res);

    // initialise
    for (i = 0; i < res * res; i++) {
        grid[i].x = 0;
        grid[i].y = 0;
        grid[i].u = 0;
        grid[i].v = 0;
        grid[i].k = 0;
    }

    while (current) {
        i = grid_index(current->x, current->y, res);
        grid[i].x += current->x;
        grid[i].y += current->y;
        grid[i].u += current->u;
        grid[i].v += current->v;
        grid[i].k += 1;
        current = current->next;
    }

    // create a file and write the data
    FILE* file;
    file = fopen(FNAME_TASK_2, "w");
    assert(file);
    fprintf(file, HEADER_TASK_2);

    for (i = 0; i < res * res; i++) {
        x = grid[i].x / grid[i].k;
        y = grid[i].y / grid[i].k;
        u = grid[i].u / grid[i].k;
        v = grid[i].v / grid[i].k;
        s = 100 * sqrt(u*u + v*v) / sqrt(x*x + y*y);
        fprintf(file, "%.6f,%.6f,%.6f,%.6f,%.6f\n",x,y,u,v,s);
    }
    fclose(file);
    free(grid);
}

// return the corresponding index of cells a point belong to
int
grid_index(float x, float y, int res) {
    assert(res > 0);
    assert(x <= XMAX);
    assert(x >= XMIN);
    assert(y <= YMAX);
    assert(y >= YMIN);
    int grid_x = 0, grid_y = 0;
    float x_mult = (XMAX - XMIN) / res;
    float y_mult = (YMAX - YMIN) / res;
    while(grid_x < (res - 1)) {
        if (x >= (grid_x * x_mult + XMIN)) {
            if (x <= ((grid_x + 1) * x_mult + XMIN)) {
                break;
            }
        }
        grid_x++;
    }
    
    while(grid_y < (res - 1)) {
        if (y >= (grid_y * y_mult + YMIN)) {
            if (y <= ((grid_y + 1) * y_mult + YMIN)) {
                break;
            }
        }
        grid_y++;
    }
    return (grid_x + res * grid_y);
}

void
velstat(List *list)
{
    float threshold = 0;
    int count = 0;
    Point* current = list->head;
    // creating the file
    FILE* file;
    file = fopen(FNAME_TASK_3, "w");
    assert(file);
    fprintf(file, HEADER_TASK_3);
    while(count != list->size) {
        count = 0;
        current = list->head;
        while(current) {
            if (current->u < threshold) {
                count++;
            }
            
            current = current->next;
        }
        if (count>0) {
            fprintf(file, "%.6f,%d,%.6f\n", threshold, count, 100.0 * count / list->size);
        }
        threshold+=0.1;
    }
    fclose(file);
}

void
wakevis(List* list)
{
    assert(list);
    int i,j,c;
    int n = 12; // Location in x for wake visualization
    float* yheight;
    yheight = (float*) calloc(n,sizeof(float));

    float x[n], max_u[n];
    Point* p[n];
    int initialised[n];
    memset(initialised, 0, n);


    Point* current = list->head;
    // set the locations in x
    for(c = 0; c < n; c++) {
        x[c] = XMIN + c * 5;
        p[c] = current;
    }

    // find the maximum value of u at each of the x location
    // initialise the value if necessary
    
    while(current) {
        for (c = 0; c < n; c++) {
            // if the x value is within the boundary
            if ((current->x >= (x[c] - BOUND)) && (current->x <= (x[c] + BOUND))) {
                // if not initialised
                if (initialised[c] != 1) {
                    // set u to the current point's value  
                    max_u[c] = current->u;
                    initialised[c] = 1;
                // if initialised
                } else {
                    if (current->u > max_u[c]) {
                        max_u[c] = current->u;
                        p[c] = current;
                    }   
                }
            }
        }
        current = current->next;
    }

    // creating and printing the data into the file
    FILE* file = fopen(FNAME_TASK_4_1, "w");
    fprintf(file,HEADER_TASK_4_1);
    for (c = 0; c < n; c++) {
        fprintf(file,"%.6f,%.6f\n", p[c]->x, p[c]->y);
    }
    fclose(file);

    // calculate and set the spacing at each location
    for (c = 0; c < n; c++) {
        yheight[c] = ceil(10 * fabs(p[c]->y));
    }

    /* Task 4: Part 2, nothing is to be changed here
       Remember to output the spacing into the array yheight
       for this to work. You also need to initialize i,j and 
       yheight so the skeleton as it stands will not compile */
     
    FILE *ft42;
    ft42 = fopen("task4_2.txt","w");
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
