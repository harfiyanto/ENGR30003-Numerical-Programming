/***************************************************************************
 *
 *   File        : tasks.h
 *   Student Id  : 772503
 *   Name        : Harfiyanto Dharma Santoso
 *
 ***************************************************************************/

#ifndef TASKS_H

void shockwave(const char* q2_file);

int is_physical(float b_l, float b_u, float theta);

double f(double b, double M, double theta, double gamma);

double fprime(double b, double M, double theta, double gamma);

double newtonraph(double initial, double theta, double M, double gamma);

float*
create_row_4(float a, float b, float c, float q);

float*
create_row_2(float a, float b);

float second_lagrange(float** data, int num_data, float xo);

float cubic_interpolation(float** data, int num_data, float x);

float* solve_tridiagonal(float** matrix, int n);

void linalgbsys(const char* q4_file);

void interp(const char* q5_file, const double xo);

void heateqn(const char* q6_file);

#endif
