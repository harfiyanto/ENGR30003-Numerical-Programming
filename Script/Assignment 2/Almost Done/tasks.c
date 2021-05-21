/***************************************************************************
 *
 *   File        : tasks.c
 *   Student Id  : 772503
 *   Name        : Harfiyanto Dharma Santoso
 *
 ***************************************************************************/
#define MAX_HEADER 100
#define MAX_ITER 1000
#define EPSILON 1e-6
#define PI 3.14159265358979323846
#define INIT_ARRAY_SIZE 1
#define QN3_COLUMN_SIZE 4
#define QN5_COLUMN_SIZE 2
#define TRIDIAGONAL_ROW 4
#define DTOR PI/180
#define RTOD 180/PI
#define MAX_B 90
#define COT(x) 1/tan(x)
#define SHOCK_OUT_FILE "out_shock.csv"
#define SHOCK_OUT_HEADER "M,theta,beta_lower,beta_upper\n"
#define LINALSYS_FILE "out_linalsys.csv"
#define LINALSYS_HEADER "x\n"
#define INTERP_FILE "out_interp.csv"
#define EXPLICIT_FE_OUT "out_heateqn_explicit_fe.csv"
#define EXPLICIT_VE_OUT "out_heateqn_explicit_ve.csv"
#define IMPLICIT_FE_OUT "out_heateqn_implicit_fe.csv"


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>
#include <assert.h>
#include "tasks.h"

void
shockwave(const char* q2_file)
{
	int i, n = 0, array_size = INIT_ARRAY_SIZE;
	char header[MAX_HEADER];
	float temp, M, theta, b_l, b_u, beta_l, beta_u, gamma;
	float* M_list = (float*) malloc(array_size * sizeof(float));
	assert(M_list);

    FILE* file = fopen(q2_file, "r");
    assert(file);
	
    fscanf(file, "%s", header);
    fscanf(file, "%f,%f,%f,%f,%f", &M, &theta, &beta_l, &beta_u, &gamma);
    fscanf(file, "%s", header);

    // reading the data dynamically
    while(fscanf(file, "%f", &temp) == 1) {
    	if (n == array_size) {
    		array_size += array_size;
    		M_list = realloc(M_list, array_size * sizeof(float*));
    	}
    	M_list[n] = temp;
    	n++;
    }

	FILE* out = fopen(SHOCK_OUT_FILE,"w");
	assert(out);

	fprintf(out, SHOCK_OUT_HEADER);

	for(i = 0; i < n; i++) {
		theta = 0;
		while (1) {
		
			if (theta == 0) {
				b_l = newtonraph(asin(1/M), M_list[i], theta*DTOR, gamma)*RTOD;
				b_u = newtonraph(MAX_B*DTOR, M_list[i], theta*DTOR, gamma)*RTOD;
			} else {
				b_l = newtonraph(beta_l*DTOR, M_list[i], theta*DTOR, gamma)*RTOD;
				b_u = newtonraph(beta_u*DTOR, M_list[i], theta*DTOR, gamma)*RTOD;
			}
			if (is_physical(b_l, b_u, theta) == 0) {
				break;
			}
			fprintf(out,"%.4f,%.0f,%.4f,%.4f\n",M_list[i],theta,b_l,b_u);
			theta++;
		}
		fprintf(out,".\n");
		fprintf(out,".\n");
	}

	free(M_list);
	fclose(out);
	//fclose(plot);
    fclose(file);
}

int
is_physical(float b_l, float b_u, float theta) {
	if (b_l == b_u) {
		return 0;
	} else if (b_l <= theta) {
		return 0;
	} else if (b_u > 90) {
		return 0;
	} else {
		return 1;
	}
}

double
f(double beta, double M, double theta, double gamma) {
	double result = 2/tan(beta) * ((M*M) * (sin(beta) * sin(beta)) - 1) / ( (M*M) * (gamma + cos(2*beta)) + 2) - tan(theta); 
	return result;
}

double
fprime(double beta, double M, double theta, double gamma) {

	double x = (4*M*M*cos(beta)/tan(beta)*sin(beta))/((gamma + cos(2*beta))*M*M + 2);
	double y = (2*((1/tan(beta))*(1/tan(beta)) + 1)*(M*M*sin(beta)*sin(beta) - 1))/((gamma + cos(2*beta))*M*M + 2);
	double z = (4*M*M*sin(2*beta)/tan(beta)*(M*M*sin(beta)*sin(beta) - 1))/((M*M*(gamma + cos(2*beta)) + 2)*(M*M*(gamma + cos(2*beta)) + 2));
	return (x - y + z);
}

double
newtonraph(double initial, double M, double theta, double gamma) {
	int i;
	double current = initial;
	for (i = 0; i < MAX_ITER; i++) {
		if (fabs(f(current, M, theta, gamma)) < EPSILON) {
			break;
		} else {
			current = current - f(current, M, theta, gamma)/fprime(current, M, theta, gamma);
		}
	}
	return current;
}

void
linalgbsys(const char* q4_file)
{
    char header[MAX_HEADER];
	int array_size = INIT_ARRAY_SIZE;
	int i, n = 0;
	float a, b, c, q;
	float** rows = (float**) malloc(array_size * sizeof(float*));
	assert(rows);

	// opening the file and throwing away the header
    FILE* file = fopen(q4_file, "r");
    assert(file);
    fscanf(file, "%s", header);
    
    // reading the data dynamically
    while(fscanf(file, "%f,%f,%f,%f", &a, &b, &c, &q) == 4) {
    	if (n == array_size) {
    		array_size += array_size;
    		rows = realloc(rows, array_size * sizeof(float*));
    	}
    	rows[n] = create_row_4(a,b,c,q);
    	n++;
    }

    float* X = solve_tridiagonal(rows, n);

    /*
    if (n < 2) {
    	printf("not enough data\n");
    	exit(EXIT_FAILURE);
    }
    */

    FILE* out = fopen(LINALSYS_FILE,"w");

    fprintf(out,LINALSYS_HEADER);
    for (i = 0; i < n; i++) {
    	fprintf(out,"%.4f\n", X[i]);
    }

    // housekeeping....
    // free up each row
    for (i = 0; i < array_size; i++) {
    	free(rows[i]);
    }

    // free the array
    //free(row);
    free(rows);
    free(X);
    fclose(out);
    fclose(file);
}

float*
solve_tridiagonal(float** matrix, int n) {
	int i;
	float* X = (float*) malloc(n * sizeof(float));
	assert(X);

	// rewrite the matrix
    for (i = 1; i < n; i++) {
    	matrix[i][0] = matrix[i][0] - (matrix[i][2] * matrix[i - 1][1] / matrix[i - 1][0]);
    	matrix[i][3] = matrix[i][3] - (matrix[i][2] * matrix[i - 1][3] / matrix[i - 1][0]); 
    	matrix[i][2] = 0;	
    }

    // calculate the solution using back substitution
    for (i = (n - 1); i >= 0; i --) {
    	if (i == n-1) {
    		X[i] = matrix[i][3] / matrix[i][0];
    	} else {
    		X[i] = (matrix[i][3] - matrix[i][1] * X[i + 1])/matrix[i][0];
    	}
    }
    return X;
}

// helper function to create a row (array of data)

float*
create_row_2(float a, float b) {
	float* row= (float*) malloc(2 * sizeof(float));
	assert(row);

	row[0] = a;
	row[1] = b;

	return row;
}

float*
create_row_4(float a, float b, float c, float q) {
	float* row= (float*) malloc(QN3_COLUMN_SIZE * sizeof(float));
	assert(row);

	row[0] = a;
	row[1] = b;
	row[2] = c;
	row[3] = q;

	return row;
}


void
interp(const char* q5_file, const double xo)
{
    
    char header[MAX_HEADER];
	int array_size = INIT_ARRAY_SIZE;
	int i, n = 0;
	float x, fx;

    float** rows = (float**) malloc(array_size * sizeof(float*));
	assert(rows);
	
    // opening the file and throwing away the header
    FILE* file = fopen(q5_file, "r");
    assert(file);

    fscanf(file, "%s", header);
    
    // reading the data dynamically
    while(fscanf(file, "%f,%f", &x, &fx) == 2) {
    	if (n == array_size) {
    		array_size += array_size;
    		rows = realloc(rows, array_size * sizeof(float*));
    	}
		rows[n] = create_row_2(x, fx);
    	n++;
    }

    FILE* out = fopen(INTERP_FILE, "w");
    assert(out);
    fprintf(out, "lagrange\n");
	fprintf(out, "%.4f\n", second_lagrange(rows, n, xo));
	fprintf(out, "cubic\n");
	fprintf(out, "%.4f\n", cubic_interpolation(rows, n, xo));

    // housekeeping....
    // free up each row
    
    for (i = 0; i < array_size; i++) {
    	free(rows[i]);
    }

    free(rows);
    fclose(file);
    fclose(out);
    // free the array

}

float
second_lagrange(float** data, int num_data, float xo) {
	int n, i;
	float result = 0.0;
	// need 3 Ls for second order lagrange
	float* L = (float*) malloc(3 * sizeof(float));
	assert(L);
    /*
	for (n = 0; n < num_data; n++) {
		if (data[n][0] >= xo) {
			if ((n == 0) || (n == 1)) {
				n = 2;
				break;
			} else {
				break;
			}
		}
	}
    */
    
    for (n = 0; n < num_data; n++) {
        if (data[n][0] > xo) {
            if ((n == 0) || (n == 1)) {
                n = 2;
                break;
            } else {
                if (fabs(data[n][0] - xo) <= fabs(data[n-1][0] - xo)) {
                    if (n < num_data - 1) {
                        n--;
                    }
                    break;
                } else {
                    break;
                }
            }
        }
    }
    

	//printf("n: %d, fx: %f\n", n, data[n][0]);

	L[0] = (xo - data[n-1][0])*(xo - data[n][0])/((data[n-2][0] - data[n-1][0])*(data[n-2][0] - data[n][0]));
	L[1] = (xo - data[n-2][0])*(xo - data[n][0])/((data[n-1][0] - data[n-2][0])*(data[n-1][0] - data[n][0]));
	L[2] = (xo - data[n-2][0])*(xo - data[n-1][0])/((data[n][0] - data[n-2][0])*(data[n][0] - data[n-1][0]));

	//printf("L[%d]: %.4f\n", 0, L[0]);
	//printf("L[%d]: %.4f\n", 1, L[1]);
	//printf("L[%d]: %.4f\n", 2, L[2]);

	for (i = n-2; i <= n; i++) {
		result += data[i][1]*L[i - n + 2];
	}
	
	free(L);
	//printf("result is: %.4f\n", result);
	return result;

}

float
cubic_interpolation(float** data, int num_data, float xo) {
	int i, n = num_data - 1;
	float result;
	float* A = (float*) malloc((n+1) * sizeof(float));
    float* B = (float*) malloc(n * sizeof(float));
    float* D = (float*) malloc(n * sizeof(float));
    float* H = (float*) malloc(n * sizeof(float));
    assert(A);
    assert(B);
    assert(D);
    assert(H);

    for (i = 0; i < n + 1; i++) {
    	A[i] = data[i][1];
    }

    for (i = 0; i < n; i++) {
    	H[i] = data[i+1][0] - data[i][0];
    }

    int matrix_size = n + 1;

    float** matrix = (float**) malloc(matrix_size * sizeof(float*));

    matrix[0] = create_row_4(1.0, 0.0, 0.0, 0.0);
    matrix[matrix_size - 1] = create_row_4(1.0, 0.0, 0.0, 0.0);

    for(i = 1; i < matrix_size - 1; i++) {
    	matrix[i] = create_row_4(2*(H[i-1]+H[i]), H[i], H[i-1], 3*(A[i+1] - A[i])/H[i] + 3*(A[i-1] - A[i])/H[i-1]);
    }

    float* C = solve_tridiagonal(matrix, matrix_size);
    assert(C);

    // calculate the corresponding B values
    for (i = 0; i < n; i++) {
    	B[i] = (A[i+1] - A[i])/H[i] - H[i]*(2*C[i] + C[i+1])/3;
    }

    // calculate the corresponding D values
    for (i = 0; i < n; i++) {
    	D[i] = (C[i+1] - C[i])/(3*H[i]);
    }

    i = 0;
    while(data[i][0] < xo) {
    	i++;
    	if(i == num_data) {
    		break;
    	}
    }

    //printf("x is %f, coefficients: %f, %f, %f, %f\n", data[i-1][0],A[i-1], B[i-1], C[i-1], D[i-1]);
    //printf("x is %f, coefficients: %f, %f, %f, %f\n", data[i-2][0],A[i-2], B[i-2], C[i-2], D[i-2]);

    result = A[i-1] + B[i-1]*(xo - data[i-1][0]) + C[i-1]*(xo - data[i-1][0])*(xo- data[i-1][0]) + 
    	D[i-1]*(xo - data[i-1][0])*(xo - data[i-1][0])*(xo - data[i-1][0]);

	for (i = 0; i < matrix_size; i++) {
    	free(matrix[i]);
    }

    free(matrix);
    free(A);
    free(B);
    free(C);
    free(D);
    free(H);
    return result;
}

void
heateqn(const char* q6_file)
{
    char header[MAX_HEADER];
    float mu = 0.0;
    int Nx, Nt, i, j;

    FILE* file = fopen(q6_file,"r");
    assert(file);

    fscanf(file,"%s", header);
    fscanf(file,"%f,%d,%d", &mu, &Nx, &Nt);

    float dx = 1.0/Nx;
    float dt = 2.0/Nt;
    float x = 0.0;

    float** fx = (float**) malloc((Nt + 1) * sizeof(float*));

    for (i = 0; i < Nt + 1; i++) {
    	fx[i] = (float*) malloc((Nx + 1) * sizeof(float));
    }

    // setting up initial values (t = 0)
    for (i = 0; i < Nx + 1; i++) {
    	x = i * dx;
    	if (x >= 0 && x < 0.125) {
    		fx[0][i] = 0;
    	} else if (x >= 0.125 && x <= 0.375) {
    		fx[0][i] = 0.5*(1 - cos(8 * PI * (x - 0.125)));
    	} else if (x > 0.375 && x <= 1) {
    		fx[0][i] = 0;
    	} else {
    		printf("something's wrong with the x value\n");
    		exit(EXIT_FAILURE);
    	}
    }

    for (i = 1; i < Nt + 1; i++) {
    	for (j = 0; j < Nx + 1; j++) {
    		if (j == 0) {
    			fx[i][j] = fx[i-1][j] + dt * mu * (fx[i-1][j] - 2*fx[i-1][j+1] + fx[i-1][j+2])/(dx*dx);
    		} else if (j == Nx) {
    			fx[i][j] = fx[i-1][j] + dt * mu * (fx[i-1][j] - 2*fx[i-1][j-1] + fx[i-1][j-2])/(dx*dx);
    		} else {
    			fx[i][j] = fx[i-1][j] + dt * mu * (fx[i-1][j+1] - 2*fx[i-1][j] + fx[i-1][j-1])/(dx*dx);
    		}
    	}
    }

    FILE* e_ve = fopen(EXPLICIT_VE_OUT,"w");
    fprintf(e_ve,"x,f(x)\n");

    for (i = 0; i < Nx + 1; i++) {
    	fprintf(e_ve,"%.4f,%.4f\n", i*dx, fx[100][i]);
    }
    fclose(e_ve);


    for (i = 1; i < Nt + 1; i++) {
    	for (j = 0; j < Nx + 1; j++) {
    		if (j == 0 || j == Nx) {
    			fx[i][j] = fx[i-1][j];
    		} else {
    			fx[i][j] = fx[i-1][j] + dt * mu * (fx[i-1][j+1] - 2*fx[i-1][j] + fx[i-1][j-1])/(dx*dx);
    		}
    	}
    }




    FILE* e_fe = fopen(EXPLICIT_FE_OUT,"w");
    fprintf(e_fe,"x,f(x)\n");

    for (i = 0; i < Nx + 1; i++) {
    	fprintf(e_fe,"%.4f,%.4f\n", i*dx, fx[100][i]);
    }
    fclose(e_fe);

    for (i = 1; i < Nt + 1; i++) {
        free(fx[i]);
    }


    // IMPLICIT
    float D = mu * dt / (dx*dx);

    float** matrix = (float**) malloc((Nx + 1) * sizeof(float*));

    for(i = 0; i < Nx + 1; i++) {
        matrix[i] = create_row_4(0.0,0.0,0.0,0.0);
    }

    for (i = 1; i < Nt + 1; i++) {
        for (j = 0; j < Nx + 1; j++) {
            if (j == 0 || j == Nx) {
                matrix[j][0] = 1.0;
                matrix[j][1] = 0.0;
                matrix[j][2] = 0.0;
                matrix[j][3] = fx[i-1][j];
            } else {
                matrix[j][0] = 2*D+1;
                matrix[j][1] = -D;
                matrix[j][2] = -D;
                matrix[j][3] = fx[i-1][j];
            }

        }
            
        fx[i] = solve_tridiagonal(matrix, Nx + 1);
        assert(fx[i]);
    }

    // print the result
    FILE* i_fe = fopen(IMPLICIT_FE_OUT,"w");
    fprintf(i_fe,"x,f(x)\n");

    for (i = 0; i < Nx + 1; i++) {
        fprintf(i_fe,"%.4f,%.4f\n", i*dx, fx[100][i]);
    }
    fclose(i_fe);

    for (i = 0; i < Nt + 1; i++) {
    	free(fx[i]);
    }

    for(i = 0; i < Nx + 1; i++) {
        free(matrix[i]);
    }
    free(fx);
    free(matrix);
  
    fclose(file);
}
