/* Skeleton code for Gauss-Seidel/Jacobi method
   Must be slightly modified for Point Jacobi method
   to sovle system of linear equation [A]{X} = {C}
   accepting input data from in_gs_jacobi.csv
*/


#include<stdio.h>
#include<stdlib.h>
#include<math.h>

// functions for iterative methods
void GS_solve();

void Jacobi_solve();

int main(int argc, char *argv[]) {
    
	// Parse the input file here
	if (argc != 2) {
		printf("wrong input. correct usage: file_name\n");
		exit(EXIT_FAILURE);
	}
	char* f_src = argv[1];
	int degree = 0;
	
	// open data file for reading coefficients
	FILE *f_in = fopen(f_src, "r");
	
	//  - then, read [A] and {C} from input file
	int i = 0, j, k;
	char header[20];
	fscanf(f_in, "%s", header);
	while(*(header+i) != '\0') {
		if (*(header+i) == ',') {
			degree++;
		}
		i++;
	}

	// building A, X, C for the equation [A]{X} = {C}
	//  - first, creating empty matrix & vectors 	
	float A[degree][degree];
    float X[degree];
    float C[degree];

	
	for (i = 0; i < degree; i++) {
		for (j = 0; j < degree; j++) {
			fscanf(f_in,"%f,", A[i]+j);
		}
		fscanf(f_in,"%f", C+i);
	}
	
	// create augmented matrix
	float AUG[degree][degree+1];
	for (i = 0; i < degree; i++) {
		for (j = 0; j < degree; j++) {
			AUG[i][j] = A[i][j];
		}
		AUG[i][degree] = C[i];
	}

	// print the augmented matrix before .
	
	printf("\nAugmented Matrix (before): \n");
	for (i = 0; i < degree; i++) {
		printf("|");
		for (j = 0; j <= degree; j++) {
			printf(" %f ", AUG[i][j]);
		}
		printf("|\n");
	}
	
	// zero each column except the leading one
	float ratio = 0;
	i = 0;
	for (i = 0; i < degree; i++) {
		for (j = i+1; j < degree; j++) {
			ratio = AUG[j][i] / AUG[i][i];
			for (k = 0; k <= degree; k++) {
				AUG[j][k] -= ratio * AUG[i][k];
			}
		}
	}
	
	/*
	}
	for (j = i+1; j < degree; j++) {
			ratio = AUG[j-1][i] / AUG[j][i];
			printf("ratio = %f\n", ratio);
			for (k = 0; k < degree; k++) {
				AUG[j][k] -= ratio * AUG[j][k];
			}
		}
	*/

	printf("\nAugmented Matrix (after): \n");
	for (i = 0; i < degree; i++) {
		printf("|");
		for (j = 0; j <= degree; j++) {
			printf(" %f ", AUG[i][j]);
		}
		printf("|\n");
	}

	
	for (i = degree - 1; i >= 0; i--) {
		X[i] = AUG[i][degree];
		for (j = i + 1; j < degree; j++) {
			X[i] -= X[j] * AUG[i][j];
		}
		X[i] = X[i] / AUG[i][i];

		printf("X[%d] : %f \n", i, X[i]);
	}
	


	// print to solution
	fclose(f_in);
    return 0;
}
