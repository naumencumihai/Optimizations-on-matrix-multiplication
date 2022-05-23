/*
	BLAS implementation for
	C = B x A x At + Bt x B,
	where At stands for the transpose form of matrix A
*/

#include "utils.h"
#include <cblas.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

double* my_solver(int N, double *A, double *B) {
	double *BA;
	double *BAAt;
	double *C;
	int dim = N * N;

	BA = (double *)calloc(dim, sizeof(*BA));
	BAAt = (double *)calloc(dim, sizeof(*BAAt));
	C = (double *)calloc(dim, sizeof(*C));

	// Error handling
	if (!BA || !BAAt || !C) {
		printf("\nERROR allocating memory for matrices\n\n");
		free(BA);
		free(BAAt);
		free(C);
		return NULL;
	}

	// writes B to BA, so that dtrmm() will modify BA not B
	memcpy(BA, B, dim * sizeof(double));

	// BA = B x A (A is an upper triangular matrix)
	cblas_dtrmm(CblasRowMajor,
				CblasRight,
				CblasUpper,
				CblasNoTrans,
				CblasNonUnit,
				N, N, 
				1.0,
				A, N, 
				BA, N);

	// writes BA to BAAt, so that dtrmm() will modify BAAt not BA
	memcpy(BAAt, BA, dim * sizeof(double));

	// BAAt = BA x At
	cblas_dtrmm(CblasRowMajor,
				CblasRight,
				CblasUpper,
				CblasTrans, 
				CblasNonUnit,
				N, N, 
				1.0,
				A, N, 
				BAAt, N);

	// writes BAAt to C, so that dgemm() will modify C not B
	memcpy(C, BAAt, dim * sizeof(double));

	// C = Bt x B + BAAt
	cblas_dgemm(CblasRowMajor, 
				CblasTrans, 
				CblasNoTrans, 
				N, N, N, 
				1.0, 
				B, N, 
				B, N, 
				1.0, 
				C, N);

	// frees allocated memory
	free(BA);
	free(BAAt);

	// returns result matrix (C)
	return C;
}
