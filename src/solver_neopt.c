/*
 * Tema 2 ASC
 * 2022 Spring
 */
#include "utils.h"
#include <stdlib.h>
#include <stdio.h>

int max(int a, int b) {
	if (a >= b){
		return a;
	} else {
		return b;
	}
}

/*
	unoptimized implementation for
	C = B x A x At + Bt x B,
	where At stands for the transpose form of matrix A
*/
double* my_solver(int N, double *A, double* B) {
	int dim = N * N;

	double *BtB = (double *)calloc(dim, sizeof(*BtB));
	if (!BtB) {
		printf("\nERROR allocating memory for matrix\n\n");
		free(BtB);
		return NULL;
	}

	// BtB = Bt x B
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < N; k++) {
				BtB[i * N + j] += B[k * N + i] * B[k * N + j];
			}
		}
	}

	double *AAt = (double *)calloc(dim, sizeof(*AAt));
	if (!AAt) {
		printf("\nERROR allocating memory for matrix\n\n");
		free(AAt);
		return NULL;
	}

	// AAt = A x At
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			for (int k = max(i, j); k < N; k++) {
				AAt[i * N + j] += A[i * N + k] * A[j * N + k];
			}
		}
	}

	double *BAAt = (double *)calloc(dim, sizeof(*BAAt));
	if (!BAAt) {
		printf("\nERROR allocating memory for matrix\n\n");
		free(BAAt);
		return NULL;
	}

	// BAAt = B x AAt
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < N; k++) {
				BAAt[i * N + j] += B[i * N + k] * AAt[k * N + j];
			}
		}
	}

	double *C = (double *)calloc(dim, sizeof(*C));
	if (!C) {
		printf("\nERROR allocating memory for matrix\n\n");
		free(C);
		return NULL;
	}

	// C = BAAt + BtB
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			C[i * N + j] = BAAt[i * N + j] + BtB[i * N + j];
		}
	}

	free(AAt);
	free(BtB);
	free(BAAt);

	return C;
}
