/*
 * Tema 2 ASC
 * 2022 Spring
 */
#include "utils.h"
#include <stdlib.h>
#include <stdio.h>

/*
	optimized implementation for
	C = B x A x At + Bt x B,
	where At stands for the transpose form of matrix A
*/
double* my_solver(int N, double *A, double* B) {
	register int dim = N * N;
	register int i, j, k;

	double *BtB = (double *)malloc(dim * sizeof(*BtB));
	double *AAt = (double *)malloc(dim * sizeof(*AAt));
	double *C = (double *)malloc(dim * sizeof(*C));

	if (!BtB || !AAt || !C) {
		printf("\nERROR allocating memory for matrix\n\n");
		free(AAt);
		free(BtB);
		free(C);
		return NULL;
	}

	// BtB = Bt x B
	for (i = 0; i < N; i++) {
		register double *BtB_pt = BtB + i * N;
		register double *B_pt = B + i;
		for (j = 0; j < N; j++, BtB_pt++) {
			register double *B_pt1 = B_pt;
			register double *B_pt2 = B + j;
			register double sum = 0;
			for (k = 0; k < N; ++k, B_pt1 += N, B_pt2 += N) {
				sum += (*B_pt1) * (*B_pt2);
			}
			*BtB_pt = sum;
		}
	}

	// AAt = A x At
	for (i = 0; i < N; i++) {
		register double *AAt_pt = AAt + i * N;
		register double *A_pt1 = A + i * N;
		for (j = 0; j < N; j++, AAt_pt++) {
			register double *A_pt = A_pt1 + i;
			register double *At_pt = A + j * N + i;
			register double sum = 0;
			for (k = i; k < N; k++, A_pt++, At_pt++) {
				sum += (*A_pt) * (*At_pt);
			}
			*AAt_pt = sum;
		}
	}

	// BAAt = B x AAt
	for (i = 0; i < N; i++) {
		register double *BAAt_pt = C + i * N;
		register double *B_pt1 = B + i * N;
		for (j = 0; j < N; j++, BAAt_pt++) {
			register double sum = 0;
			register double *B_pt2 = B_pt1;
			register double *AAt_pt = AAt + j;
			for (k = 0; k < N; ++k, ++B_pt2, AAt_pt += N) {
				sum += (*B_pt2) * (*AAt_pt);
			}
			*BAAt_pt = sum;
		}
	}

	// C = BAAt + BtB
	for (i = 0; i < N; i++) {
		register double *C_pt = C + i * N;
		register double *BtB_pt = BtB + i * N;
		for (j = 0; j < N; j++, C_pt++, BtB_pt++) {
			*C_pt += (*BtB_pt);
		}
	}

	free(AAt);
	free(BtB);

	return C;
}
