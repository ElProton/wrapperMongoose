#ifndef _SPASM_H
#define _SPASM_H

#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include <stdint.h>
#include <assert.h>
#include <stdbool.h>

/* --- primary SpaSM routines and data structures --- */

typedef int spasm_GFp;

typedef struct {                /* matrix in compressed-sparse row format */
	int nzmax;                    /* maximum number of entries */
	int n;                        /* number of rows */
	int m;                        /* number of columns */
	int *p;                       /* row pointers (size n+1) */
	int *j;                       /* column indices, size nzmax */
	spasm_GFp *x;                 /* numerical values, size nzmax (optional) */
	int prime;
}      spasm;

typedef struct {                /* matrix in triplet form */
	int nzmax;                    /* maximum number of entries */
	int nz;                       /* # entries */
	int n;                        /* number of rows */
	int m;                        /* number of columns */
	int *i;                       /* row indices, size nzmax */
	int *j;                       /* column indices (size nzmax) */
	spasm_GFp *x;                 /* numerical values, size nzmax (optional) */
	int prime;
}      spasm_triplet;


/* example (this is Matrix/t1)

		[ 4.5  0.0  3.2  0.0 ]
		[ 3.1  2.9  0.0  0.9 ]
A = [ 0.0  1.7  3.0  0.0 ]
		[ 3.5  0.4  0.0  1.0 ]

Triplet form (nz != -1) :

i = {   2,   1,   3,   0,   1,   3,   3,   1,   0,   2 }
j = {   2,   0,   3,   2,   1,   0,   1,   3,   0,   1 }
x = { 3.0, 3.1, 1.0, 3.2, 2.9, 3.5, 0.4, 0.9, 4.5, 1.7 }

the coefficients may appear in any order.

Compressed Row form :

p = {   0,             3,             6,        8,      10 }
i = {   0,   1,   3,   1,   2,   3,   0,   2,   1,   3 }
x = { 4.5, 3.1, 3.5, 2.9, 1.7, 0.4, 3.2, 3.0, 0.9, 1.0 }

In particular, the actual number of nnz is p[n]. Coefficients of a row need not be sorted by column index.

The numerical values are optional (useful for storing a sparse graph, or the pattern of a matrix). */

int spasm_nnz(const spasm * A);
void *spasm_malloc(size_t size);
void *spasm_calloc(size_t count, size_t size);
void *spasm_realloc(void *ptr, size_t size);

spasm_triplet *spasm_triplet_alloc(int n, int m, int nzmax, int prime, int with_values);
void spasm_triplet_free(spasm_triplet * A);
void spasm_add_entry(spasm_triplet * T, int i, int j, spasm_GFp x);
spasm *spasm_compress(const spasm_triplet * T);

spasm *spasm_csr_alloc(int n, int m, int nzmax, int prime, int with_values);
void spasm_csr_free(spasm * A);

spasm_triplet *spasm_load_mm(FILE * f, int prime);
void spasm_save_csr(FILE * f, const spasm * A);

/* utilities */
static inline int spasm_max(int a, int b) {
	return (a > b) ? a : b;
}

static inline int spasm_min(int a, int b) {
	return (a < b) ? a : b;
}

static inline void spasm_swap(int *a, int i, int j) {
	int x = a[i];
	a[i] = a[j];
	a[j] = x;
}

static inline int spasm_row_weight(const spasm * A, int i) {
	int *Ap = A->p;
	return Ap[i + 1] - Ap[i];
}
#endif
