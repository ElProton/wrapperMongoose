#include <assert.h>
#include <math.h>
#include <err.h>
#include <sys/time.h>

#include "mini_spasm.h"
#include "mmio.h"

double spasm_wtime() {
	struct timeval ts;
	gettimeofday(&ts, NULL);
	return (double)ts.tv_sec + ts.tv_usec / 1E6;
}


int spasm_nnz(const spasm * A) {
	return A->p[A->n];
}


void *spasm_malloc(size_t size) {
	void *x = malloc(size);
	if (x == NULL)
		err(1, "malloc failed");
	return x;
}

void *spasm_calloc(size_t count, size_t size) {
	void *x = calloc(count, size);
	if (x == NULL)
		err(1, "calloc failed");
	return x;
}

void *spasm_realloc(void *ptr, size_t size) {
	void *x = realloc(ptr, size);
	if (ptr != NULL && x == NULL && size != 0)
		err(1, "realloc failed");
	return x;
}


/* return a string representing n in 4 bytes */
static void spasm_human_format(int64_t n, char *target) {
	if (n < 1000) {
		sprintf(target, "%d", (int) n);
		return;
	}
	if (n < 1000000) {
		sprintf(target, "%.1fk", n / 1e3);
		return;
	}
	if (n < 1000000000) {
		sprintf(target, "%.1fm", n / 1e6);
		return;
	}
	if (n < 1000000000000ll) {
		sprintf(target, "%.1fg", n / 1e9);
		return;
	}
	if (n < 1000000000000000ll) {
		sprintf(target, "%.1ft", n / 1e12);
		return;
	}
}

/*
 * change the max # of entries in a sparse matrix. If nzmax < 0, then the
 * matrix is trimmed to its current nnz.
 */
void spasm_csr_realloc(spasm * A, int nzmax) {
	if (nzmax < 0)
		nzmax = spasm_nnz(A);
	A->j = spasm_realloc(A->j, nzmax * sizeof(int));
	if (A->x != NULL)
		A->x = spasm_realloc(A->x, nzmax * sizeof(spasm_GFp));
	A->nzmax = nzmax;
}

/*
 * change the max # of entries in a sparse matrix. If nzmax < 0, then the
 * matrix is trimmed to its current nnz.
 */
void spasm_triplet_realloc(spasm_triplet * A, int nzmax) {
	if (nzmax < 0)
		nzmax = A->nz;
	A->i = spasm_realloc(A->i, nzmax * sizeof(int));
	A->j = spasm_realloc(A->j, nzmax * sizeof(int));
	if (A->x != NULL)
		A->x = spasm_realloc(A->x, nzmax * sizeof(spasm_GFp));
	A->nzmax = nzmax;
}

/* free a sparse matrix */
void spasm_csr_free(spasm * A) {
	if (A == NULL)
		return;
	free(A->p);
	free(A->j);
	free(A->x);		/* trick : free does nothing on NULL pointer */
	free(A);
}

void spasm_triplet_free(spasm_triplet * A) {
	free(A->i);
	free(A->j);
	free(A->x);		/* trick : free does nothing on NULL pointer */
	free(A);
}


/* allocate a sparse matrix (triplet form) */
spasm_triplet *spasm_triplet_alloc(int n, int m, int nzmax, int prime, int with_values) {
	spasm_triplet *A;

	A = spasm_malloc(sizeof(spasm_triplet));
	A->m = m;
	A->n = n;
	A->nzmax = nzmax;
	A->prime = prime;
	A->nz = 0;
	A->i = spasm_malloc(nzmax * sizeof(int));
	A->j = spasm_malloc(nzmax * sizeof(int));
	A->x = (with_values ? spasm_malloc(nzmax * sizeof(spasm_GFp)) : NULL);
	return A;
}

/*
 * Load a matrix in MatrixMarket sparse format.
 * Heavily inspired by the example program:
 *     http://math.nist.gov/MatrixMarket/mmio/c/example_read.c
 */
spasm_triplet *spasm_load_mm(FILE * f, int prime) {
	MM_typecode matcode;
	int n, m, nnz;

	double start = spasm_wtime();
	if (mm_read_banner(f, &matcode) != 0) 
		errx(1, "Could not process Matrix Market banner.\n");

	if (!mm_is_matrix(matcode) || !mm_is_sparse(matcode))
		errx(1, "Matrix Market type: [%s] not supported", mm_typecode_to_str(matcode));
	
	int symmetric = mm_is_symmetric(matcode);
	int skew = mm_is_skew(matcode);

	if (!mm_is_general(matcode) && !symmetric && !skew)
		errx(1, "Matrix market type [%s] not supported",  mm_typecode_to_str(matcode));

	if (mm_read_mtx_crd_size(f, &n, &m, &nnz) != 0)
		errx(1, "Cannot read matrix size");

	fprintf(stderr, "[IO] loading %d x %d MTX [%s] modulo %d, %d nnz...", n, m, mm_typecode_to_str(matcode), prime, nnz);
	fflush(stderr);
	
	if (mm_is_pattern(matcode))
		prime = -1;

	spasm_triplet *T = spasm_triplet_alloc(n, m, nnz, prime, prime != -1);

	for (int i = 0; i < nnz; i++) {
		int u, v, w;
		double x, y;

		if (mm_is_pattern(matcode)) {
			if (2 != fscanf(f, "%d %d\n", &u, &v))
				errx(1, "parse error entry %d\n", i);
			spasm_add_entry(T, u - 1, v - 1, 1);
		} else if (mm_is_integer(matcode)) {
			if (3 != fscanf(f, "%d %d %d\n", &u, &v, &w))
				errx(1, "parse error entry %d\n", i);
			spasm_add_entry(T, u - 1, v - 1, w);
		} else if (mm_is_real(matcode)) {
			if (3 != fscanf(f, "%d %d %lg\n", &u, &v, &x))
				errx(1, "parse error entry %d\n", i);
			spasm_add_entry(T, u - 1, v - 1, (int) (100000 * x));
		} else if (mm_is_complex(matcode)) {
			if (4 != fscanf(f, "%d %d %lg %lg\n", &u, &v, &y, &x))
				errx(1, "parse error entry %d\n", i);
			spasm_add_entry(T, u - 1, v - 1, (int) (1000 * (y + 100 * x)));
		} else {
			errx(1, "Don't know how to read matrix");
		}
	}

	if (symmetric || skew)
		nnz *= 2;

	if (symmetric) {
		int mult = skew ? -1 : 1;
		int nz = T->nz;
		for (int px = 0; px < nz; px++)
			if (T->j[px] != T->i[px])
				spasm_add_entry(T, T->j[px], T->i[px], (T->x != NULL) ? (mult * T->x[px]) : 1);
	}

	char s_nnz[16];
	spasm_human_format(T->nz, s_nnz);
	fprintf(stderr, "%s NNZ [%.1fs]\n", s_nnz, spasm_wtime() - start);
	return T;
}


/*
 * save a matrix in SMS format. TODO : change name to spasm_csr_save
 */
void spasm_save_csr(FILE * f, const spasm * A) {
	assert(f != NULL);

	int *Aj = A->j;
	int *Ap = A->p;
	spasm_GFp *Ax = A->x;
	int n = A->n;
	int m = A->m;
	int prime = A->prime;

	fprintf(f, "%d %d M\n", n, m);
	for (int i = 0; i < n; i++)
		for (int px = Ap[i]; px < Ap[i + 1]; px++) {
			spasm_GFp x = (Ax != NULL) ? Ax[px] : 1;
			x = (x > prime / 2) ? x - prime : x;
			fprintf(f, "%d %d %d\n", i + 1, Aj[px] + 1, x);
		}
	fprintf(f, "0 0 0\n");
}

/* add an entry to a triplet matrix; enlarge it if necessary */
void spasm_add_entry(spasm_triplet * T, int i, int j, spasm_GFp x) {
	spasm_GFp x_p;

	assert((i >= 0) && (j >= 0));
	int prime = T->prime;

	if (T->nz == T->nzmax)
		spasm_triplet_realloc(T, 1 + 2 * T->nzmax);
	if (T->x != NULL) {
		x_p = ((x % prime) + prime) % prime;
		if (x_p == 0)
			return;
		T->x[T->nz] = x_p;
	}

	// fprintf(stderr, "Adding (%d, %d, %d)\n", i, j, x);

	T->i[T->nz] = i;
	T->j[T->nz] = j;
	T->nz += 1;
	T->n = spasm_max(T->n, i + 1);
	T->m = spasm_max(T->m, j + 1);
}

/* allocate a sparse matrix (compressed-row form) */
spasm *spasm_csr_alloc(int n, int m, int nzmax, int prime, int with_values) {
	spasm *A;

	if (prime > 46337) {
		prime = 46337;
		fprintf(stderr, "WARNING: modulus has been set to 46337.\n");
	}
	A = spasm_malloc(sizeof(spasm));	/* allocate the cs struct */
	A->m = m;		/* define dimensions and nzmax */
	A->n = n;
	A->nzmax = nzmax;
	A->prime = prime;
	A->p = spasm_malloc((n + 1) * sizeof(int));
	A->j = spasm_malloc(nzmax * sizeof(int));
	A->x = with_values ? spasm_malloc(nzmax * sizeof(spasm_GFp)) : NULL;
	return A;

}





/* in-place */
void spasm_deduplicate(spasm * A) {
	int m = A->m;
	int n = A->n;
	int *Ap = A->p;
	int *Aj = A->j;
	spasm_GFp *Ax = A->x;
	int prime = A->prime;

	int *v = spasm_malloc(m * sizeof(*v));
	for (int j = 0; j < m; j++)
		v[j] = -1;

	int nz = 0;
	for (int i = 0; i < n; i++) {
		int p = nz;
		for (int it = Ap[i]; it < Ap[i + 1]; it++) {
			int j = Aj[it];
			if (v[j] < p) { /* occurs in previous row */
				v[j] = nz;
				Aj[nz] = j;
				if (Ax)
					Ax[nz] = Ax[it];
				nz++;
			} else {
				if (Ax)
					Ax[v[j]] = (Ax[v[j]] + Ax[it]) % prime;
			}
		}
		Ap[i] = p;
	}
	Ap[n] = nz;
	free(v);
	spasm_csr_realloc(A, -1);
}

/* C = compressed-row form of a triplet matrix T */
spasm *spasm_compress(const spasm_triplet * T) {
	int m = T->m;
	int n = T->n;
	int nz = T->nz;
	int *Ti = T->i;
	int *Tj = T->j;
	spasm_GFp *Tx = T->x;
	
	double start = spasm_wtime();
	fprintf(stderr, "[CSR] Compressing... ");
	fflush(stderr);

	/* allocate result */
	spasm *C = spasm_csr_alloc(n, m, nz, T->prime, Tx != NULL);

	/* get workspace */
	int *w = spasm_calloc(n, sizeof(int));
	int *Cp = C->p;
	int *Cj = C->j;
	spasm_GFp *Cx = C->x;

	/* compute row counts */
	for (int it = 0; it < nz; it++)
		w[Ti[it]]++;

	/* compute row pointers (in both Cp and w) */
	int sum = 0;
	for (int k = 0; k < n; k++) {
		Cp[k] = sum;
		sum += w[k];
		w[k] = Cp[k];
	}
	Cp[n] = sum;

	/* dispatch entries */
	for (int k = 0; k < nz; k++) {
		int px = w[Ti[k]]++;
		Cj[px] = Tj[k];
		if (Cx != NULL)
			Cx[px] = Tx[k];
	}
	free(w);
	spasm_deduplicate(C);

	/* success; free w and return C */
	char mem[16];
	int size = sizeof(int) * (n + nz) + sizeof(spasm_GFp) * ((Cx != NULL) ? nz : 0);
	spasm_human_format(size, mem);
	fprintf(stderr, "%d actual NZ, Mem usage = %sbyte [%.2fs]\n", spasm_nnz(C), mem, spasm_wtime() - start);
	return C;
}