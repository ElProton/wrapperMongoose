#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <err.h>
#include <time.h>


#include "Mongoose_EdgeCut_Connector.h"
#include "mini_spasm.h"

#include "modules.h"

int main(int argc, char **argv)
{
	int ch,msec;
	struct option longopts[6] = {
		{"A", no_argument, NULL, 'A'},
		{"B", no_argument, NULL, 'B'},
		{NULL, 0, NULL, 0}
	};

	char mode = -1;
	while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
		switch (ch) {
		case 'A':
		case 'B':
			mode = ch;
			break;
		default:
			errx(1, "Unknown option");
		}
	}
	if (mode == -1) {
		warnx("Choosing classic mode by default");
		mode = 'A';
	}

	clock_t begin_time = clock();
	spasm_triplet *T = spasm_load_mm(stdin, -1);
	assert(T->n == T->m);
	
	
	int *Ti = T->i;
	int *Tj = T->j;
	for (int px = 0; px < T->nz; px++)
		if (Ti[px] == Tj[px]) {
			spasm_swap(Ti, px, T->nz - 1);
			spasm_swap(Tj, px, T->nz - 1);
			T->nz--;
		}

	spasm *A = spasm_compress(T);
	spasm_triplet_free(T);

	clock_t loading_matrix_time = clock() - begin_time;
	msec = loading_matrix_time * 1000 / CLOCKS_PER_SEC;
	printf("Loading : %d s  %d ms \n",msec/1000,msec%1000);


	if (mode=='B') {
		struct modular_partition_t *partition = modular_partition(A);
		spasm *Q = partition->Q;


		A = partition -> Q;
	}

	clock_t search_modules_time = clock() - loading_matrix_time;
	msec = search_modules_time * 1000 / CLOCKS_PER_SEC;
	printf("Search Modules : %d s  %d ms \n", msec/1000,msec%1000);

	int n = A -> n;
	int m = A -> m;
	int *Ap = A->p;
	int *Aj = A->j;
	int *Ax = A->x;
	int nzmax = A->nzmax;

	int64_t *ap64 = spasm_calloc(n + 1, sizeof(int64_t));
	int64_t *aj64 = spasm_calloc(spasm_nnz(A), sizeof(int64_t));

	for (int i = 0; i < n; i++)
		ap64[i] = Ap[i];

	for (int i = 0; i < spasm_nnz(A); i++)
		aj64[i] = Aj[i];
	
	
	/*double *axd = spasm_calloc(spasm_nnz(A), sizeof(double));
	for(int i = 0; i < spasm_nnz(A); i++){
		axd[i] = Ax[i];
	}*/

	GraphC *g =  spasm_malloc(sizeof(*g));
		
	g->n = n;
	g->m = n;
	g->nz = spasm_nnz(A);
	g->p = ap64;
	g->i = aj64;
    	g->x = NULL;
	g->w = NULL;

	clock_t wrapping_time = clock() - search_modules_time;
	msec = wrapping_time * 1000 / CLOCKS_PER_SEC;
	printf("Wrapping : %d s  %d ms \n", msec/1000,msec%1000);

	EdgeCutC* ec = connector_edge_cut(g);

	clock_t edge_cut_time = clock() - wrapping_time;
	msec = edge_cut_time * 1000 / CLOCKS_PER_SEC;
	printf("Edge Cut : %d s  %d ms \n", msec/1000,msec%1000);

	free(g);
	free(aj64);
	free(ap64);
	//free(axd);
	spasm_csr_free(A);

	printf("cut cost: %f \n",  ec->cut_cost);
	printf("cut size: %li \n",  ec->cut_size);
	printf("imbalance: %f \n",  ec->imbalance);

	return EXIT_SUCCESS;
}
