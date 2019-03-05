#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <err.h>


#include "Mongoose_EdgeCut_Connector.h"
#include "mini_spasm.h"

#include "modules.h"

int main(int argc, char **argv)
{
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


	spasm_triplet *T = spasm_load_mm(stdin, -1);
	assert(T->n == T->m);
	
	/* élimine les entrées sur la diagonale */
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

	if (mode=='B') {
		struct modular_partition_t *partition = modular_partition(A);
		A = partition -> M;
	}

	int n = A -> n;
	int *Ap = A->p;
	int *Aj = A->j;

	int64_t *ap64 = spasm_calloc(n + 1, sizeof(int64_t));
	int64_t *aj64 = spasm_calloc(spasm_nnz(A), sizeof(int64_t));

	for (int i = 0; i < n; i++)
		ap64[i] = Ap[i];

	for (int i = 0; i < spasm_nnz(A); i++)
		aj64[i] = Aj[i];
	
	/*n = A -> m;
	double *axd = spasm_calloc(n, sizeof(double));
	for(i = 0; i <= n; i++){
		val = A->x[i];
		vald = (double) val;
		axd[i] = vald;
	}*/

	GraphC *g =  spasm_malloc(sizeof(*g));
		
	g->n = n;
	g->m = n;
	g->nz = spasm_nnz(A);
	g->p = ap64;
	g->i = aj64;
    	g->x = NULL;
	g->w = NULL;

	connector_edge_cut(g);

	free(g);
	free(aj64);
	free(ap64);
	spasm_csr_free(A);

	return 0;
}