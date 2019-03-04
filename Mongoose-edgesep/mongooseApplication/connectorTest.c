#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <err.h>


#include "Mongoose_EdgeCut_Connector.h"
#include "mini_spasm.h"

#include "modules.h"

int main(int argc, char **argv)
{
	int ch,i,px,val,n;
	int64_t val64;
	double vald;
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
	int *Ti = T->i;
	int *Tj = T->j;
	for (px = 0; px < T->nz; px++) {
		if (Ti[px] == Tj[px]) {
			spasm_swap(Ti, px, T->nz - 1);
			spasm_swap(Tj, px, T->nz - 1);
			T->nz--;
		}
	}

	spasm *A = spasm_compress(T);
	spasm_triplet_free(T);

	if(mode=='B'){
	struct modular_partition_t *partition = modular_partition(A);
	A = partition -> M;
	}

	n = A -> n;
	int64_t *ap64 = spasm_calloc(n, sizeof(int64_t));
	for(i = 0; i <= n; i++){
		val = A->p[i];
		val64 = (int64_t) val;
		ap64[i] = val64;
	}

	n = A -> m;
	int64_t *aj64 = spasm_calloc(n, sizeof(int64_t));
	for(i = 0; i <= n; i++){
		val = A->j[i];
		val64 = (int64_t) val;
		aj64[i] = val64;
	}

	/*n = A -> m;
	double *axd = spasm_calloc(n, sizeof(double));
	for(i = 0; i <= n; i++){
		val = A->x[i];
		vald = (double) val;
		axd[i] = vald;
	}*/

	GraphC *g =  spasm_malloc(sizeof(*g));
		

	g->n = A->n;
	g->m = A->m;
	g->nz = A->p[A->n];
	g->p = ap64;    /* problème de type */
	g->i = aj64;    /* problème de type */
    	g->x = NULL;     /* problème de type */
	g->w = NULL;

	//spasm_csr_free(A);
	connector_edge_cut(g);

	free(g);
	free(aj64);
	free(ap64);

	return 0;
}
