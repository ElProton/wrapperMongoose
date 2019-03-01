#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <err.h>


#include "../Include/Mongoose_EdgeCut_Connector.h"
#include "mini_spasm.h"
#include "modules.c"




int main(int argc, char **argv)
{
	int ch;
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
	for (int px = 0; px < T->nz; px++) {
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

	GraphC* g = (GraphC*) malloc(sizeof(GraphC));
	

	g->n = A->n;
	g->m = A->m;
	g->nz = A->p[A->n];
	g->p = (int64_t*) A->p;
	g->i = (int64_t*) A->j;
    	g->x = (double*) A->x;
	g->w = NULL;

	
	connector_edge_cut(g);
/*

	connector_edge_cut(argv[1]);*/
	return 0;
}
