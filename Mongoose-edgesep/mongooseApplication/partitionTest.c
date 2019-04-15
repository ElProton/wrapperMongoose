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
	for (int px = 0; px < T->nz; px++){
		if (Ti[px] == Tj[px]) {
			spasm_swap(Ti, px, T->nz - 1);
			spasm_swap(Tj, px, T->nz - 1);
			T->nz--;
		}
	}

	spasm *A = spasm_compress(T);
	spasm_triplet_free(T);

	clock_t end_time = clock();
	clock_t loading_matrix_time = end_time - begin_time;
	begin_time = end_time;
	msec = loading_matrix_time * 1000 / CLOCKS_PER_SEC;
	printf("Loading : %f \n",msec/1000.0);

	GraphC *g =  spasm_malloc(sizeof(*g));

	g->w = NULL;
	g->x = NULL;

	int64_t *ap64;
	double *w,*x;
	int n_module = A -> n;
	if (mode=='B') {
		struct modular_partition_t *partition = modular_partition(A);
		spasm *Q = partition->Q;

		if(Q->n == n_module){
			n_module = 0;
		}else{

			spasm *M = partition -> M;
			int Mn = M->n;
			int *Mp = M->p;
			int *Mj = M->j;
			int *Qj = Q->j;
			int *Qp = Q->p;
			int n_edges = spasm_nnz(Q);

			w = spasm_calloc(Mn+1, sizeof(double));
			x = spasm_calloc(n_edges+1, sizeof(double));
			for (int i = 0; i < Mn; i++){
				w[i] = Mp[i+1]-Mp[i];
				/*printf("w%d : %f  [",i,w[i]);
				for (int k=Mp[i]; k < Mp[i+1]; k++)
					printf("%d ", Mj[k]);
				printf("]\n");*/
			}
			int limit = 0;
			int a = 0;
			for (int i = 0; i < Mn; i++){
			
				int n_neighbors = Qp[i+1]-Qp[i];
				int j = 0;
				while(j < n_neighbors){
					double neigh_weight = w[Qj[limit+j]];
					double mod_weight = w[i];
					x[a] = neigh_weight*mod_weight;
					j++;
					a++;
				}
			limit = limit + n_neighbors;
			}

			A = partition -> Q;
			n_module = A -> n;

			g->w = w;
			g->x = x;
		}

	}

	end_time = clock();
	clock_t search_modules_time = end_time - begin_time;
	begin_time = end_time;
	msec = search_modules_time * 1000 / CLOCKS_PER_SEC;
	printf("Search Modules : %f \n", msec/1000.0);

	int n = A -> n;
	int nnz = spasm_nnz(A);
	printf("n : %d \n", n_module);
	printf("nnz : %d \n", nnz);
	
	int *Aj = A->j;
	int64_t *aj64 = spasm_calloc(spasm_nnz(A), sizeof(int64_t));
	for (int i = 0; i < spasm_nnz(A); i++)
		aj64[i] = Aj[i];


	int *Ap = A->p;
	ap64 = spasm_calloc(n + 1, sizeof(int64_t));
	for (int i = 0; i <= n; i++)
		ap64[i] = Ap[i];
		
	g->n = n;
	g->m = n;
	g->nz = nnz;
	g->p = ap64;
	g->i = aj64;
    	

	end_time = clock();
	clock_t wrapping_time = end_time - begin_time;
	begin_time = end_time;
	msec = wrapping_time * 1000 / CLOCKS_PER_SEC;
	printf("Wrapping : %f \n", msec/1000.0);

	EdgeCutC* ec = connector_edge_cut(g);

	end_time = clock();
	clock_t edge_cut_time = end_time - begin_time;
	begin_time = end_time;
	msec = edge_cut_time * 1000 / CLOCKS_PER_SEC;
	printf("Edge Cut : %f \n", msec/1000.0);

	free(g);
	free(aj64);
	free(ap64);
	spasm_csr_free(A);

	bool *cut = ec->partition;
	for(int i=0;i<ec->n;i++){
		printf("cut %d : %d\n",i,cut[i]);
	}

	printf("cut cost: %f \n",  ec->cut_cost);
	printf("cut size: %li \n",  ec->cut_size);
	printf("imbalance: %f \n",  ec->imbalance);

	return EXIT_SUCCESS;
}
