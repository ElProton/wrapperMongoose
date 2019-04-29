#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <err.h>
#include <time.h>

#include "partitionTest.h"

struct connected_component_partition_t *component_partition(spasm * A)
{
	int An = A -> n;
	double *w,*x;
	struct modular_partition_t *partition = modular_partition(A);
	spasm *Q = partition->Q;

	struct connected_component_partition_t *result = spasm_malloc(sizeof(*result));

	if(Q->n == An){
		result->w = NULL;
		result->x = NULL;
		result->A = A;
	}else{
		
		spasm *M = partition -> M;
		int Mn = M->n;
		int *Mp = M->p;
		int *Mj = M->j;
		//int *Mj = M->j;
		int *Qj = Q->j;
		int *Qp = Q->p;
		int n_edges = spasm_nnz(Q);
		w = spasm_calloc(Mn+1, sizeof(*w));
		x = spasm_calloc(n_edges+1, sizeof(*x));
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
				double neighbors_weight = w[Qj[limit+j]];
				double module_weight = w[i];
				x[a] = neighbors_weight * module_weight;
				j++;
				a++;
			}

		limit = limit + n_neighbors;
		}

		result->A = partition -> Q;

		result->w = w;
		result->x = x;
	}
	
	return result;

}

int main(int argc, char **argv)
{
	int ch,time_sec;
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
	clock_t time = end_time - begin_time;
	begin_time = end_time;
	time_sec = time * 1000 / CLOCKS_PER_SEC;
	printf("Loading : %f \n",time_sec/1000.0);

	GraphC *g =  spasm_malloc(sizeof(*g));

	g->w = NULL;
	g->x = NULL;

	double *w,*x;
	if (mode=='B') {
        struct connected_component_t *test = spasm_malloc(sizeof(*test));



		//int components = get_components_from_complement(A);

		if(test->n > 1){
			spasm *matrix_component = spasm_calloc(test->n,sizeof(*matrix_component));
			for(int i = 0; i < test->n; i++){
				
			}
			for(int i = 0; i < test->n; i++){
				spasm *Z = matrix_component[i];
				struct connected_component_partition_t  *compo_partition = component_partition(Z);
			}
		}
		else{
			struct connected_component_partition_t *compo_partition = component_partition(A);
		}

		A = compo_partition->A;
		w = compo_partition->w;
		x = compo_partition->x;

		g->w = w;
		g->x = x;
	}


	end_time = clock();
	time = end_time - begin_time;
	begin_time = end_time;
	time_sec = time * 1000 / CLOCKS_PER_SEC;
	printf("Search Modules : %f \n", time_sec/1000.0);

	int nnz = spasm_nnz(A);
	int An = A->n;
	printf("n : %d \n", An);
	printf("nnz : %d \n", nnz);

	int n = A -> n;	
	int *Aj = A->j;
	int64_t *aj64 = spasm_calloc(spasm_nnz(A), sizeof(int64_t));
	for (int i = 0; i < spasm_nnz(A); i++)
		aj64[i] = Aj[i];


	int *Ap = A->p;
	int64_t *ap64 = spasm_calloc(n + 1, sizeof(int64_t));
	for (int i = 0; i <= n; i++)
		ap64[i] = Ap[i];
		
	g->n = n;
	g->m = n;
	g->nz = nnz;
	g->p = ap64;
	g->i = aj64;
    	

	end_time = clock();
	time = end_time - begin_time;
	begin_time = end_time;
	time_sec = time * 1000 / CLOCKS_PER_SEC;
	printf("Wrapping : %f \n", time_sec/1000.0);

	EdgeCutC* ec = connector_edge_cut(g);

	end_time = clock();
	time = end_time - begin_time;
	begin_time = end_time;
	time_sec = time * 1000 / CLOCKS_PER_SEC;
	printf("Edge Cut : %f \n", time_sec/1000.0);

	free(g);
	free(aj64);
	free(ap64);
	spasm_csr_free(A);

	/*bool *cut = ec->partition;
	for(int i=0;i<ec->n;i++){
		printf("cut %d : %d\n",i,cut[i]);
	}*/

	printf("cut cost: %f \n",  ec->cut_cost);
	printf("cut size: %li \n",  ec->cut_size);
	printf("imbalance: %f \n",  ec->imbalance);

	return EXIT_SUCCESS;
}
