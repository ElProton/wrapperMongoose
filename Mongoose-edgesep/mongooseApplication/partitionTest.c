#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <err.h>
#include <time.h>


#include "Mongoose_EdgeCut_Connector.h"
#include "mini_spasm.h"

#include "partitionTest.h"

struct component_partition_t *component_partition(spasm *A){
		spasm *B = A;

		double *w,*x;
		int n_module = A -> n;
		struct modular_partition_t *partition = modular_partition(A);
		struct component_partition_t *result = malloc(sizeof(*result));
		spasm *Q = partition->Q;

		if(Q->n == n_module){
			result->w = NULL;
			result->x = NULL;
			result->A = A;
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

			B = partition -> Q;
			n_module = B -> n;

			result->A = B;
			result->w = w;
			result->x = x;
		}

		return result;
}

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
	if (mode=='B') {
		struct component_partition_t *quotient = component_partition(A);

		A = quotient -> A;
		g->w = quotient -> w;
		g->x = quotient -> x;
	}

	end_time = clock();
	clock_t search_modules_time = end_time - begin_time;
	begin_time = end_time;
	msec = search_modules_time * 1000 / CLOCKS_PER_SEC;
	printf("Search Modules : %f \n", msec/1000.0);

	int n = A -> n;
	int nnz = spasm_nnz(A);
	printf("n : %d \n", n);
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


	printf("cut cost: %f \n",  ec->cut_cost);
	printf("cut size: %li \n",  ec->cut_size);
	printf("imbalance: %f \n",  ec->imbalance);


	/*struct permutation_factorisante_t *perm = malloc(sizeof(*perm));
	perm->n = 11;
	struct permutation_item_t **items = calloc(11,sizeof(*items));
	struct permutation_item_t *item1 = malloc(sizeof(*item1));
	struct permutation_item_t *item2 = malloc(sizeof(*item1));
	struct permutation_item_t *item3 = malloc(sizeof(*item1));
	struct permutation_item_t *item4 = malloc(sizeof(*item1));
	struct permutation_item_t *item5 = malloc(sizeof(*item1));
	struct permutation_item_t *item6 = malloc(sizeof(*item1));
	struct permutation_item_t *item7 = malloc(sizeof(*item1));
	struct permutation_item_t *item8 = malloc(sizeof(*item1));
	struct permutation_item_t *item9 = malloc(sizeof(*item1));
	struct permutation_item_t *item10 = malloc(sizeof(*item1));
	struct permutation_item_t *item11 = malloc(sizeof(*item1));

	item1->vertex = 0;
	item2->vertex = 1;
	item3->vertex = 2;
	item4->vertex = 3;
	item5->vertex = 4;
	item6->vertex = 5;
	item7->vertex = 6;
	item8->vertex = 7;
	item9->vertex = 8;
	item10->vertex = 9;
	item11->vertex = 10;

			item1->left_parentheses = 0;
			item1->right_parentheses = 0;
			item2->left_parentheses = 0;
			item2->right_parentheses = 0;
			item3->left_parentheses = 0;
			item3->right_parentheses = 0;
			item4->left_parentheses = 0;
			item4->right_parentheses = 0;
			item5->left_parentheses = 0;
			item5->right_parentheses = 0;
			item6->left_parentheses = 0;
			item6->right_parentheses = 0;
			item7->left_parentheses = 0;
			item7->right_parentheses = 0;
			item8->left_parentheses = 0;
			item8->right_parentheses = 0;
			item9->left_parentheses = 0;
			item9->right_parentheses = 0;
			item10->left_parentheses = 0;
			item10->right_parentheses = 0;
			item11->left_parentheses = 0;
			item11->right_parentheses = 0;


	items[0] = item11;
	items[1] = item10;
	items[2] = item9;
	items[3] = item8;
	items[4] = item7;
	items[5] = item6;
	items[6] = item5;
	items[7] = item1;
	items[8] = item2;
	items[9] = item3;
	items[10] = item4;

	perm->items = items;

	perm = search_fractures(perm, A);
	print_permutation(perm);

	free(item1);
	free(item2);
	free(item3);
	free(item4);
	free(item5);
	free(item6);
	free(item7);
	free(item8);
	free(item9);
	free(item10);
	free(item11);
	free(items);
	free(perm);*/


	return EXIT_SUCCESS;
}
