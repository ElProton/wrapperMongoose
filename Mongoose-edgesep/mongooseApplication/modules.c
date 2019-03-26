#include <assert.h>
#include <stdio.h>
#include <getopt.h>
#include <err.h>

#include "modules.h"

struct class_t *class_new()
{
	struct class_t *C = spasm_malloc(sizeof(*C));
	struct node_t *head = spasm_malloc(sizeof(*head));
	C->nodes = head;
	head->next = head;
	head->prev = head;
	// head->vertex = -1;
	head->class = C;
	C->size = 0;
	C->marks = 0;
	C->split = 0;

	C->Lpos = -1;
	C->Kpos = -1;

	return C;
}

void class_destroy(struct class_t *C)
{
	free(C->nodes);
	free(C);
}

void print_class(struct class_t *X)
{
	for (struct node_t * x = X->nodes->next; x != X->nodes; x = x->next)
		printf("%d ", x->vertex + 1);
}

void print_partition(struct class_t *class_head)
{
	printf("P = ");
	for (struct class_t * X = class_head->next; X != class_head;
	     X = X->next) {
		print_class(X);
		printf("| ");
	}
	printf("\n");
}

void class_insert(struct class_t *C, struct class_t *X)
{
	X->prev = C;
	X->next = C->next;
	X->prev->next = X;
	X->next->prev = X;
}

void class_insert_gauche(struct class_t *C, struct class_t *X)
{
	X->next = C;
	X->prev = C->prev;
	X->prev->next = X;
	X->next->prev = X;
}

void class_remove(struct class_t *X)
{
	X->next->prev = X->prev;
	X->prev->next = X->next;
}

void node_insert(struct node_t *x, struct class_t *C)
{
	x->class = C;
	C->size++;
	struct node_t *head = C->nodes;
	x->prev = head;
	x->next = head->next;
	x->next->prev = x;
	x->prev->next = x;
}

void node_insert_gauche(struct node_t *x, struct class_t *C)
{
	x->class = C;
	C->size++;
	struct node_t *head = C->nodes;
	x->next = head;
	x->prev = head->prev;
	x->next->prev = x;
	x->prev->next = x;
}

void node_remove(struct node_t *x)
{
	x->class->size--;
	x->next->prev = x->prev;
	x->prev->next = x->next;
}

void transfer(struct class_t *X, struct class_t *Y, struct node_t *x)
{
	if (x->class != X)
		return;
	node_remove(x);
	node_insert(x, Y);
}

void transfer_gauche(struct class_t *X, struct class_t *Y, struct node_t *x)
{
	if (x->class != X)
		return;
	node_remove(x);
	node_insert_gauche(x, Y);
}

void refine(struct module_ctx_t *ctx, struct node_t *x)
{
	int *Ap = ctx->A->p;
	int *Aj = ctx->A->j;
	int Nx_start = Ap[x->vertex];
	int Nx_end = Ap[x->vertex + 1];
	struct class_t *X = x->class;
	struct node_t *nodes = ctx->nodes;

	// printf("Refine with x = %d / N(x) = ", x->vertex + 1);
	// for (int px = Nx_start; px < Nx_end; px++)
	//      printf("%d ", Aj[px] + 1);
	// printf("\n");
	for (int it = Nx_start; it < Nx_end; it++) {
		struct node_t *y = &nodes[Aj[it]];
		y->class->marks++;
	}

	for (int it = Nx_start; it < Nx_end; it++) {
		struct node_t *y = nodes + Aj[it];
		struct class_t *Y = y->class;
		if (Y != X && Y->marks > 0 && Y->marks < Y->size) {
			//printf(" is properly split by %d\n", y->vertex);
			// printf(" is properly split by %d\n", x->vertex);
			struct class_t *Ya = class_new();
			class_insert(Y, Ya);
			Y->split = 1;
			
		}
		Y->marks = 0;
	}

	for (int it = Nx_start; it < Nx_end; it++) {
		struct node_t *y = &nodes[Aj[it]];
		struct class_t *Y = y->class;
		if (!Y->split)
			continue;

		struct class_t *Ya = Y->next;
		//printf("s : %d\n" , y->vertex);
		transfer(Y, Ya, y);
	}

	for (int it = Nx_start; it < Nx_end; it++) {

		struct node_t *y = nodes + Aj[it];
		struct class_t *Y = y->class->prev;
		if (Y->split) {
			struct class_t *Ya = Y->next;
			Y->split = 0;
			struct class_t **L = ctx->L;
			struct class_t **K = ctx->K;
			if (Y->Lpos >= 0) {
				Ya->Lpos = ctx->L_sp;
				L[ctx->L_sp++] = Ya;
				// printf("Split class was already in L. Adding to L: ");
				// print_class(Ya);
				// printf("\n");
				continue;
			}

			struct class_t *Ymin, *Ymax;
			if (Y->size < Ya->size) {
				Ymin = Y;
				Ymax = Ya;
			} else {
				Ymin = Ya;
				Ymax = Y;
			}
			// printf("Ymin = ");
			// print_class(Ymin);
			// printf("\nYmax = ");
			// print_class(Ymax);
			// printf("\n");
			// printf("Adding Ymin to L\n");
			Ymin->Lpos = ctx->L_sp;
			L[ctx->L_sp++] = Ymin;

			if (Y->Kpos >= 0) {
				if (Y != Ymax) {
					Ymax->Kpos = Y->Kpos;
					K[Y->Kpos] = Ymax;
					Y->Kpos = -1;
					// printf("Split class was in K, and the Ymax part replaces the smallest part in K\n");
				} else {
					// printf("Split class was in K, and the largest part remains in K\n");
				}
			} else {
				// printf("Split class was not in K. Adding Ymax to K\n");
				Ymax->Kpos = ctx->K_hi;
				K[ctx->K_hi++] = Ymax;
			}

		}
	}

}

int fraction(int *potentiel,spasm *M, spasm *A, int n, int i, int j){
	int *Ap = A->p;
	int *Aj = A->j;
	int *Mj = M->j;
	int i_index = Mj[potentiel[i]]; //L'indice dans Ap du sommet à l'indice i dans Mp
	int j_index = Mj[potentiel[j]];

	for(int a = 0; a < i; a++){
		int a_index = Mj[potentiel[a]];//le no de sommet dans A

		if(a_index != -1){
			int a_ind = Ap[a_index];
			int edge_i = 0;
			int edge_j = 0;
			for(int p = 0; p < Ap[a_index +1] - a_ind; p++) {
				if(Aj[a_ind+p] == i_index){
					edge_i++;
				}
				if(Aj[a_ind+p] == j_index){
					edge_j++;
				}
			}
			if(edge_i != edge_j){
				//printf("casseur %d entre %d et %d\n",a_index,i,j);
				goto casseur;
			}
		}
	}

	for(int a = j+1; a < n; a++){
		int a_index = Mj[potentiel[a]];//le no de sommet dans A

		if(a_index != -1){
			int a_ind = Ap[a_index];
			int edge_i = 0;
			int edge_j = 0;
			for(int p = 0; p < Ap[a_index +1] - a_ind; p++) {
				if(Aj[a_ind+p] == i_index){
					edge_i++;
				}
				if(Aj[a_ind+p] == j_index){
					edge_j++;
				}
			}
			if(edge_i != edge_j){
				//printf("casseur %d entre %d et %d\n",a_index,i,j);
				goto casseur;
			}
		}
	}
	/*for(int ind = 0; ind < n; ind++){
		int s = Ap[ind];
		if(s < i_ind || s > j_ind){
			int nb_voisin = Ap[ind+1] - s;
			int edge_i = 0;
			int edge_j = 0;
			for(int p = 0; p < nb_voisin; p++){
				if(Aj[s+p] == i){
					edge_i++;
				}
				if(Aj[s+p] == j){
					edge_j++;
				}
			}
			if(edge_i != edge_j){
				goto casseur;
			}
		}
	}*/
	return 0;

casseur:
	return 1;
}

void clean_decomposition(spasm *M, spasm *A, int ind_v){
	int n = M->n;
	int *Mp = M->p;
	int *Mj = M->j;
	int *potentiel = spasm_calloc(n,sizeof(int));

	int del = 0;
	for(int i=0; i<n; i++){
		potentiel[i] = Mp[i];
	}
	int b_d = n;
	int b_g = ind_v;
	for(int i=b_g; i > 0; i--){
		//printf("g : %d d : %d\n",i,b_d);
		for(int j=b_d-1; j > i; j--){
			//printf("i: %d et j : %d\n",i,j);
			if(!fraction(potentiel,M,A,n,i,j)){

				for(int k = i+1; k<=j; k++){
					//printf("k: %d \n",k);
					potentiel[k] = -1;
					del++;
				}
				b_g = i;
				if(j == b_d){
					b_d = i;
					potentiel[b_d] = -1;
					del++;
				}
			
				break;
			}
		}
	}

	int *Mptest = calloc(M->n,sizeof(int));
	int z = 0;
	for(int i=0;i<M->n;i++){
		if(potentiel[i] != -1){
			Mptest[z] = potentiel[i];
			z++;
		}
	}
	Mptest[n - del] = Mp[M->n];

	M->p = Mptest;
	M->n = n - del;
	//printf("%d\n",del);
}

/*void clean_decomposition(spasm *M,spasm *A){
	int n = M->n;
	int *Mp = M->p;
	int *Mj = M->j;
	int *res = spasm_calloc(n+1,sizeof(int));
	int *potentiel = spasm_calloc(2,sizeof(int));
	int b_p = 1;
	int c = 0;
	int first = 1;
	for(int ind = n; ind > 0; ind--){
		int i = Mp[ind];
		int a_index = Mj[i-1];
		potentiel[b_p--] = a_index;
		if(b_p < 0){
			if(fraction(A,potentiel[0],potentiel[1])){
				if(first==1){
					res[c++] = Mp[ind+1];
					first = 0;
				}
				res[c++] = i;
			}
			potentiel[1] = potentiel[0];
			b_p = 0;
		}

	}
	free(potentiel);
	res[c] = 0;
	int *Mptest = calloc(c+1,sizeof(int));
	for(int i = 0; i <= c; i++){
		Mptest[c-i] = res[i];
	}
	M->p = Mptest;
	M->n = c;
	free(res);
}*/

struct modular_partition_t *modular_partition(spasm * A)
{
	struct class_t *class_head;
	struct node_t *nodes;
	int *Ap = A->p;
	int *Aj = A->j;
	int n = A->n;

	class_head = class_new();
	class_head->next = class_head;
	class_head->prev = class_head;
	// class_head->size = -1;
	struct class_t *initial_class = class_new();
	class_insert(class_head, initial_class);

	nodes = spasm_malloc(n * sizeof(*nodes));
	for (int i = 0; i < n; i++) {
		nodes[i].vertex = i;
		node_insert(&nodes[i], initial_class);
	}

	int *queue = spasm_malloc(n * sizeof(int));
	int *mark = spasm_calloc(n, sizeof(int));
	int lo = 0, hi = 0;
	//print_partition(class_head);


	struct node_t *pivot = &nodes[0];
	int u = pivot->vertex;
	/*struct class_t *class_d = class_new();
	class_insert(initial_class, class_d);
	for (int it = Ap[u]; it < Ap[u + 1]; it++) {
		int v = Aj[it];
		struct node_t *s_d = &nodes[v];
		transfer(initial_class, class_d, s_d);
	}

	struct class_t *class = class_new();
	class_insert(initial_class, class);
	transfer(initial_class, class, pivot);*/
	for (int i = 0; i < n; i++) {

		if (mark[i])
			continue;
		int start = hi;
		queue[hi++] = i;
		mark[i] = 1;;
		while (lo < hi) {
			int u = queue[lo++];
			for (int it = Ap[u]; it < Ap[u + 1]; it++) {
				int v = Aj[it];
				if (mark[v])
					continue;
				queue[hi++] = v;
				mark[v] = 1;
			}
		}
		if(hi-start<=2)
			continue;
		struct node_t *x = &nodes[i];
		struct class_t *class = class_new();
		class_insert(initial_class,class);
		transfer(initial_class,class,x);
	}

	free(queue);
	free(mark);

	if (initial_class->size == 0) {
		class_remove(initial_class);
		class_destroy(initial_class);
	}

	struct module_ctx_t ctx;
	ctx.A = A;
	ctx.nodes = nodes;
	ctx.L = spasm_malloc(n * sizeof(struct class_t *));
	ctx.K = spasm_malloc(n * sizeof(struct class_t *));
	ctx.L_sp = 0;
	ctx.K_lo = 0;
	ctx.K_hi = 0;
	struct class_t *Z = class_head->next;
	for (struct class_t * X = Z; X != class_head; X = X->next)
		if (X->size > Z->size) {
			ctx.L[ctx.L_sp++] = Z;
			Z = X;
		} else {
			ctx.L[ctx.L_sp++] = X;
		}
	ctx.K[ctx.K_hi++] = Z;

	while (ctx.L_sp > 0 || ctx.K_lo < ctx.K_hi) {
		//print_partition(class_head);
		if (ctx.L_sp == 0) {
			struct class_t *X = ctx.K[ctx.K_lo++];
			X->Kpos = -1;

			/*printf("Dequeuing from K : ");
			print_class(X);
			printf("\n");*/
			struct node_t *x = X->nodes->next;
			refine(&ctx, x);
		} else {
			struct class_t *X = ctx.L[--ctx.L_sp];
			X->Lpos = -1;

			/*printf("Popped from L : ");
			print_class(X);
			printf("\n");*/
			for (struct node_t * x = X->nodes->next; x != X->nodes;
			     x = x->next)
				refine(&ctx, x);
		}
	}
	free(ctx.L);
	free(ctx.K);

	print_partition(class_head);
	int m = 0;
	int *module = spasm_malloc(sizeof(int) * n);
	// debugging purposes
	for (int i = 0; i < n; i++)
		module[i] = -1;
	spasm_triplet *M = spasm_triplet_alloc(0, 0, n, -1, 0);
	for (struct class_t * X = class_head->next; X != class_head;
	     X = X->next) {
		for (struct node_t * u = X->nodes->next; u != X->nodes;
		     u = u->next) {
			module[u->vertex] = m;
			spasm_add_entry(M, m, u->vertex, 1);
		}
		m++;
	}

	spasm_triplet *S = spasm_triplet_alloc(n, n, spasm_nnz(A), -1, 0);
	spasm_triplet *Q = spasm_triplet_alloc(m, m, n + m, -1, 0);
	int *module_reached = spasm_malloc(m * sizeof(int));
	for (int i = 0; i < m; i++)
		module_reached[i] = -1;
	for (struct class_t * X = class_head->next; X != class_head;
	     X = X->next) {
		int first = 1;
		for (struct node_t * x = X->nodes->next; x != X->nodes;
		     x = x->next) {
			int u = x->vertex;
			assert(module[u] >= 0);
			for (int it = Ap[u]; it < Ap[u + 1]; it++) {
				int v = Aj[it];
				assert(module[v] >= 0);
				int Mu = module[u];
				int Mv = module[v];
				if (Mu == Mv) {
					// fprintf(stderr, "%d -> %d, inside module\n", u + 1, v + 1);
					spasm_add_entry(S, u, v, 1);
				} else {
					if (first && module_reached[Mv] != Mu) {
						// fprintf(stderr, "%d -> %d, from module %d to other module %d\n", u + 1, v + 1, Mu, Mv);
						spasm_add_entry(Q, Mu, Mv, 1);
						module_reached[Mv] = Mu;
					}
				}
			}
			first = 0;
		}
	}
	free(module_reached);
	free(module);

	struct class_t *X = class_head->next;
	while (X != class_head) {
		struct class_t *Y = X->next;
		class_destroy(X);
		X = Y;
	}
	class_destroy(class_head);
	free(nodes);

	struct modular_partition_t *R = spasm_malloc(sizeof(*R));
	R->S = spasm_compress(S);
	spasm_triplet_free(S);
	R->Q = spasm_compress(Q);
	spasm_triplet_free(Q);
	spasm *M_comp = spasm_malloc(sizeof(M_comp));
	M_comp = spasm_compress(M);
	int *Mptest = M_comp -> p;

	/*int v;
	for(v = 0; v <= M_comp->n; v++){
		if(Mptest[v]==1)
			break;
	}

	for(int i = 0; i <= M_comp->n; i++){
		printf("Mp%d : %d\n",i,Mptest[i]);
	}


	clean_decomposition(M_comp,A,v);
	Mptest = M_comp -> p;
	for(int i = 0; i <= M_comp->n; i++){
		printf("Mpclean%d : %d\n",i,Mptest[i]);
	}*/
	R->M = M_comp;

	spasm_triplet_free(M);
	return R;

}
