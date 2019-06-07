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
	C->insert = 0;

	C->Lpos = -1;
	C->Kpos = -1;

	return C;
}

void class_destroy(struct class_t *C)
{
	free(C->nodes);
	free(C);
}

void print_permutation(struct permutation_factorisante_t *permutation){
	int n = permutation->n;
	struct permutation_item_t **items = permutation->items;

	for(int i=0;i<n;i++){
		struct permutation_item_t *item = items[i];
		for(int j=0; j < item->left_parentheses; j++){
			printf("(");
		}

		printf(" %d ",item->vertex);
		for(int j=0; j < item->right_parentheses; j++){
			printf(")");
		}
	}
	printf("\n");
}


void print_class(struct class_t *X)
{
	for (struct node_t * x = X->nodes->next; x != X->nodes; x = x->next)
		printf("%d ", x->vertex + 1);
}

void print_class_index(struct class_t *X)
{
	for (struct node_t * x = X->nodes->next; x != X->nodes; x = x->next)
		printf("%d ", X->begin_index);
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

void print_partition_index(struct class_t *class_head)
{
	printf("i = ");
	for (struct class_t * X = class_head->next; X != class_head;
	     X = X->next) {
		print_class_index(X);
		printf("| ");
	}
	printf("\n");
}


void class_insert(struct class_t *C, struct class_t *X)
{
	X->begin_index = C->begin_index;
	X->prev = C;
	X->next = C->next;
	X->prev->next = X;
	X->next->prev = X;
}

void class_insert_gauche(struct class_t *C, struct class_t *X)
{
	X->begin_index = C->begin_index + C->size;
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
	X->begin_index++;
}

void transfer_gauche(struct class_t *X, struct class_t *Y, struct node_t *x)
{
	if (x->class != X)
		return;
	node_remove(x);
	node_insert_gauche(x, Y);
	Y->begin_index--;
}

void refine(struct module_ctx_t *ctx, struct node_t *x, int ind_pivot)
{
	int *Ap = ctx->A->p;
	int *Aj = ctx->A->j;
	int Nx_start = Ap[x->vertex];
	int Nx_end = Ap[x->vertex + 1];
	struct class_t *X = x->class;
	struct node_t *nodes = ctx->nodes;

	 printf("Refine with x = %d / N(x) = ", x->vertex + 1);
	 for (int px = Nx_start; px < Nx_end; px++)
	      printf("%d ", Aj[px] + 1);
	 printf("\n");
	for (int it = Nx_start; it < Nx_end; it++) {
		struct node_t *y = &nodes[Aj[it]];
		y->class->marks++;
	}

	for (int it = Nx_start; it < Nx_end; it++) {
		struct node_t *y = nodes + Aj[it];
		struct class_t *Y = y->class;
		if (Y != X && Y->marks > 0 && Y->marks < Y->size) {
			//printf(" is properly split by %d\n", y->vertex);
			//printf(" x_ind : %d v_ind : %d y_ind : %d\n", X->begin_index,ind_pivot,Y->begin_index);
			struct class_t *Ya = class_new();
			if(Y->insert == 0){
				if((X->begin_index <= ind_pivot &&
					ind_pivot <= Y->begin_index)
				|| (Y->begin_index <= ind_pivot &&
					ind_pivot <= X->begin_index)){

					//printf("insertion a droite\n");
					class_insert(Y,Ya);
					Y->insert = -1;
				}
				else{
					//printf("insertion a gauche\n");
					class_insert_gauche(Y, Ya);
					Y->insert = 1;
				}
			}
			Y->split = 1;
			
		}
		Y->marks = 0;
	}

	for (int it = Nx_start; it < Nx_end; it++) {
		struct node_t *y = &nodes[Aj[it]];
		struct class_t *Y = y->class;
		if (!Y->split)
			continue;

		if(Y->insert == -1){
			//printf("transfert a droite\n");
			struct class_t *Ya = Y->next;
			transfer(Y, Ya, y);
		}
		else if(Y->insert == 1){
			//printf("transfert a gauche\n");
			struct class_t *Ya = Y->prev;
			transfer_gauche(Y, Ya, y);
		}
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

spasm *create_matrix(struct node_t *nodes, int *new_vertex, int *old_vertex, spasm *A, int n){
	int *Ap = A->p;
	int *Aj = A->j;
	int p_borne = 0;

	spasm *result = spasm_malloc(sizeof(*result));
	int *p = spasm_calloc(n+1,sizeof(*p));
	int *j = spasm_calloc(n*n,sizeof(*j));

	p[p_borne++] = 0;
	int nz = 0;
	for(int i=0;i<n;i++){
		struct node_t *node = &nodes[i];
		struct class_t *class_node = node->class;
		int n_p = 0;
		printf("x: %d class : %d\n",old_vertex[i],class_node->begin_index);
		for(int it = Ap[old_vertex[i]]; it < Ap[old_vertex[i] + 1];it++){

			struct node_t *neigh = &nodes[Aj[it]];
			printf("v: %d class : %d\n",neigh->vertex,neigh->class->begin_index);
			if(neigh->class->begin_index == node->class->begin_index){
				j[nz++] = neigh->vertex;
				n_p++;
			}
		}
		p[p_borne + 1] = p[p_borne] + n_p;
		p_borne++;
	}

	result->prime = A->prime;
	result->p = p;
	result->j = j;
	result->n = n;
	result->m = n;
	result->nzmax = nz;
	
	return result;
}

/*struct permutation_factorisante_t *make_permutation(struct class_t *class_x, struct node_t *nodes, spasm *A){
	struct permutation_factorisante_t *permutation = malloc(sizeof(*permutation));
	permutation->n = class_x->size;
	permutation->index = 0;
	struct permutation_item_t **items = calloc(permutation->n,sizeof(**items));
	permutation->items = items;
	
	struct class_t *initial_class = class_x->next;

	initial_class = class_decomposition(initial_class,nodes,A);

	for (struct class_t * X = class_x->next; X!= class_x; X = X->next) {
		if(X->size == 1){
			struct permutation_item_t *item = malloc(sizeof(*item));
			struct node_t *node = X->nodes->next;

			item->vertex = node->vertex;
			item->left_parentheses = 0;
			item->right_parentheses = 0;
			items[permutation->index] = item;
			permutation->index++;
			
		}
		else{

			//Ajouter ici la création de classe de matrice et des noeuds temporaire quand elle sera fonctionnel

			make_permutation_rec(X_tmp,tmp_nodes,B,permutation);
		}
	}
	
	return permutation;
}

void *make_permutation_rec(struct class_t *class_x, struct node_t *nodes, spasm *A, struct permutation_factorisante_t *permutation){
	struct class_t *initial_class = class_x;

	initial_class = class_decomposition(initial_class,nodes,A);

	for (struct class_t * X = class_x->next; X!= class_x; X = X->next) {
		if(X->size == 1){
			struct permutation_item_t *item = malloc(sizeof(*item));
			struct node_t *node = X->nodes->next;

			item->vertex = node->vertex;
			item->left_parentheses = 0;
			item->right_parentheses = 0;
			items[permutation->index] = item;
			permutation->index++;
			
		}
		else{
			//Ajouter ici la création de classe de matrice et des noeuds temporaire quand elle sera fonctionnel
			
			make_permutation_rec(X_tmp,tmp_nodes,B,permutation);
		}
	}
}*/

struct permutation_factorisante_t *search_fractures(struct permutation_factorisante_t *permutation, spasm *A){
	int *index_vertex = calloc(permutation->n,sizeof(*index_vertex));
	struct permutation_item_t **items = permutation->items;
	for(int i=0; i < permutation->n; i++){
		struct permutation_item_t *item = items[i];
		index_vertex[item->vertex] = i;
	}

	int *fractures = calloc(A->n,sizeof(*fractures));
	for(int i=0; i+1 < permutation->n; i++){
		struct permutation_item_t *a = items[i];
		struct permutation_item_t *b = items[i+1];
		printf("i : %d\n",a->vertex);
		printf("i+1 : %d\n",b->vertex);

		int *Ap = A->p;
		int ind_a_ap = Ap[a->vertex];
		int ind_b_ap = Ap[b->vertex];

		int nb_neigh_a = Ap[a->vertex + 1] - Ap[a->vertex];
		int nb_neigh_b = Ap[b->vertex + 1] - Ap[b->vertex];

		int *Aj = A->j;


		for(int x=0; x < A->n; x++)
			fractures[x] = 0;

		printf("neigh for a : ");
		for(int x=0; x < nb_neigh_a; x++){
			int neigh = Aj[ind_a_ap + x];
			if(neigh != b->vertex){
				printf("%d, ",neigh);
				fractures[neigh]++;
			}
		}
		printf("\n");
		printf("neigh for b : ");
		for(int x=0; x < nb_neigh_b; x++){
			int neigh = Aj[ind_b_ap + x];
			if(neigh != a->vertex){
				printf("%d, ",neigh);
				fractures[neigh]--;
			}
		}
		printf("\n");
		int x = 0;
		int flag_frac_g = 0;
		while(flag_frac_g == 0 && x < i){
			if(fractures[x] != 0){
				a->right_parentheses++;
				int ind_item = index_vertex[x];
				struct permutation_item_t *item = items[ind_item];
				item->left_parentheses++;
				flag_frac_g++;
			}
		x++;
		}
		int more_right = -1;
		for(x = i+2; x <= A->n; x++){
			if(fractures[x] != 0){
				more_right = x;
			}
		}
		if(more_right != -1){
			b->left_parentheses++;
			int ind_item = index_vertex[x];
			struct permutation_item_t *item = items[ind_item];
			item->right_parentheses++;
		}

printf("[");
for(int x=0; x < A->n; x++)
	printf("%d, ",fractures[x]);
printf("]\n");
print_permutation(permutation);
	}
	free(fractures);
	free(index_vertex);

	return permutation;

}

struct class_t *class_decomposition(struct class_t *initial_class, struct node_t *nodes, spasm * A){
	int *Ap = A->p;
	int *Aj = A->j;
	int n = A->n;
	int *queue = spasm_malloc(n * sizeof(int));
	int *mark = spasm_calloc(n, sizeof(int));
	int lo = 0, hi = 0;


	struct node_t *pivot = &nodes[0];
	int u = pivot->vertex;



	struct class_t *class_d = class_new();
	class_insert(initial_class, class_d);
	for (int it = Ap[u]; it < Ap[u + 1]; it++) {
		int v = Aj[it];
		struct node_t *s_d = &nodes[v];
		transfer(initial_class, class_d, s_d);
	}

	struct class_t *class = class_new();
	class_insert(initial_class, class);
	transfer(initial_class, class, pivot);


	/* Partie de code du départ (ne fait aucune separation dans le cas d'un matrice 	   sans arretes */

	/*for (int i = 0; i < n; i++) {

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
	}*/

	free(queue);
	free(mark);

	if (initial_class->size == 0) {
		class_remove(initial_class);
		class_destroy(initial_class);
	}

	print_partition(initial_class);

	struct module_ctx_t ctx;
	ctx.A = A;
	ctx.nodes = nodes;
	ctx.L = spasm_malloc(n * sizeof(struct class_t *));
	ctx.K = spasm_malloc(n * sizeof(struct class_t *));
	ctx.L_sp = 0;
	ctx.K_lo = 0;
	ctx.K_hi = 0;
	struct class_t *Z = initial_class->next;
	for (struct class_t * X = Z; X != initial_class; X = X->next)
		if (X->size > Z->size) {
			ctx.L[ctx.L_sp++] = Z;
			Z = X;
		} else {
			ctx.L[ctx.L_sp++] = X;
		}
	ctx.K[ctx.K_hi++] = Z;

	int ind_pivot = pivot->class->begin_index;
	while (ctx.L_sp > 0 || ctx.K_lo < ctx.K_hi) {
		if (ctx.L_sp == 0) {
			struct class_t *X = ctx.K[ctx.K_lo++];
			X->Kpos = -1;

			/*printf("Dequeuing from K : ");
			print_class(X);
			printf("\n");*/
			struct node_t *x = X->nodes->next;
			refine(&ctx, x, ind_pivot);
		} else {
			struct class_t *X = ctx.L[--ctx.L_sp];
			X->Lpos = -1;

			/*printf("Popped from L : ");
			print_class(X);
			printf("\n");*/
			for (struct node_t * x = X->nodes->next; x != X->nodes;
			     x = x->next)
				refine(&ctx, x, ind_pivot);
		}
	}
	free(ctx.L);
	free(ctx.K);

	print_partition(initial_class);
	
	return initial_class;
}


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
	class_head->begin_index = 0;
	// class_head->size = -1;
	struct class_t *initial_class = class_new();
	class_insert(class_head, initial_class);

	nodes = spasm_malloc(n * sizeof(*nodes));
	for (int i = 0; i < n; i++) {
		nodes[i].vertex = i;
		node_insert(&nodes[i], initial_class);
	}	

	print_partition(class_head);	
	initial_class = class_decomposition(initial_class,nodes,A);
	struct class_t * X_rec = class_head->next;
		
			//Insertion des noeud dans la nouvelle classe et création de la matrice correspondante

			struct class_t *tmp_class_head = class_new();
			tmp_class_head->next = tmp_class_head;
			tmp_class_head->prev = tmp_class_head;
			tmp_class_head->begin_index = 0;

			struct class_t *tmp_class = class_new();
			class_insert(tmp_class_head, tmp_class);
	
			struct node_t *nodes_tmp = spasm_malloc(X_rec->size * sizeof(*nodes_tmp));

			spasm *B = spasm_malloc(sizeof(*B));
			int *p = spasm_calloc(X_rec->size + 1, sizeof(*p));
			int *j = spasm_calloc(n*n,sizeof(*j));

			int nz = 0;
			int p_borne = 0;
			p[p_borne] = 0;
			int *new_vertex = calloc(n, sizeof(*new_vertex));

			for (int i = 0; i < n; i++){
				new_vertex[i] = -1;
			}

			int ind_new_node = 0;
			for (struct node_t * x = X_rec->nodes->next; x != X_rec->nodes; x = x->next) {
				nodes_tmp[ind_new_node].vertex = ind_new_node;
				nodes_tmp[ind_new_node].real_vertex_value = x->vertex;
				node_insert(&nodes_tmp[ind_new_node], tmp_class);
				ind_new_node++;
			}

			int ind_new_vertex = 0;
			for (struct node_t * x = X_rec->nodes->next; x != X_rec->nodes; x = x->next) {
				new_vertex[x->vertex] = nodes_tmp[ind_new_vertex++].vertex;
			}

			for (struct node_t * x = X_rec->nodes->next; x != X_rec->nodes; x = x->next) {

				struct class_t *class_node = x->class;
				int n_p = 0;
				for(int it = Ap[x->vertex]; it < Ap[x->vertex + 1];it++){

					struct node_t *neigh = &nodes[Aj[it]];
					if(neigh->class->begin_index == x->class->begin_index){
						j[nz++] = new_vertex[neigh->vertex];
						n_p++;
					}
				}
				p[p_borne + 1] = p[p_borne] + n_p;
				p_borne++;
			}

			B->prime = A->prime;
			B->p = p;
			B->j = j;
			B->n = X_rec->size;
			B->m = X_rec->size;
			B->nzmax = nz;

			printf("Bn : %d\n",B->n);
			printf("bnz : %d\n", B->nzmax);
	
			int *Bp = B->p;
			for(int i = 0; i<=B->n;i++){
				printf("%d  ",Bp[i]);
			}
			printf("\n");

			int *Bj = B->j;
			for(int i = 0; i<B->nzmax;i++){
				printf("%d  ",Bj[i]);
			}
			printf("\n");
			print_partition(tmp_class_head);

			for (struct node_t * x = tmp_class->nodes->next; x != tmp_class->nodes; x = x->next){
				printf("node vertex : %d\n ", x->vertex + 1);
				printf("real vertex : %d\n",x->real_vertex_value + 1);
			}

			printf("\n");

			for (struct node_t * x = tmp_class_head->next->nodes->next; x != tmp_class->nodes; x = x->next) {
				printf("tmp node vertex : %d\n ", x->vertex + 1);
			}

			for (int i = 0; i < n; i++)
				printf("%d ",new_vertex[i]);
			tmp_class = class_decomposition(tmp_class,tmp_class->nodes,B);
			//print_partition(class_head);


	int m = 0;
	int *module = spasm_malloc(sizeof(int) * n);
	// debugging purposes
	for (int i = 0; i < n; i++)
		module[i] = -1;
	spasm_triplet *M = spasm_triplet_alloc(0, 0, n, -1, 0);
	for (struct class_t * X = class_head->next; X != class_head;
	     X = X->next) {
		if(X->size != 0) {
			for (struct node_t *u = X->nodes->next; u != X->nodes;
				 u = u->next) {
				module[u->vertex] = m;
				spasm_add_entry(M, m, u->vertex, 1);
			}
			m++;
		}
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
	int *Mjtest = M_comp -> j;
	/*printf("%d\n",ind_pivot);
	int v;
	for(v = 0; v <= M_comp->n; v++){
		//printf("mp : %d\n",Mjtest[Mptest[v]]);
		if(Mjtest[Mptest[v]]==0)
			break;
	}

	for(int i = 0; i <= M_comp->n; i++){
		printf("Mp%d : %d\n",i,Mptest[i]);
	}


	clean_decomposition(M_comp,A,M_comp->n - v - 1);
	Mptest = M_comp -> p;
	for(int i = 0; i <= M_comp->n; i++){
		printf("Mpclean%d : %d\n",i,Mptest[i]);
	}*/
	R->M = M_comp;

	spasm_triplet_free(M);
	return R;

}
