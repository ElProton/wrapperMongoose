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

void print_tree_node_info(struct tree_node_t *node){
	if(node->father == NULL){
		printf("{root  ");
	}else{
		printf("{node  ");
	}

	printf(" n_c : %d, n_d : %d  ",node->n_class,node->n_vertex);
	printf("old_vertex : [");
	int *old_vertex = node->old_vertex;
	for(int i=0; i<node->n_vertex;i++){
		printf("%d, ",old_vertex[i]);
	}
	printf("]\n");
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

	 /*printf("Refine with x = %d / N(x) = ", x->vertex + 1);
	 for (int px = Nx_start; px < Nx_end; px++)
	      printf("%d ", Aj[px] + 1);
	 printf("\n");*/
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

struct tree_node_t *make_root(struct class_t *initial_class){
	struct tree_node_t *root = malloc(sizeof(*root));

	int n_class = 0;
	int n_vertex = 0;
	
	for(struct class_t *class = initial_class; class->check != 1; class = class->next){
		n_class++;
		class->check = 1;
		n_vertex = n_vertex + class->size;
	}
	root->n_class = n_class;

	int *old_vertex = calloc(n_vertex, sizeof(*old_vertex));
	struct class_t **classes = calloc(n_class, sizeof(**classes));

	int ind_class = 0;

	for(struct class_t *class = initial_class; class->check != 0; class = class->next){
		classes[ind_class++] = class;
		class->check = 0;
	}

	for(int i=0; i < n_vertex; i++)
		old_vertex[i] = -1;

	root->old_vertex = old_vertex;
	root->n_vertex = n_vertex;
	root->father = NULL;

	return root;
}

struct tree_node_t **make_childs(struct tree_node_t *father){
	int ind_child = 0;
	struct tree_node_t **childs = calloc(father->n_class,sizeof(**childs));

	struct class_t **classes = father->classes;
	
	for(int i=0; i < father->n_class; i++){
		struct class_t *actual_class = classes[i];
		struct tree_node_t *child = malloc(sizeof(*child));
		child->father = father;
		
		int *old_vertex = calloc(actual_class->size, sizeof(*old_vertex));
		int ind_old_vertex = 0;

		struct node_t *first = actual_class->nodes;
		old_vertex[ind_old_vertex++] = first->vertex;

		for(struct node_t *x = first->next; x->vertex != first->vertex; x = x->next){
			old_vertex[ind_old_vertex++] = x->vertex;
		}
		
		child->old_vertex = old_vertex;
		child->n_vertex = actual_class->size;
		//childs[ind_child++] = child;
	}

	return childs;
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

	struct tree_node_t *root = make_root(initial_class);
	root->matrix = A;

	struct tree_node_t **childs = make_childs(root);

	struct class_t **child_classes = root->classes;
	struct class_t *effective_class = child_classes[0];
	//print_class(effective_class);

	//Insertion des noeud dans la nouvelle classe et création de la matrice correspondante

	/*struct class_t *tmp_class_head = class_new();
	tmp_class_head->next = tmp_class_head;
	tmp_class_head->prev = tmp_class_head;
	tmp_class_head->begin_index = 0;

	struct class_t *tmp_class = class_new();
	class_insert(tmp_class_head, tmp_class);
	
	struct node_t *nodes_tmp = spasm_malloc(n * sizeof(*nodes_tmp));

	spasm *B = spasm_malloc(sizeof(*B));
	int *p = spasm_calloc(n+1,sizeof(*p));
	int *j = spasm_calloc(n*n,sizeof(*j));

	int nz = 0;
	int p_borne = 0;
	p[p_borne] = 0;

	for (int i = n - effective_class->size; i < n; i++) {

		struct node_t *node = &nodes[i];
		nodes_tmp[i].vertex = nodes[i].vertex;
		node_insert(&nodes_tmp[i], tmp_class);

		struct class_t *class_node = node->class;
		int n_p = 0;
		for(int it = Ap[node->vertex]; it < Ap[node->vertex + 1];it++){

			struct node_t *neigh = &nodes[Aj[it]];
			if(neigh->class->begin_index == node->class->begin_index){
				j[nz++] = neigh->vertex;
				n_p++;
			}
		}
		p[p_borne + 1] = p[p_borne] + n_p;
		p_borne++;
	}

	B->prime = A->prime;
	B->p = p;
	B->j = j;
	B->n = effective_class->size;
	B->m = effective_class->size;
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

	struct node_t *node = &nodes_tmp[0];*/
	//printf("%d\n",node->vertex);
	//printf("%d\n",node->next->vertex);

	//tmp_class = class_decomposition(tmp_class,tmp_class->nodes,B);
	//print_partition(class_head);	
	//print_class(tmp_class);


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
