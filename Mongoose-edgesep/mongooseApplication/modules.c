#include <assert.h>
#include <stdio.h>
#include <getopt.h>
#include <err.h>

#include "mini_spasm.h"

// modules.c:(.text+0x12ed): undefined reference to `spasm_compress'


struct modular_partition_t {
	spasm *S, *M, *Q;
};

struct node_t;
struct class_t;

struct node_t {
	struct node_t *prev, *next;
	int vertex;
	struct class_t *class;
};

struct class_t {
	struct class_t *prev, *next;
	struct node_t *nodes;
	int size;
	int marks;
	int split;

	int Lpos;
	int Kpos;

};

struct module_ctx_t {
	spasm *A;
	struct node_t *nodes;
	struct class_t **L;
	struct class_t **K;
	int L_sp, K_lo, K_hi;

};

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
			// print_class(Y);
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
		transfer(Y, Ya, y);
	}

	for (int it = Nx_start; it < Nx_end; it++) {
		struct node_t *y = nodes + Aj[it];
		struct class_t *Y = y->class->prev;
		if (Y->split) {
			struct class_t *Ya = Y->next;
			// printf("Processing update to classes Yb = ");
			// print_class(Y);
			// printf(" and Ya = ");
			// print_class(Ya);
			// printf("\n");
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
	for (int i = 0; i < n; i++) {
		if (mark[i])
			continue;
		int start = hi;
		queue[hi++] = i;
		mark[i] = 1;
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

		if (hi - start <= 2)
			continue;
		struct node_t *x = &nodes[i];
		struct class_t *class = class_new();
		class_insert(initial_class, class);
		transfer(initial_class, class, x);

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
		// print_partition(class_head);
		if (ctx.L_sp == 0) {
			struct class_t *X = ctx.K[ctx.K_lo++];
			X->Kpos = -1;

			//printf("Dequeuing from K : ");
			//print_class(X);
			//printf("\n");
			struct node_t *x = X->nodes->next;
			refine(&ctx, x);
		} else {
			struct class_t *X = ctx.L[--ctx.L_sp];
			X->Lpos = -1;

			//printf("Popped from L : ");
			//print_class(X);
			//printf("\n");
			for (struct node_t * x = X->nodes->next; x != X->nodes;
			     x = x->next)
				refine(&ctx, x);
		}
	}
	free(ctx.L);
	free(ctx.K);

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
	R->M = spasm_compress(M);
	spasm_triplet_free(M);
	return R;

}

/**
int main(int argc, char **argv)
{
	int ch;
	struct option longopts[6] = {
		{"tabulated", no_argument, NULL, 't'},
		{"verbose", no_argument, NULL, 'v'},
		{"M", no_argument, NULL, 'M'},
		{"Q", no_argument, NULL, 'Q'},
		{"S", no_argument, NULL, 'S'},
		{NULL, 0, NULL, 0}
	};

	char mode = -1;
	while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
		switch (ch) {
		case 't':
		case 'v':
		case 'M':
		case 'Q':
		case 'S':
			if (mode != -1)
				errx(1, "incompatible options");
			mode = ch;
			break;
		default:
			errx(1, "Unknown option");
		}
	}
	if (mode == -1) {
		warnx("Choosing verbose mode by default");
		mode = 'v';
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

	struct modular_partition_t *partition = modular_partition(A);
	switch (mode) {
	case 'Q':
		spasm_save_csr(stdout, partition->Q);
		break;
	case 'M':
		spasm_save_csr(stdout, partition->M);
		break;
	case 'S':
		spasm_save_csr(stdout, partition->S);
		break;
	case 'v':
	case 't':
		{
			int trivial = 0;
			int nontrivial = 0;
			int nontrivial_size = 0;
			int largest = 0;
			int module_edges = spasm_nnz(partition->S);
			int quotient_edges = spasm_nnz(partition->Q);
			spasm *M = partition->M;
			int m = M->n;
			for (int i = 0; i < m; i++) {
				int size = spasm_row_weight(M, i);
				if (size == 1) {
					trivial++;
				} else {
					nontrivial++;
					nontrivial_size += size;
				}
				largest = spasm_max(largest, size);
			}
			// we have not lost vertices
			assert(nontrivial_size + trivial == A->n);
			// we have not lost edges
			int edges = module_edges;
			spasm *Q = partition->Q;
			int *Qp = Q->p;
			int *Qj = Q->j;
			for (int i = 0; i < m; i++)
				for (int it = Qp[i]; it < Qp[i + 1]; it++) {
					int j = Qj[it];
					edges +=
					    spasm_row_weight(M,
							     i) *
					    spasm_row_weight(M, j);
				}
			assert(edges == spasm_nnz(A));

			if (mode == 't')
				printf("%d; %d; %d; %d; %d; %d; %d; %d\n", A->n,
				       spasm_nnz(A), trivial, nontrivial,
				       nontrivial_size, largest, module_edges,
				       quotient_edges);
			if (mode == 'v')
				printf("NotImplemented\n");
		};
		break;
	}

	spasm_csr_free(A);

	spasm_csr_free(partition->S);
	spasm_csr_free(partition->M);
	spasm_csr_free(partition->Q);
	free(partition);

	return 0;
}*/
