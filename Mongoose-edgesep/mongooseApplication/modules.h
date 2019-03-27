#include "mini_spasm.h"

struct modular_partition_t {
	spasm *S, *M, *Q;
};

struct node_t;
struct class_t;

struct node_t {
	struct node_t *prev, *next;
	int vertex;
	struct class_t *class;
	int index;
};

struct class_t {
	struct class_t *prev, *next;
	struct node_t *nodes;
	int begin_index;
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

struct modular_partition_t *modular_partition(spasm * A);
