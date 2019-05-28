#include "mini_spasm.h"

struct modular_partition_t {
	spasm *S, *M, *Q;
};

struct node_t;
struct class_t;
struct tree_node_t;

struct node_t {
	struct node_t *prev, *next;
	int vertex;
	struct class_t *class;
};

struct class_t {
	struct class_t *prev, *next;
	struct node_t *nodes;
	int begin_index;
	int size;
	int marks;
	int split;
	int insert;

	int Lpos;
	int Kpos;
	int check;
};

struct module_ctx_t {
	spasm *A;
	struct node_t *nodes;
	struct class_t **L;
	struct class_t **K;
	int L_sp, K_lo, K_hi;

};

struct tree_node_t {
	int n_class;
	int n_vertex;
	struct class_t **classes;
	int *old_vertex;
	struct tree_node_t *father;
	struct tree_node_t *childs;
	spasm *matrix;
};

struct modular_partition_t *modular_partition(spasm * A);

struct tree_node_t *make_root(struct class_t *initial_class);
struct tree_node_t **make_childs(struct tree_node_t *root);
