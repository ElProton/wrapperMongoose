#include "mini_spasm.h"

struct modular_partition_t {
	spasm *S, *M, *Q;
};

struct node_t;
struct class_t;
struct tree_node_t;
struct permutation_factorisante_t;
struct permutation_item_t;

struct node_t {
	struct node_t *prev, *next;
	int vertex;
	int real_vertex_value;
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

struct permutation_factorisante_t {
	struct permutation_item_t **items;
	int n;
	int index;
};

struct permutation_item_t {
	int vertex;
	int left_parentheses;
	int right_parentheses;
};

void print_permutation(struct permutation_factorisante_t *permutation);
struct modular_partition_t *modular_partition(spasm * A);
struct permutation_factorisante_t *make_permutation(struct class_t *class_x, struct node_t *nodes, spasm *A);
void *make_permutation_rec(struct class_t *class_x, struct permutation_factorisante_t *permutation, struct node_t *nodes, spasm *A);
struct permutation_factorisante_t *search_fractures(struct permutation_factorisante_t *permutation, spasm *A);

struct tree_node_t *make_root(struct class_t *initial_class);
struct tree_node_t **make_childs(struct tree_node_t *root);
