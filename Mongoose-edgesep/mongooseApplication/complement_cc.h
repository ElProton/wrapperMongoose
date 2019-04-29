#include "mini_spasm.h"

struct connected_component_t {
    int *component_list;
    int n;
};

int get_components_from_complement(spasm *A);
