#include "Mongoose_EdgeCut_Connector.h"
#include "mini_spasm.h"

#include "modules.h"
#include "complement_cc.h"

struct connected_component_partition_t {
	spasm *A;
	double *w,*x;
};

struct connected_component_partition_t *component_partition(spasm * A);
