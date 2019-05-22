#include "Mongoose_EdgeCut_Connector.h"
#include "mini_spasm.h"

#include "modules.h"

struct component_partition_t {
	spasm *A;
	double *w,*x;
};

struct component_partition_t *component_partition(spasm * A);
