#include <cstdlib>

#include "../Include/Mongoose_EdgeCut.hpp"
#include "../Include/Mongoose_Graph.hpp"
#include "../Include/Mongoose_EdgeCut_Connector.h"
#include "../SuiteSparse_config/SuiteSparse_config.h"
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

EdgeCutC* connector_edge_cut(const GraphC *g){
    Mongoose::EdgeCut* ec;

    EdgeCutC* res = (EdgeCutC*)SuiteSparse_malloc(1, sizeof(EdgeCutC));
    Mongoose::Graph *g_obj = Mongoose::Graph::create(fmax(g->n,g->m),g->nz,g->p,g->i,g->x,g->w);
    


    ec = edge_cut(g_obj);

    res->partition =  ec->partition;

    res->n = (int64_t) ec->n;
    res->cut_cost = ec->cut_cost;
    res->cut_size = (int64_t) ec->cut_size;
    res->w0 = ec->w0;
    res->w1 = ec->w1;
    res->imbalance = ec->imbalance;

    return res;
}



#ifdef __cplusplus
}}
#endif
