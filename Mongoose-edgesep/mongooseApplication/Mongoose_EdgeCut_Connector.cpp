#include <cstdlib>
#include <cmath>

#include "Mongoose_EdgeCut.hpp"
#include "Mongoose_Graph.hpp"
#include "Mongoose_EdgeCut_Connector.h"
#include "SuiteSparse_config.h"



#ifdef __cplusplus
extern "C" {
#endif

EdgeCutC* connector_edge_cut(const GraphC *g){
    Mongoose::EdgeCut* ec;

    EdgeCutC *res = (EdgeCutC*)SuiteSparse_malloc(1, sizeof(*res));

    Mongoose::Graph *g_obj = Mongoose::Graph::create(fmax(g->n,g->m),g->nz,g->p,g->i,g->x,g->w);
    Mongoose::EdgeCut_Options *options = Mongoose::EdgeCut_Options::create();

    options->matching_strategy = (Mongoose::MatchingStrategy) HEMSRdeg;
    options->initial_cut_type = (Mongoose::InitialEdgeCutType) InitialEdgeCut_QP;


    ec = edge_cut(g_obj,options);

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
