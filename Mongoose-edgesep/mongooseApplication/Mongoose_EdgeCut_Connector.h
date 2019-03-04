#ifndef MONGOOSE_EDGECUT_CONNECTOR
#define MONGOOSE_EDGECUT_CONNECTOR

#ifndef __cplusplus
typedef unsigned char bool;
static const bool False = 0;
static const bool True = 1;
#endif

#ifdef __cplusplus
extern "C" {
#endif

enum MatchingStrategy
{
    Random   = 0,
    HEM      = 1,
    HEMSR    = 2,
    HEMSRdeg = 3
};

enum InitialEdgeCutType
{
    InitialEdgeCut_QP           = 0,
    InitialEdgeCut_Random       = 1,
    InitialEdgeCut_NaturalOrder = 2
};


typedef struct {
    bool *partition;
    int64_t n;  

    double cut_cost;
    int64_t  cut_size;
    double w0;         
    double w1;         
    double imbalance;  
} EdgeCutC;

typedef struct {
    int64_t  n;
    int64_t m; 
    int64_t  nz; 
    int64_t  *p;
    int64_t  *i;
    double *x;
    double *w;
} GraphC;

EdgeCutC* connector_edge_cut(const GraphC* g);

#endif
