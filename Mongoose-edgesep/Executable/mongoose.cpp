#include "Mongoose.hpp"
#include <iostream>
#include <iomanip>
#include <math.h>

using namespace Mongoose;
using namespace std;

int main(int argn, const char** argv){
EdgeCut_Options* options = EdgeCut_Options::create();
if(!options) return EXIT_FAILURE;
options->matching_strategy = HEMSRdeg;
options->initial_cut_type = InitialEdgeCut_QP;
Graph* graph = read_graph(argv[1]);
if(!graph){
options->~EdgeCut_Options();
return EXIT_FAILURE;
}
EdgeCut* result = edge_cut(graph,options);
cout << "partitionning end" << endl;
cout << "cut cost:" << setprecision(2) << fabs(100 * result->imbalance) << "%" << endl;
options->~EdgeCut_Options();
graph->~Graph();
result->~EdgeCut();
return EXIT_SUCCESS;
}
