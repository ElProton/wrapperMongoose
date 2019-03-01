/* ========================================================================== */
/* === Source/Mongoose_Waterdance.cpp ======================================= */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * Mongoose Graph Partitioning Library  Copyright (C) 2017-2018,
 * Scott P. Kolodziej, Nuri S. Yeralan, Timothy A. Davis, William W. Hager
 * Mongoose is licensed under Version 3 of the GNU General Public License.
 * Mongoose is also available under other licenses; contact authors for details.
 * -------------------------------------------------------------------------- */

#include "../Include/Mongoose_Waterdance.hpp"
#include "../Include/Mongoose_ImproveFM.hpp"
#include "../Include/Mongoose_ImproveQP.hpp"
#include "../Include/Mongoose_Internal.hpp"

namespace Mongoose
{

void waterdance(EdgeCutProblem *graph, const EdgeCut_Options *options)
{
    Int numDances = options->num_dances;
    for (Int i = 0; i < numDances; i++)
    {
        improveCutUsingFM(graph, options);
        improveCutUsingQP(graph, options);
    }
}

} // end namespace Mongoose
