#ifndef BUSH_H
#define BUSH_H

#include <limits.h>
#include <math.h>
#include "tap.h"
#include "networks.h"
#include "datastructures.h"
#include "utils.h"

typedef struct {
   gap_type gapFunction;
   double   convergenceGap;
   double   maxTime;
   long     maxIterations;
   int      innerIterations;
   double   minLinkFlow;
   double   minDerivative;
   double   newtonStep;
   double   nodeGapTolerance;
} algorithmBParameters_type;

/* Algorithm B functions */
algorithmBParameters_type initializeAlgorithmBParameters();
void AlgorithmB(network_type *network, algorithmBParameters_type *parameters);
void updateBushB(long origin, network_type *network, bool *inBush, double *bushFlows, double *SPcosts, double *LPcosts, long *bushOrder, long *SPtree, long *LPtree, algorithmBParameters_type *parameters);
void updateFlowsB(long origin, network_type *network, bool* inBush, long* bushOrder, double *bushFlows, algorithmBParameters_type *parameters);
void initializeOriginB(long origin, network_type *network, bool *inBush, long *bushOrder, double*bushFlows, algorithmBParameters_type *parameters);

/* General bush functions */
void findBushSP(network_type *network, bool *inBush, long *bushOrder, long *SPtree, double *SPcosts);
void findBushLP(network_type *network, bool *inBush, long *bushOrder, long *SPtree, double *SPcosts);
void findBushUsedLP(network_type *network, bool *inBush, long *bushOrder, double *bushFlows, long *LPtree, double *LPcosts, algorithmBParameters_type *parameters);
void recalculateBushFlows(network_type *network, bool *inBush, double *bushFlows, long *bushOrder, long *SPtree, algorithmBParameters_type *parameters);
void bushTopologicalOrder(network_type *network, bool *inBush, long *bushOrder);

#endif
