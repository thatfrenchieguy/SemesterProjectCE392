#ifndef TAP_H
#define TAP_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "fileio.h"
#include "networks.h"
#include "utils.h"
#include "datastructures.h"

#define IS_MISSING	-1

typedef enum {
	BPR,
	BPR_DER,
	BPR_INT,
	SO_BPR
} cost_type;

typedef enum {
	RELATIVE_GAP_1,
	RELATIVE_GAP_2,
	AEC,
	AEC_OB,
	MEC
} gap_type;


//////////////////////////////////
// General network calculations //
//////////////////////////////////

double BeckmannFunction(network_type *network);
double calculateGap(network_type *network, gap_type gapFunction);
void finalizeNetwork(network_type *network);
double linkCost(arc_type *arc, cost_type costFunction);
double SPTT(network_type *network);
double TSTT(network_type *network);
void updateAllCosts(network_type *network, cost_type costFunction);
void updateAllCostDers(network_type *network, cost_type costFunction);

inline long arcNumber(network_type *network, arc_type *arc);
double averageExcessCost(network_type *network);
double relativeGap1(network_type *network);
double relativeGap2(network_type *network);

#endif
