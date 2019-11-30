/*
   This file contains routines applicable across many types of algorithms for traffic assignment:
   calculation of Beckmann function, gap functions (relative gap, average excess cost), evaluation
   of BPR functions, total system travel time, and so forth.
*/

#include "tap.h"

//////////////////////////////////
// General network calculations //
//////////////////////////////////

/*
arcNumber converts a pointer to an arc to the index number for that arc; to do this, the network struct
needs to be passed along with the arc pointer.
*/
inline long arcNumber(network_type *network, arc_type *arc) {
	return arc - network->arcs;
}

/*
BeckmannFunction calculates the Beckmann function for a network given its current flows.
*/
double BeckmannFunction(network_type *network) {
	double Beckmann = 0;
	long i;
	for (i = 0; i < network->numArcs; i++) {
		Beckmann += linkCost(&(network->arcs[i]), BPR_INT); }
	return Beckmann;
}

/*
calculateGap is a wrapper function for the various gap calculation functions.  Argument gapFunction
indicates which function to call.  RELATIVE_GAP_1 is what I taught you in class.  RELATIVE_GAP_2 is
an alternative definition of relative gap which appears sometimes in the literature.  AEC is average
excess cost.  AEC_OB and MEC are other gap functions which I have not yet implemented.
*/
double calculateGap(network_type *network, gap_type gapFunction) {
	updateAllCosts(network, BPR);
	switch (gapFunction) {
		case RELATIVE_GAP_1: return relativeGap1(network); break;
		case RELATIVE_GAP_2: return relativeGap2(network); break;
		case AEC: return averageExcessCost(network); break;
		case AEC_OB: fatalError("Gap type AEC_OB not yet implemented."); break;
		case MEC: fatalError("Gap type MEC not yet implemented."); break;
		default: fatalError("Unknown gap type %d\n", gapFunction);
	}
	return IS_MISSING; /* This should never be reached; provided to eliminate compiler warnings */
}

/*
linkCost performs a link-specific calculation, specified by the argument costFunction.  Most commonly used
to evaluate a link performance function such as the BPR function, but passing different costFunctions allows
you to calculate other quantities with the same function.  The following costFunctions have been implemented
thus far:
   BPR -- evaluates a BPR function using the arc'link s alpha and beta values, its capacity, and its current flow
   BPR_DER -- evaluates the *derivative* of the BPR function at the current flow (useful for Newton's method)
   BPR_INT -- evaluates the *integral* of the BPR function, between zero and the current flow (useful for Beckmann function!)
   SO_BPR -- evaluates the marginal cost of the link (flow + externality).  By using this instead of the BPR function, the
             system optimal assignment can be solved with the existing code and no other changes.
*/
double linkCost(arc_type *arc, cost_type costFunction) {
   if (arc->capacity <= 0) return INFINITY; /* Protect against division by zero */
	switch (costFunction) {
	case BPR:
      if (arc->flow <= 0) return arc->freeFlowTime + arc->fixedCost; // Protect against negative flow values and 0^0 errors
		return arc->fixedCost + arc->freeFlowTime * (1 + arc->alpha * pow(arc->flow / arc->capacity, arc->beta));
	case BPR_DER:
      if (arc->flow <= 0) { /* Protect against negative flow values and 0^0 errors */
         if (arc->beta != 0)
            return 0;
         else
            return arc->freeFlowTime * arc->alpha / arc->capacity;
      }
		return arc->freeFlowTime * arc->alpha * arc->beta / arc->capacity
			* pow(arc->flow / arc->capacity, arc->beta - 1);
	case BPR_INT:
      if (arc->flow <= 0) return 0; /* Protect against negative flow values and 0^0 errors */
      return arc->flow * (arc->fixedCost + arc->freeFlowTime * (1 + arc->alpha / (arc->beta + 1) * pow(arc->flow / arc->capacity, arc->beta)));
	case SO_BPR:
		return arc->freeFlowTime * (1 + arc->alpha * (1 + arc->beta) * pow(arc->flow / arc->capacity, arc->beta));
	default:
		fatalError("Unknown cost function.");
		return IS_MISSING;
	}
}


/*
SPTT calculates the shortest-path travel time on the network, that is, the total system travel time if everyone could be
loaded on the current shortest paths without changing the travel times.  This function actually re-solves shortest paths
for each origin, rather than using the OD cost field in the network struct -- this allows SPTT to be calculated without
changing any values in the network, at the expense of a little more run time.
*/
double SPTT(network_type *network) {
	long i, j;
	double sptt = 0;
	declareVector(double, SPcosts, network->numNodes);
	declareVector(long, backnode, network->numNodes);
	for (i = 0; i < network->numZones; i++) {
		BellmanFord(i, SPcosts, backnode, network, DEQUE);
		for (j = 0; j < network->numZones; j++) {
			sptt += network->OD[i][j].demand * SPcosts[j];
		}
	}
	deleteVector(SPcosts);
	deleteVector(backnode);
	return sptt;
}

/*
TSTT calculates the total system travel time on the network, that is, the dot product of link flow and travel time.
*/
double TSTT(network_type *network) {
	double sum = 0;
	long i;
	for (i = 0; i < network->numArcs; i++)
		sum += network->arcs[i].flow * linkCost(&network->arcs[i], BPR);
   if (isnan(sum)) displayNetwork(DEBUG, network); /* Oops.  Indicates some kind of numerical error (likely division by zero or 0^0) */
	return sum;
}

/*
updateAllCosts recalculates and updates all link costs in the network, according to costFunction
*/
void updateAllCosts(network_type *network, cost_type costFunction) {
	long i;
	for (i = 0; i < network->numArcs; i++)
		network->arcs[i].cost = linkCost(&network->arcs[i], costFunction);
}

/*
updateAllCostDers recalculates and updates all link derivatives in the network, according to derFunction
*/
void updateAllCostDers(network_type *network, cost_type derFunction) {
	long i;
	for (i = 0; i < network->numArcs; i++)
		network->arcs[i].der = linkCost(&network->arcs[i], derFunction);
}

///////////////////
// GAP FUNCTIONS //
///////////////////

/*
averageExcessCost calculates the difference between TSTT and SPTT, normalized by total demand in the network.
*/
double averageExcessCost(network_type *network) {
	double sptt = SPTT(network), tstt = TSTT(network);
	if (tstt < sptt) warning(LOW_NOTIFICATIONS, "Negative gap.  TSTT and SPTT are %f %f\n", tstt, sptt);
	return ((tstt - sptt) / network->totalODFlow);
}

/*
relativeGap1 calculates the ratio between (TSTT - SPTT) and SPTT.
*/
double relativeGap1(network_type *network) {
	double sptt = SPTT(network), tstt = TSTT(network);
	displayMessage(DEBUG, "Current relative gap:\nCurrent TSTT: %f\nShortest path TSTT: %f\n", tstt, sptt);
	if (tstt < sptt) warning(LOW_NOTIFICATIONS, "Negative gap.  TSTT and denom are %f %f\n", tstt, sptt);
	return (tstt / sptt - 1);
}

/*
relativeGap2 calculates the ratio between the current value of the Beckmann function and a lower bound on the
optimal value of the Beckmann function, which is calculated using convexity: the tangent plane at any point
always underestimates the true value of the Beckmann function, and one such value is given by the current value
of the Beckmann function, minus (TSTT - SPTT).  (This is because the gradient of the Beckmann function is just
the current link travel times; evaluating the tangent plane approximation at the target link flow solution x*
gives this formula.)  Since any lower bound will do, this function stores the best lower bound seen thus far.
*/
double relativeGap2(network_type *network) {
	double sptt = SPTT(network), tstt = TSTT(network);

	network->beckmann = BeckmannFunction(network);
	network->beckmannLB = min(network->beckmannLB, network->beckmann + sptt - tstt);
	displayMessage(DEBUG, "Current relative gap:\nCurrent TSTT: %f\nShortest path TSTT: %f\n", tstt, sptt);
	return (network->beckmann / network->beckmannLB - 1);
}

