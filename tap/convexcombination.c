/*
   This file contains routines for MSA and Frank-Wolfe (both are "convex combination" types of
   algorithms, and supporting helper functions.

   Both of these implementations use an efficient network loading procedure for calculating the
   target solutions x*.  Rather than finding the shortest path for each OD pair, and then
   loading the demand on these paths one at a time, it loads demand on *all* shortest paths with the
   same origin simultaneously, which is more efficient -- O(#links) instead of O(#links^2).  SPtree
   is a network struct which assists with this process.

*/

#include "convexcombination.h"

/////////////////////
// MAIN ALGORITHMS //
/////////////////////

/*
FrankWolfe takes two arguments: a network for equilibrium assignment, and a struct of algorithm parameters
*/
void FrankWolfe(network_type *network, FrankWolfeParameters_type *parameters) {
	double gap;
	double elapsedTime = 0;
	double lambda = 1;
	long i, iteration = 0;
	/* Allocate memory */
	declareVector(double, oldFlows, network->numArcs);
	declareVector(double, SPFlows, network->numArcs);
	declareVector(double, tempFlows, network->numArcs);

	/* Build an auxiliary network which will store shortest paths from each origin, to speed up x* calculation */
	network_type SPtree;
	SPtree.numNodes = network->numNodes;
	SPtree.numArcs = network->numNodes - 1;
	SPtree.nodes = newVector(network->numNodes, node_type);
	SPtree.arcs = newVector(network->numNodes - 1, arc_type);

	clock_t stopTime = clock();

	if (parameters->warmStart == FALSE) /* Must initialize link flows using shortest paths unless warm-starting */
		initializeFlows(network, &SPtree, tempFlows);

	do {
		iteration++;
		elapsedTime += ((double)(clock() - stopTime)) / CLOCKS_PER_SEC; /* Exclude gap calculations from run time */
		gap = calculateGap(network, parameters->gapFunction);
		displayMessage(LOW_NOTIFICATIONS, "Iteration %ld: gap %.15f, Beckmann %.13g, time %.3f s.\n", iteration, gap, BeckmannFunction(network), elapsedTime);
		if (gap < parameters->convergenceGap) break;
		stopTime = clock();

		/* Store x values */
		for (i = 0; i < network->numArcs; i++) {
			oldFlows[i] = network->arcs[i].flow;
		}

      /* Find x* */
		getSPFlows(network, SPFlows, &SPtree, tempFlows);

		/* Find step size */
      lambda = FrankWolfeLineSearch(network, parameters, oldFlows, SPFlows);

      /* Update x and t values */
		for (i = 0; i < network->numArcs; i++) {
			network->arcs[i].flow = lambda * SPFlows[i] + (1 - lambda) * oldFlows[i];
		}
		updateAllCosts(network, BPR);
	} while (elapsedTime < parameters->maxTime && iteration < parameters->maxIterations);

   /* Clean up memory */
	deleteVector(oldFlows);
	deleteVector(tempFlows);
	deleteVector(SPFlows);
	deleteVector(SPtree.nodes);
	deleteVector(SPtree.arcs);
}


/*
initializeFrankWolfeParameters -- creates a default set of parameters for the Frank-Wolfe algorithm,
which can be adjusted later.  See main.c for documentation of these parameters.  In particular, be sure
to set at least one of convergenceGap, maxTime, or maxIterations, or the code will run forever.
*/
FrankWolfeParameters_type initializeFrankWolfeParameters() {
   FrankWolfeParameters_type parameters;

   parameters.gapFunction = RELATIVE_GAP_1;
   parameters.convergenceGap = 0;
   parameters.maxTime = INFINITY;
   parameters.maxIterations = LONG_MAX;
   parameters.warmStart = FALSE;
   parameters.lineSearchIterations = 11;
   parameters.lineSearch = BISECTION;

   return parameters;
}

/*
FrankWolfeLineSearch -- identifies the lambda value according to the parameters struct (either BISECTION or NEWTON's method).
Other arguments are the underlying network struct, the existing link flows x (oldFlows), and the target links flows x* (SPFlows)
*/
double FrankWolfeLineSearch(network_type *network, FrankWolfeParameters_type *parameters, double *oldFlows, double *SPFlows) {
   switch (parameters->lineSearch) {
      case BISECTION:   return FrankWolfeBisection(network, parameters, oldFlows, SPFlows);
      case NEWTON:      return FrankWolfeNewton(network, parameters, oldFlows, SPFlows);
   }
   fatalError("Unknown line search type %d\n", parameters->lineSearch);
   return IS_MISSING; /* Should never be reached; avoids compiler warning */
}

/*
FrankWolfeBisection -- Implements the bisection method for a fixed number of iterations, as specified in the parameters struct.
Other arguments are the underlying network struct, the existing link flows x (oldFlows), and the target links flows x* (SPFlows)
*/
double FrankWolfeBisection(network_type *network, FrankWolfeParameters_type *parameters, double *oldFlows, double *SPFlows) {
   int i, iteration;
   double lambda, lmax, lmin, der;

   lmin = 0; lmax = 1;
   for (iteration = 0; iteration < parameters->lineSearchIterations; iteration++) {
      lambda = (lmax + lmin) / 2;
      der = 0;
      for (i = 0; i < network->numArcs; i++) {
         network->arcs[i].flow = (1 - lambda) * oldFlows[i] + lambda * SPFlows[i];
         der += linkCost(&network->arcs[i], BPR) * (SPFlows[i] - oldFlows[i]);
      }
      if (der > 0) lmax = lambda; else lmin = lambda;
   }

   return lambda;
}

/*
FrankWolfeBisection -- Implements Newton's method for a fixed number of iterations, as specified in the parameters struct.
Other arguments are the underlying network struct, the existing link flows x (oldFlows), and the target links flows x* (SPFlows)
*/
double FrankWolfeNewton(network_type *network, FrankWolfeParameters_type *parameters, double *oldFlows, double *SPFlows) {
   int i, iteration;
   double num, denom, lambda;

   for (iteration = 0; iteration < parameters->lineSearchIterations; iteration++) {
      num = 0;
      denom = 0;
      for (i = 0; i < network->numArcs; i++) {
         num += linkCost(&network->arcs[i], BPR) * (SPFlows[i] - oldFlows[i]);
         denom += linkCost(&network->arcs[i], BPR_DER) * (SPFlows[i] - oldFlows[i]) * (SPFlows[i] - oldFlows[i]);
      }
      if (denom != 0) { /* Typical case */
         lambda = - num / denom;
         lambda = min(lambda, 1);
         lambda = max(lambda, 0);
      } else { /* Denominator is zero; indicates constant link performance functions so all flow should be shifted to shorter path */
         lambda = 1;
      }
      for (i = 0; i < network->numArcs; i++) { /* Update based on new value of lambda for next iteration of Newton's method */
         network->arcs[i].flow = (1 - lambda) * oldFlows[i] + lambda * SPFlows[i];
      }
   }

   return lambda;
}

/*
MSA takes two arguments: a network for equilibrium assignment, and a struct of algorithm parameters
*/
void MSA(network_type *network, MSAparameters_type *parameters) {

	double gap;
	double elapsedTime = 0;
	long i, iteration = 0;

	/* Allocate memory */
	declareVector(double, SPFlows, network->numArcs);
	declareVector(double, tempFlows, network->numArcs);
	network_type SPtree;

	/* Build an auxiliary network which will store shortest paths from each origin, to speed up x* calculation */
	SPtree.numNodes = network->numNodes;
	SPtree.numArcs = network->numNodes - 1;
	SPtree.nodes = newVector(network->numNodes, node_type);
	SPtree.arcs = newVector(network->numNodes - 1, arc_type);

	clock_t stopTime = clock();
	if (parameters->warmStart == FALSE) /* Must initialize link flows using shortest paths unless warm-starting */
		initializeFlows(network, &SPtree, tempFlows);

	do {
		iteration++;
		elapsedTime += ((double)(clock() - stopTime)) / CLOCKS_PER_SEC; /* Exclude gap calculations from run time */
		gap = calculateGap(network, parameters->gapFunction);
		displayMessage(LOW_NOTIFICATIONS, "Iteration %ld: gap %.15f, Beckmann %.13g, time %.3f s.\n", iteration, gap, BeckmannFunction(network), elapsedTime);
		if (gap < parameters->convergenceGap) break;
		stopTime = clock();

      /* Find x* */
		getSPFlows(network, SPFlows, &SPtree, tempFlows);

		/* Update flow values using MSA rule, lambda = 1/iteration */
		for (i = 0; i < network->numArcs; i++) {
			network->arcs[i].flow = (1 - 1.0 / iteration) * network->arcs[i].flow + (1.0 / iteration) * SPFlows[i]; }
		updateAllCosts(network, BPR);
	} while (elapsedTime < parameters->maxTime && iteration < parameters->maxIterations);

   /* Clean up memory */
	deleteVector(tempFlows);
	deleteVector(SPFlows);
	deleteVector(SPtree.nodes);
	deleteVector(SPtree.arcs);
}

/*
initializeMSAParameters -- creates a default set of parameters for the Frank-Wolfe algorithm,
which can be adjusted later.  See main.c for documentation of these parameters.  In particular, be sure
to set at least one of convergenceGap, maxTime, or maxIterations, or the code will run forever.
*/
MSAparameters_type initializeMSAparameters() {
   MSAparameters_type parameters;

   parameters.gapFunction = RELATIVE_GAP_1;
   parameters.convergenceGap = 0;
   parameters.maxTime = INFINITY;
   parameters.maxIterations = LONG_MAX;
   parameters.warmStart = FALSE;

   return parameters;
}

//////////////////////
// HELPER FUNCTIONS //
//////////////////////

/*
getSPFlows calculates the target flow vector x* by iterating over each origin; loadOrigin finds the
contribution to x* from each origin (tempFlows), and these flows are summed together to get the total
x* vector.
Arguments: the underlying network; an array SPFlows which will store the new x* values;
   and the SPtree network and tempFlows array.  The latter two are passed by reference from
   the high-level MSA/FW function to avoid overhead of repeatedly allocating and deallocating dynamic
   memory.
*/
void getSPFlows(network_type *network, double *SPFlows, network_type *SPtree, double *tempFlows) {
	long i, j;
	for (i = 0; i < network->numArcs; i++) {
		SPFlows[i] = 0;
	}
	for (i = 0; i < network->numZones; i++) {
		loadOrigin(i, network, SPtree, tempFlows);
		for (j = 0; j < network->numArcs; j++) {
			SPFlows[j] += tempFlows[j];
		}
	}
}


/*
initializeFlows finds an initial solution for MSA or Frank-Wolfe, by assigning all flow to
shortest paths.  Implementation is quite similar to getSPFlows, except that we are working
directly with the network flows x (rather than a temporary array of SPFlows for x*), and
recalculating travel times after loading each origin.
Arguments: the underlying network, and the SPtree network and tempFlows array.  The latter two are passed
   by reference from the high-level MSA/FW function to avoid overhead of repeatedly allocating and
   deallocating dynamic memory.
*/
void initializeFlows(network_type *network, network_type *SPtree, double *newArcFlows) {
	long i, j;
	for (i = 0; i < network->numArcs; i++) {
		network->arcs[i].flow = 0;
		network->arcs[i].cost = network->arcs[i].freeFlowTime;
	}

	for (i = 0; i < network->numZones; i++) {
		updateAllCosts(network, BPR);
		loadOrigin(i, network, SPtree, newArcFlows);
		for (j = 0; j < network->numArcs; j++) {
			network->arcs[j].flow += newArcFlows[j];
		}
	}

}

/*
loadOrigin Loads all demand from an origin onto its shortest path tree.  Using the fact that the shortest path tree is
acyclic, we can do this efficiently by proceeding in reverse topological order, identifying the total number of vehicles
from this origin using each link.  Because of Bellman's principle, there is significant overlap among the shortest paths
which start from the same origin.  This implementation avoids the redundant effort of loading flow among overlapping paths.

SPtree is a network_type consisting of exactly numNodes nodes and numNodes - 1 arcs... run topologicalOrder on it to assign demand.
It is passed by reference to loadOrigin so it can be allocated/deallocated properly in whatever higher-level function calls loadOrigin.
Arguments: the origin to load, the network struct, and the SPtree network and newArcFlows arrays.
*/
void loadOrigin(long origin, network_type *network, network_type *SPtree, double *newArcFlows) {
	long curnode, i, j;

	declareVector(double, label, network->numNodes);
	declareVector(arc_type *, backarc, network->numNodes);

	/* Find all-to-one shortest paths from origin */
	arcBellmanFord(origin, label, backarc, network, DEQUE);

	/* Update network struct with OD shortest path travel times */
	for (i = 0; i < network->numZones; i++) {
		network->OD[origin][i].cost = label[i];
	}

	/* Update the SPtree network based on shortest paths found from this origin.
	   As a tree we immediately know each node (except the origin) has exactly one
	   incoming link; as a shortest path tree we know this link must be the
	   predecessor given by the shortest path algorithm */
	j = 0;
	for (i = 0; i < network->numNodes; i++) {
		if (i == origin) continue;
		SPtree->arcs[j].tail = backarc[i]->tail;
		SPtree->arcs[j].head = i;
		j++;
	}

	/* Redo the forward/reverse star lists for the network */
	finalizeNetwork(SPtree);

	/* Find topological order of the shortest path tree */
	declareVector(long, invNodeOrder, network->numNodes);
	declareVector(long, nodeOrder, network->numNodes);
	topologicalOrder(SPtree, invNodeOrder, nodeOrder);

   /* Now load vehicles onto this tree in reverse topological order -- only one sweep of the tree
      is needed to find all flow from this origin.  remainingVehicles gives the number of vehicles
      at each node which have not yet been fully assigned to their path.
      In the TNTP file format, origins are numbered first.  The second loop thus picks up where
      the first one left off.
      */
	declareVector(double, remainingVehicles, network->numNodes);
	for (i = 0; i < network->numZones; i++) remainingVehicles[i] = network->OD[origin][i].demand;
	for (; i < network->numNodes; i++) remainingVehicles[i] = 0;
	for (i = 0; i < network->numArcs; i++) newArcFlows[i] = 0;

	/* Here is the main loop, in reverse topological order */
	for (i = network->numNodes - 1; i > 0; i--) {
		curnode = invNodeOrder[i];
		if (curnode == origin) break;
		if (backarc[curnode] == NULL && remainingVehicles[curnode] > 0) { /* Uh oh! -- this shouldn't happen if the network is built properly. */
			if (verbosity >= FULL_DEBUG) { /* Dump extra info to help debug if the verbosity level is high enough */
				displayMessage(DEBUG, "Node\tbackarc\tremaining vehicles\tSP cost\n");
				for (j = 0; j < network->numNodes; j++) {
					curnode = nodeOrder[j];
					if (backarc[curnode] != NULL)
						displayMessage(DEBUG, "%ld\t(%ld,%ld)\t%f\t%f\n", curnode, backarc[curnode]->tail, backarc[curnode]->head, remainingVehicles[curnode], label[curnode]);
					else
						displayMessage(DEBUG, "%ld\t(---,---)\t%f\t%f\n", curnode, remainingVehicles[curnode], label[curnode]);
				}
			}
			fatalError("LoadOrigin -- no path from %d to %d even though positive demand %f exists there!", origin, curnode, remainingVehicles[curnode]);
		}
		newArcFlows[backarc[curnode] - network->arcs] = remainingVehicles[curnode];
		remainingVehicles[backarc[curnode]->tail] += remainingVehicles[curnode];
		remainingVehicles[curnode] = 0;
	}

   /* Clean up memory */
	deleteVector(label);
	deleteVector(backarc);
	deleteVector(invNodeOrder);
	deleteVector(nodeOrder);
	deleteVector(remainingVehicles);
	for (i = 0; i < network->numNodes; i++) {
		clearArcList(&(SPtree->nodes[i].forwardStar));
		clearArcList(&(SPtree->nodes[i].reverseStar));
	}
}
