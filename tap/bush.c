/*
   This file contains routines for Algorithm B and supporting helper functions.  Some of these may be
   useful in implementing other bush-based algorithms such as OBA or LUCE.

*/
#include "bush.h"

/*
AlgorithmB: this function steers the algorithm at a high level; the real work is done
by the updateBushB and updateFlowsB functions.

bushFlows, inBush, and bushOrder store the bushes themselves, as 2D arrays. inBush[r][ij] is TRUE if
link (i,j) is in the bush for origin r, and FALSE otherwise; bushFlows[r][ij] shows the flow from
origin r on link (i,j); and bushOrder[r][i] gives the number of the i-th node in topological order for
origin r.  (This is *not* the order of the i-th node; bushOrder is implemented this way so that
iterating over bushOrder[r][0], bushOrder[r][1] ... bushOrder[r][network->numNodes - 1] will cover
all nodes in topological order.  This is more useful for B than looking up the index for a particular node.)

A few arrays (LPtree, SPcosts, etc.) are defined here; these will hold the backnode and cost labels
for the longest (LP) and shortest (SP) path computations.  These calculations are run separately for
each origin, so the arrays are only defined once and passed by reference to updateFlowsB and updateBushB,
saving memory and time.
*/
void AlgorithmB(network_type *network, algorithmBParameters_type *parameters) {

	/* Strong connectivity check */
	makeStronglyConnectedNetwork(network);

   /* Allocate memory for bushes */
	long origin, i, iteration = 0;
	declareMatrix(double, bushFlows, network->numZones, network->numArcs);
	declareMatrix(bool, inBush, network->numZones, network->numArcs);
	declareMatrix(long, bushOrder, network->numZones, network->numNodes);
	declareVector(long, LPtree, network->numNodes);
	declareVector(long, SPtree, network->numNodes);
	declareVector(double, LPcosts, network->numNodes);
	declareVector(double, SPcosts, network->numNodes);

	double elapsedTime = 0, gap = INFINITY;

   /* Initialize */
	clock_t stopTime = clock(); /* used for timing */
	for (i = 0; i < network->numArcs; i++) network->arcs[i].flow = 0;
	for (origin = 0; origin < network->numZones; origin++) {
		updateAllCosts(network, BPR);
		initializeOriginB(origin, network, inBush[origin], bushOrder[origin], bushFlows[origin], parameters);
	}
	for (i = 0; i < network->numArcs; i++) {
		network->arcs[i].cost = linkCost(&network->arcs[i], BPR);
		network->arcs[i].der = linkCost(&network->arcs[i], BPR_DER);
	}

   /* Iterate */
	do {
	   iteration++;

      /* Update bushes */
		for (origin = 0; origin < network->numZones; origin++) {
			updateBushB(origin, network, inBush[origin], bushFlows[origin], SPcosts, LPcosts, bushOrder[origin], SPtree, LPtree, parameters);
			updateFlowsB(origin, network, inBush[origin], bushOrder[origin], bushFlows[origin], parameters);
		}

      /* Shift flows */
		for (i = 0; i < parameters->innerIterations; i++) {
			for (origin = 0; origin < network->numZones; origin++) {
				updateFlowsB(origin, network, inBush[origin], bushOrder[origin], bushFlows[origin], parameters);
			}
		}

      /* Check gap and report progress */
		elapsedTime += ((double)(clock() - stopTime)) / CLOCKS_PER_SEC; /* Exclude gap calculations from run time */
		gap = calculateGap(network, parameters->gapFunction);
		displayMessage(LOW_NOTIFICATIONS, "Iteration %ld: gap %.15f, Beckmann %.13g, time %.3f s.\n", iteration, gap, BeckmannFunction(network), elapsedTime);
		stopTime = clock();
	} while (elapsedTime < parameters->maxTime && iteration < parameters->maxIterations && gap > parameters->convergenceGap);

   /* Clean up */
	deleteMatrix(bushFlows, network->numZones);
	deleteMatrix(inBush, network->numZones);
	deleteMatrix(bushOrder, network->numZones);
	deleteVector(LPtree);
	deleteVector(SPtree);
	deleteVector(LPcosts);
	deleteVector(SPcosts);

}

/*
initializeAlgorithmBParameters sets up the parameters struct for B and fills in default values.
*/
algorithmBParameters_type initializeAlgorithmBParameters() {
   algorithmBParameters_type parameters;

   parameters.gapFunction = RELATIVE_GAP_1;
   parameters.convergenceGap = 0;
   parameters.maxTime = INFINITY;
   parameters.maxIterations = LONG_MAX;

   parameters.innerIterations = 20;
   parameters.minLinkFlow = 1e-14;
   parameters.minDerivative = 1e-6;
   parameters.newtonStep = 1;
   parameters.nodeGapTolerance = 1e-14;

   return parameters;
}


/*
findBushSP implements shortest-path finding on bushes, taking advantage of the topological order.
Arguments are a pointer to the underlying network, the inBush and bushOrder arrays for the relevant
origin, and the SPtree and SPcosts arrays which will return the backnode and cost labels.
*/
void findBushSP(network_type *network, bool *inBush, long *bushOrder, long *SPtree, double *SPcosts) {
   arcListElt *curArc;
	long i, ij, curnode;

	for (i = 0; i < network->numNodes; i++) {
		SPcosts[i] = INFINITY;
		SPtree[i] = NO_PATH_EXISTS;
	}

	SPcosts[bushOrder[0]] = 0;
	for (i = 0; i < network->numNodes; i++) {
		curnode = bushOrder[i];
		for (curArc = network->nodes[curnode].forwardStar.head; curArc != NULL; curArc = curArc->next) {
         ij = ptr2arc(network, curArc->arc);
			if (inBush[ij] == TRUE && SPcosts[curnode] + network->arcs[ij].cost < SPcosts[network->arcs[ij].head]) {
				SPcosts[network->arcs[ij].head] = SPcosts[curnode] + network->arcs[ij].cost;
				SPtree[network->arcs[ij].head] = ij;
			}
		}
	}
}

/*
findBushLP implements longest-path finding on bushes, taking advantage of the topological order.
Arguments are a pointer to the underlying network, the inBush and bushOrder arrays for the relevant
origin, and the SPtree and SPcosts arrays which will return the backnode and cost labels.
*/
void findBushLP(network_type *network, bool *inBush, long *bushOrder, long *LPtree, double *LPcosts) {
   arcListElt *curArc;
	long i, ij, curnode;

	for (i = 0; i < network->numNodes; i++) {
		LPcosts[i] = -INFINITY;
		LPtree[i] = NO_PATH_EXISTS;
	}

	LPcosts[bushOrder[0]] = 0;
	for (i = 0; i < network->numNodes; i++) {
		curnode = bushOrder[i];
		for (curArc = network->nodes[curnode].forwardStar.head; curArc != NULL; curArc = curArc->next) {
         ij = ptr2arc(network, curArc->arc);
         if (inBush[ij] == TRUE && LPcosts[curnode] + network->arcs[ij].cost > LPcosts[network->arcs[ij].head]) {
				LPcosts[network->arcs[ij].head] = LPcosts[curnode] + network->arcs[ij].cost;
				LPtree[network->arcs[ij].head] = ij;
			}
		}
	}

	for (i = 0; i < network->numNodes; i++) {
	   if (LPtree[i] == NO_PATH_EXISTS && i != bushOrder[0]) LPcosts[i] = INFINITY;
	}
}

/*
findBushUsedLP finds the longest *used* paths on bushes (those with positive flow), taking advantage
of the topological order.  Arguments are a pointer to the underlying network, the inBush, bushOrder,
and bushFlows  arrays for the relevant origin, and the SPtree and SPcosts arrays which will return
the backnode and cost labels.  The parameters struct is also passed in case you wish to only count
"used" links as those with flows greater than a small epsilon, to guard against numerical errors.
*/
void findBushUsedLP(network_type *network, bool *inBush, long *bushOrder, double *bushFlows, long *LPtree, double *LPcosts, algorithmBParameters_type *parameters) {
   arcListElt *curArc;
	long i, ij, curnode;

	for (i = 0; i < network->numNodes; i++) {
		LPcosts[i] = -INFINITY;
		LPtree[i] = NO_PATH_EXISTS;
	}

	LPcosts[bushOrder[0]] = 0;
	for (i = 0; i < network->numNodes; i++) {
		curnode = bushOrder[i];
		if (LPcosts[curnode] == -INFINITY) continue;
		for (curArc = network->nodes[curnode].forwardStar.head; curArc != NULL; curArc = curArc->next) {
         ij = ptr2arc(network, curArc->arc);
         if (LPcosts[curnode] + network->arcs[ij].cost > LPcosts[network->arcs[ij].head] && bushFlows[ij] > parameters->minLinkFlow) {
				LPcosts[network->arcs[ij].head] = LPcosts[curnode] + network->arcs[ij].cost;
				LPtree[network->arcs[ij].head] = ij;
			}
		}
	}
}

/*
updateBushB updates bushes by removing unused links, finding shortest/longest paths, and adding shortcuts
*/
void updateBushB(long origin, network_type *network, bool *inBush, double *bushFlows, double *SPcosts, double *LPcosts, long *bushOrder, long *SPtree, long *LPtree, algorithmBParameters_type *parameters) {
	long i, ij, j, newarcs = 0;

	/* Find shortest bush path and check reachability */
	findBushSP(network, inBush, bushOrder, SPtree, SPcosts);
	for (i = 0; i < network->numNodes; i++) {
	   if (SPcosts[i] == INFINITY && i != origin) {
	      warning(FULL_NOTIFICATIONS, "%lu unreachable from %lu on bush shortest path\n", i, origin);
      }
	}

	/* Remove arcs with near-zero flow, but keep the SP tree arc to retain connectivity */
	for (ij = 0; ij < network->numArcs; ij++) {
		if (inBush[ij] == TRUE && bushFlows[ij] < parameters->minLinkFlow && SPtree[network->arcs[ij].head] != ij) {
			bushFlows[ij] = 0;
			inBush[ij] = FALSE;
		}
	}

	/* Find longest bush path and check reachability */
	findBushLP(network, inBush, bushOrder, LPtree, LPcosts);
	for (i = 0; i < network->numNodes; i++) {
	   if (LPcosts[i] == INFINITY && i != origin) {
         for (j = 0; j < network->numNodes; j++) {
            displayMessage(FULL_DEBUG, "%lu %f %lu %f %lu\n", j, SPcosts[j], SPtree[j], LPcosts[j], LPtree[j]);
         }
         fatalError("%lu unreachable from %lu on bush longest path\n", i, origin);
      }
	}

	/* Add shortcut arcs using strict criteria (shortcut according to both L and U labels); avoid centroid connectors */
	for (ij = 0; ij < network->numArcs; ij++) {
		if (SPcosts[network->arcs[ij].tail] + network->arcs[ij].cost < SPcosts[network->arcs[ij].head]
				 && LPcosts[network->arcs[ij].tail] + network->arcs[ij].cost < LPcosts[network->arcs[ij].head]
				 && (network->arcs[ij].tail == origin || network->arcs[ij].tail >= network->firstThroughNode))
      {
			if (inBush[ij] == FALSE) newarcs++;
			inBush[ij] = TRUE;
		}
	}

   /* If none added using strict criteria, use looser criteria (shortcut according to U labels only); avoid centroid connectors */
	if (newarcs == 0) {
		for (ij = 0; ij < network->numArcs; ij++) {
			if (LPcosts[network->arcs[ij].tail] + network->arcs[ij].cost < LPcosts[network->arcs[ij].head]
				&& (network->arcs[ij].tail == origin || network->arcs[ij].tail >= network->firstThroughNode))
         {
            inBush[ij] = TRUE;
         }
		}
	}

   /* Update topological order */
	bushTopologicalOrder(network, inBush, bushOrder);

	/* Extra descending pass to ensure consistent flow conservation */
	recalculateBushFlows(network, inBush, bushFlows, bushOrder, SPtree, parameters);
}

/*
updateFlowsB performs flow shifts on a bush, following Nie (2010).  Moving through the bush in forward topological order,
identify divergence nodes, longest/shortest path segments, and perform the Newton shift.
*/
void updateFlowsB(long origin, network_type *network, bool* inBush, long* bushOrder, double *bushFlows, algorithmBParameters_type *parameters) {
	long i, ij, curnode, diverge, SPNode, LPNode;
	bool found_diverge;
	double gap = 0, dx, dx2, der1, der2, cost1, cost2, stepsize;

	declareVector(long, mark, network->numNodes);
	declareVector(long, LPtree, network->numNodes);
	declareVector(long, SPtree, network->numNodes);
	declareVector(double, LPcosts, network->numNodes);
	declareVector(double, SPcosts, network->numNodes);
	declareVector(double, tempFlows, network->numArcs);

   /* Update longest and shortest paths */
	findBushSP(network, inBush, bushOrder, SPtree, SPcosts);
	findBushUsedLP(network, inBush, bushOrder, bushFlows, LPtree, LPcosts, parameters);

	linkedList *segment1 = createLinkedList();
	linkedList *segment2 = createLinkedList();
	linkedListElt* currentarc;

	/* Update tempFlows array; this stores the amount of flow on (i,j) NOT from this bush.
	   These flows are treated as fixed for now.  */
	for (ij = 0; ij < network->numArcs; ij++) {
		tempFlows[ij] = network->arcs[ij].flow - bushFlows[ij];
	}

	/* The mark array helps with finding divergence nodes */
	for (i = 0; i < network->numNodes; i++) mark[i] = NO_PATH_EXISTS;

	/* Ascending pass to calculate flow shifts (cf.  Nie, 2010) */
	for (i = 1; i < network->numNodes; i++) {
		curnode = bushOrder[i];
		if (SPtree[curnode] == NO_PATH_EXISTS || LPtree[curnode] == NO_PATH_EXISTS) continue;
		gap = fabs(LPcosts[curnode] - SPcosts[curnode]);
		displayMessage(FULL_DEBUG, "Gap for node %d is %f\n", curnode, gap);
		/* Can skip node if the difference between L and U labels is small, or if the SP and LP backnodes agree */
		if (gap < parameters->nodeGapTolerance || LPtree[curnode] == SPtree[curnode]) continue;

		/* Find divergence node (curnode) */
		found_diverge = FALSE; SPNode = curnode; LPNode = curnode;
		mark[curnode] = curnode;
		while (found_diverge == FALSE) {
			if (SPNode != origin) {
				SPNode = network->arcs[SPtree[SPNode]].tail;
				if (mark[SPNode] == curnode) {diverge = SPNode; found_diverge = TRUE;}
				mark[SPNode] = curnode;
			}
			if (LPNode != origin) {
				LPNode = network->arcs[LPtree[LPNode]].tail;
				if (mark[LPNode] == curnode) {diverge = LPNode; found_diverge = TRUE;}
				mark[LPNode] = curnode;
			}
		}

		/* Form the longest-path and shortest-path segments */
		dx = INFINITY; dx2 = INFINITY; der1 = 0; der2 = 0; cost1 = 0; cost2 = 0;
		ij = LPtree[curnode];
		do {
			if (bushFlows[ij] < dx) dx = bushFlows[ij];
			insertLinkedList(segment1, ij, NULL);
			LPNode = network->arcs[ij].tail;
			if (LPNode == origin || LPNode == diverge) break;
			ij = LPtree[LPNode];
		} while (TRUE);
		ij = SPtree[curnode];
		do {
			if (bushFlows[ij] < dx2) dx2 = bushFlows[ij];
			insertLinkedList(segment2, ij, NULL);
			SPNode = network->arcs[ij].tail;
			if (SPNode == origin || SPNode == diverge) break;
			ij = SPtree[SPNode];
		} while (TRUE);

		/* Perform Newton shift */
		do {
			currentarc = segment1->head;
         /* Calculate step size */
			do {
				ij = (*currentarc).value;
				displayMessage(FULL_DEBUG, "Arc (%d,%d) has cost,der %f %f\n", network->arcs[ij].tail+1, network->arcs[ij].head+1, network->arcs[ij].cost, network->arcs[ij].der);
				cost1 += network->arcs[ij].cost;
				der1 += network->arcs[ij].der;
				currentarc = (*currentarc).next;
			} while (currentarc != NULL);
			currentarc = segment2->head;
			do {
				ij = (*currentarc).value;
				displayMessage(FULL_DEBUG, "Arc (%d,%d) has cost,der %f %f\n", network->arcs[ij].tail+1, network->arcs[ij].head+1, network->arcs[ij].cost, network->arcs[ij].der);
				cost2 += network->arcs[ij].cost;
				der2 += network->arcs[ij].der;
				currentarc = (*currentarc).next;
			} while (currentarc != NULL);
			if (der1 + der2 == 0) der1 = parameters->minDerivative;
			stepsize = parameters->newtonStep * (cost1 - cost2) / (der1 + der2);
			if (stepsize < 0 && -stepsize > dx2) stepsize = -dx2; else if (stepsize > dx) stepsize = dx;

			/* Shift flow and update link performance functions */
			dx = stepsize;
			currentarc = segment1->head;
			do {
				ij = (*currentarc).value;
				bushFlows[ij] -= dx;
				displayMessage(FULL_DEBUG, "Subtracting %f from (%d,%d)\n", dx, network->arcs[ij].tail+1, network->arcs[ij].head+1);
				network->arcs[ij].flow = tempFlows[ij] + bushFlows[ij];
				if (network->arcs[ij].flow < 0) network->arcs[ij].flow = 0;
				network->arcs[ij].cost = linkCost(&network->arcs[ij], BPR);
				network->arcs[ij].der = linkCost(&network->arcs[ij], BPR_DER);
				currentarc = (*currentarc).next;
			} while (currentarc != NULL);
			currentarc = segment2->head;
			do {
				ij = (*currentarc).value;
				bushFlows[ij] += dx;
				displayMessage(FULL_DEBUG, "Adding %f to (%d,%d)\n", dx, network->arcs[ij].tail+1, network->arcs[ij].head+1);
				network->arcs[ij].flow = tempFlows[ij] + bushFlows[ij];
				if (network->arcs[ij].flow < 0) network->arcs[ij].flow = 0;
				network->arcs[ij].cost = linkCost(&network->arcs[ij], BPR);
				network->arcs[ij].der = linkCost(&network->arcs[ij], BPR_DER);
				currentarc = (*currentarc).next;
			} while (currentarc != NULL);
		} while (0); /* Can add convergence criteria here if you want more inner-inner iterations.  In this case need to reset the 'mark' array */
		clearLinkedList(segment1);
		clearLinkedList(segment2);
	}

	if (verbosity == FULL_DEBUG) { printf("Beckmann %f\n", BeckmannFunction(network)); exit(-1); }

   deleteLinkedList(segment1);
   deleteLinkedList(segment2);
	deleteVector(LPcosts);
	deleteVector(SPcosts);
	deleteVector(LPtree);
	deleteVector(SPtree);
	deleteVector(mark);
	deleteVector(tempFlows);
}

/*
initializeOriginB initializes bushes to shortest paths at free-flow
*/
void initializeOriginB(long origin, network_type *network, bool *inBush, long *bushOrder, double*bushFlows, algorithmBParameters_type *parameters) {
	long i;

	declareVector(double, label, network->numNodes);
	declareVector(arc_type *, backarc, network->numNodes);
	declareVector(double, remainingVehicles, network->numNodes);
	declareVector(long, SPtree, network->numNodes);

	arcBellmanFord(origin, label, backarc, network, DEQUE);

	for (i = 0; i < network->numArcs; i++) {
		inBush[i] = FALSE;
		bushFlows[i] = 0;
	}

   SPtree[origin] = NO_PATH_EXISTS;
	for (i = 0; i < network->numNodes; i++) {
		if (i != origin) {
			if (backarc[i] == NULL) fatalError("No path from %lu to %lu!\n", origin, i);
			inBush[ptr2arc(network, backarc[i])] = TRUE;
			SPtree[i] = ptr2arc(network, backarc[i]);
		}
	}

	bushTopologicalOrder(network, inBush, bushOrder);

	recalculateBushFlows(network, inBush, bushFlows, bushOrder, SPtree, parameters);

	deleteVector(label);
	deleteVector(backarc);
	deleteVector(remainingVehicles);
	deleteVector(SPtree);
}

/*
recalculateBushFlows processes a bush's flows to address various numerical issues without attempting to shift flow.
Assumes flows which are close to zero should actually be zero; then makes a pass in descending topological order,
allocating flow proportional to what it was before.
*/
void recalculateBushFlows(network_type *network, bool *inBush, double *bushFlows, long *bushOrder, long *SPtree, algorithmBParameters_type *parameters) {
   arcListElt *curArc;
   long i, ij, curnode;
   double outflow, inflow, delta;

	/* Wipe near-zero flows */
   for (i = 0; i < network->numArcs; i++) {
	   if (bushFlows[i] < parameters->minLinkFlow) bushFlows[i] = 0;
	}

	/* Perform a backward pass */
	for (i = network->numNodes - 1; i >= 1; i--) {
		curnode = bushOrder[i];

		/* Calculate node outflow and inflow */
		if (curnode < network->numZones) {
		   outflow = network->OD[bushOrder[0]][curnode].demand;
      } else {
         outflow = 0;
      }
      for (curArc = network->nodes[curnode].forwardStar.head; curArc != NULL; curArc = curArc->next) {
         ij = ptr2arc(network, curArc->arc);
         if (inBush[ij] == TRUE) outflow += bushFlows[ij];
      }
  		inflow = 0;
      for (curArc = network->nodes[curnode].reverseStar.head; curArc != NULL; curArc = curArc->next) {
         ij = ptr2arc(network, curArc->arc);
         if (inBush[ij] == TRUE) inflow += bushFlows[ij];
      }
		if (outflow < parameters->minLinkFlow) outflow = 0;

		/* Redistribute flow among incoming links according to proportions */
		if (outflow > 0) {
			if (inflow > 0) { /* Normal case: redistribute proportionally */
            for (curArc = network->nodes[curnode].reverseStar.head; curArc != NULL; curArc = curArc->next) {
               ij = ptr2arc(network, curArc->arc);
					if (inBush[ij] == TRUE) {
						delta = network->arcs[ij].flow - bushFlows[ij];
						bushFlows[ij] *= outflow / inflow;
						network->arcs[ij].flow = delta + bushFlows[ij];
					}
				}

			} else { /* Exceptional case: positive outflow, zero inflow.  Put all inflow on shortest path */
			   delta = network->arcs[SPtree[curnode]].flow - bushFlows[SPtree[curnode]];
				bushFlows[SPtree[curnode]] = outflow;
				network->arcs[SPtree[curnode]].flow = delta + bushFlows[SPtree[curnode]];
			}
		} else { /* Exceptional case: no outflow, zero all incoming flows */
         for (curArc = network->nodes[curnode].reverseStar.head; curArc != NULL; curArc = curArc->next) {
            ij = ptr2arc(network, curArc->arc);
            if (inBush[ij] == TRUE) {
					network->arcs[ij].flow -= bushFlows[ij];
					if (network->arcs[ij].flow < 0) network->arcs[ij].flow = 0;
					bushFlows[ij] = 0;
				}
			}
		}
	}

}

/*
bushTopologicalOrder is a specialized algorithm for finding topological order on bushes.
*/
void bushTopologicalOrder(network_type *network, bool *inBush, long *bushOrder) {
   arcListElt *curArc;
	long i, ij, j, next;

	declareVector(long, indegree, network->numNodes);
	for (i = 0; i < network->numNodes; i++) {
	   indegree[i] = 0;
	   bushOrder[i] = NO_PATH_EXISTS;
	}
	for (i = 0; i < network->numArcs; i++) if (inBush[i] == TRUE) indegree[network->arcs[i].head]++;
	queue_type LIST = createQueue(network->numNodes, network->numNodes);
	next = 0;
	for (i = 0; i < network->numNodes; i++) if (indegree[i] == 0) enQueue(&LIST, i);
	while (LIST.curelts > 0) {
		i = deQueue(&LIST);
		bushOrder[next++] = i;
      for (curArc = network->nodes[i].forwardStar.head; curArc != NULL; curArc = curArc->next) {
         ij = ptr2arc(network, curArc->arc);
			if (inBush[ij] == FALSE) continue;
			j = network->arcs[ij].head;
			indegree[j]--;
			if (indegree[j] == 0) enQueue(&LIST, j);
		}
	}
	if (next < network->numNodes) fatalError("Graph given to bushTopologicalOrder contains a cycle.");

	deleteQueue(&LIST);
	deleteVector(indegree);
}
