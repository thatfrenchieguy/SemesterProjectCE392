/*
   This file contains (1) implementations of general-purpose network algorithms and (2) supporting infrastructure
   for the network data structure.

   Regarding (1), this file includes two variants of label-correcting shortest path algorithms (BellmanFord
   and arcBellmanFord), one returning the previous node in the shortest path and the other returning the previous
   link in the shortest path.  A heap-based implementation of Dijkstra's (label-setting) shortest path algorithm
   is also included, as is a connectivity checker (search), topological order finder, a routine creating
   artificial links to ensure the network is strongly connected (every node reachable from every other), and
   a procedure for sorting a list of links into forward star order.

   Regarding (2), this file contains code for displaying network data in human-readable format, and implementations
   of linked lists for links and paths.

*/

#include "networks.h"

/*
BellmanFord is the one-to-all shortest path algorithm taught in class, with the use of a scan-eligible list to handle
cycles in the network.  The label and backnode arguments are pointers to arrays which will contain the L and q labels
upon completion of the algorithm.  origin gives the starting point for these paths.  The network structure is passed
via reference (pointer) to reduce overhead.  Argument q takes one of three values, specifying how the scan eligible
list is maintained: FIFO (first-in first-out order, so nodes are scanned in the order added to the list), LIFO
(last-in first-out, nodes scanned in reverse order of being added to the list), and DEQUE (double-ended queue; nodes
generally go to the end of the queue, as in FIFO, but if a node has been scanned before it is added to the front end
as in LIFO).
*/
void BellmanFord(long origin, double *label, long *backnode, network_type *network, queueDiscipline q) {
	long j, curnode;
	arcListElt *i;
	double tempLabel;
   /* My implementation of a queue uses a "circular" array to store elements;
      queue is empty iff the "read" and "write" pointers are identical.
      See datastructures.c for more details. */
	queue_type SEL = createQueue(network->numNodes, network->numNodes);

	/* Initialize */
	for (j = 0; j < network->numNodes; j++) {
		label[j] = INFINITY;
		backnode[j] = NO_PATH_EXISTS;
	}
	label[origin] = 0;
	enQueue(&SEL, origin);

   /* Iterate */
	while (SEL.readptr != SEL.writeptr) {  /* See comment above about read and write pointers */
		curnode = deQueue(&SEL);
		for (i = network->nodes[curnode].forwardStar.head; i != NULL; i = i->next) {
			tempLabel = label[curnode] + i->arc->cost;
			j = i->arc->head;
			if (tempLabel < label[j]) { /* We have found a better path to node j */
				label[j] = tempLabel;
				backnode[j] = curnode;
				if (j >= network->firstThroughNode) { /* Ensure we do not use centroids/centroid connectors as "shortcuts" */
					switch (q) {
						case FIFO: 	enQueue(&SEL, j);	break;
						case DEQUE:
							switch (SEL.history[j]) {
								case NEVER_IN_QUEUE: enQueue(&SEL, j); break;
								case WAS_IN_QUEUE: frontQueue(&SEL, j); break;
							}
							break;
						case LIFO:  frontQueue(&SEL, j); break;
						default: fatalError("bellmanFord: Unsupported queue structure");
						break;
					}
				}
			}
		}
	}

	/* Clean up */
	deleteQueue(&SEL);
}

/*
arcBellmanFord is virtually identical to BellmanFord, except that it returns pointers to the previous *link* in the
shortest paths, rather than the index of the previous *node*.  (backarc is an array of pointers to links).  For explanations
of other arguments, see description of BellmanFord.
*/
void arcBellmanFord(long origin, double *label, arc_type **backarc, network_type *network, queueDiscipline q) {
	long j, curnode;
	arcListElt *i;
	double tempLabel;
	queue_type SEL = createQueue(network->numNodes, network->numNodes);

	for (j = 0; j < network->numNodes; j++) {
		label[j] = INFINITY;
		backarc[j] = NULL;
	}
	label[origin] = 0;
	enQueue(&SEL, origin);

	while (SEL.readptr != SEL.writeptr) {
		curnode = deQueue(&SEL);
		for (i = network->nodes[curnode].forwardStar.head; i != NULL; i = i->next) {
			tempLabel = label[curnode] + i->arc->cost;
			j = i->arc->head;
			if (tempLabel < label[j]) {
				label[j] = tempLabel;
				backarc[j] = i->arc;
				if (j >= network->firstThroughNode) {
					switch (q) {
						case FIFO: 	enQueue(&SEL, j);	break;
						case DEQUE:
							switch (SEL.history[j]) {
								case NEVER_IN_QUEUE: enQueue(&SEL, j); break;
								case WAS_IN_QUEUE: frontQueue(&SEL, j); break;
							}
							break;
						case LIFO:
						default: fatalError("arcBellmanFord: Unsupported queue structure");
						break;
					}
				}
			}
		}
	}

	deleteQueue(&SEL);

}


/*
heapDijsktra is an implementation of Dijkstra's algorithm using a binary heap.  This is a one-to-all
shortest path algorithm, starting at 'origin'.  *label and *backnode are arrays which return the cost
and previous node labels upon completion.  *network is a pointer to the underlying network data structure,
passing by reference is faster.
*/
void heapDijkstra(long origin, double *label, long *backnode, network_type *network) {

	long j;
	arcListElt *i;
	long curnode;

   /* Initialize heap */
	double tempLabel;
	heap_type *dijkstraHeap = createHeap(network->numNodes, network->numNodes);

   /* Initialize Dijkstra's */
	for (j = 0; j < network->numNodes; j++) {
		dijkstraHeap->valueFn[j] = INFINITY; /* valueFn in the heap stores the cost labels */
		backnode[j] = NO_PATH_EXISTS;
	}

   /* Now iterate until the heap is empty */
	insertHeap(dijkstraHeap, origin, 0);
	while (dijkstraHeap->last > 0) {
		curnode = findMinHeap(dijkstraHeap);
		deleteMinHeap(dijkstraHeap);
		for (i = network->nodes[curnode].forwardStar.head; i != NULL; i = i->next) {
			j = i->arc->head;
			tempLabel = dijkstraHeap->valueFn[curnode] + i->arc->cost;
			if (tempLabel < dijkstraHeap->valueFn[j]) {
				backnode[j] = curnode;
				if (j < network->firstThroughNode) {
					dijkstraHeap->valueFn[j] = tempLabel;
					continue;
				}
				if (dijkstraHeap->valueFn[j] == INFINITY)
					insertHeap(dijkstraHeap, j, tempLabel);
				else
					decreaseKey(dijkstraHeap, j, tempLabel);
			}
		}
	}

   /* Now copy labels to return, and clean up memory */
	memcpy(label + 1, dijkstraHeap->valueFn + 1, sizeof(double) * (network->numNodes));
	deleteHeap(dijkstraHeap);
}


/*
finalizeNetwork: After adding the links and nodes to the network struct, this function
generates the forward and reverse star lists, and also calculates fixed costs for each
link (which only need to be done this one time, rather than at each call to linkCost).
*/
void finalizeNetwork(network_type *network) {
	long i;

	for (i = 0; i < network->numNodes; i++) {
		initializeArcList(&(network->nodes[i].forwardStar));
		initializeArcList(&(network->nodes[i].reverseStar));
	}
	for (i = 0; i < network->numArcs; i++) {
		insertArcList(&(network->nodes[network->arcs[i].tail].forwardStar), &(network->arcs[i]), network->nodes[network->arcs[i].tail].forwardStar.tail);
		insertArcList(&(network->nodes[network->arcs[i].head].reverseStar), &(network->arcs[i]), network->nodes[network->arcs[i].head].reverseStar.tail);
      network->arcs[i].fixedCost = network->arcs[i].length * network->distanceFactor + network->arcs[i].toll * network->tollFactor;
	}
}

/*
search: Given an initial node (the origin argument), performs a search to identify all nodes reachable from origin, or from
which origin can be reached, depending on the argument 'd' (FORWARD = nodes reachable from origin; REVERSE = nodes from which
origin can be reached).  Argument 'q' indicates the search order (FIFO = breadth-first search, LIFO = depth-first search).
Upon termination, the arrays order and backnode are returned.  backnode indicates the previous/next node on the path from/to
origin (depending on d); if backnode[i] is the symbolic constant NO_PATH_EXISTS then node i is not connected.  The order
array indicates the order in which nodes are found.
*/
void search(long origin, long* order, long *backnode, network_type *network, queueDiscipline q, direction_type d) {
	long i, j, next;
	arcListElt *curarc;

	/* Initialize; any node for which backnode[i] remains at NO_PATH_EXISTS is not connected from/to origin */
	for(i = 0; i < network->numNodes; i++) {
		backnode[i] = NO_PATH_EXISTS;
	}
	backnode[origin] = 0;
	next = 1;
	order[origin] = next;

   /* List of visited nodes is maintained as a queue with discipline q */
	queue_type LIST = createQueue(network->numNodes, network->numNodes);
	enQueue(&LIST, origin);

	/* This code uses a circular queue implementation; queue is empty iff readptr and writeptr are identical */
	while (LIST.readptr != LIST.writeptr) {
		i = deQueue(&LIST);
		/* Identify the proper list (forward or reverse) ... */
		switch (d) {
			case FORWARD: curarc = network->nodes[i].forwardStar.head; break;
			case REVERSE: curarc = network->nodes[i].reverseStar.head; break;
			default: fatalError("Unknown direction in search."); break;
		}
		/* ...and now iterate through all its elements */
		while (curarc != NULL) {
			switch (d) {
				case FORWARD: j = curarc->arc->head; break;
				case REVERSE: j = curarc->arc->tail; break;
			}
			if (backnode[j] == NO_PATH_EXISTS) { /* Is admissible; arc discovers a new node */
				backnode[j] = i;
				displayMessage(FULL_DEBUG, "Next node found is %d-%d\n", j / 4, j % 4);
				order[j] = ++next;
				if (j >= network->firstThroughNode) {
					switch (q) {
						case FIFO: enQueue(&LIST, j); break;
						case LIFO: frontQueue(&LIST, j); break;
						case DEQUE:
							switch (LIST.history[j]) {
								case NEVER_IN_QUEUE: enQueue(&LIST, j); break;
								case WAS_IN_QUEUE: frontQueue(&LIST, j); break;
							}
					default: fatalError("Unsupported queue type in search."); break;
					}
				}
			}
			curarc = curarc->next;
		}
	}
	deleteQueue(&LIST);
	if (verbosity >= FULL_DEBUG) waitForKey();
	return;
}

/*
topologicalOrder takes a network and finds a topological order (if one exists).
Argument *network is a pointer to the network struct, arguments *sequence and
*order are arrays which this function fills in.  order is the topological order
itself (maps each node to its index); sequence gives the inverse function (index -> node)
which is helpful if you want to iterate over the nodes in order.
*/
void topologicalOrder(network_type *network, long* sequence, long* order) {
	long i, j, cur, next;
	arcListElt *ij;

	declareVector(long, indegree, network->numNodes);
	for (i = 0; i < network->numNodes; i++) {
	   indegree[i] = 0;
	   order[i] = NO_PATH_EXISTS;
	}
	for (i = 0; i < network->numArcs; i++) indegree[network->arcs[i].head]++;
	queue_type LIST = createQueue(network->numNodes, network->numNodes);
	next = 0; cur = 0;
	for (i = 0; i < network->numNodes; i++) if (indegree[i] == 0) enQueue(&LIST, i);
	while (LIST.curelts > 0) {
		i = deQueue(&LIST);
		order[i] = ++next; sequence[cur++] = i;
		for (ij = network->nodes[i].forwardStar.head; ij != NULL; ij = ij->next) {
			j = ij->arc->head;
			indegree[j]--;
			if (indegree[j] == 0) enQueue(&LIST, j);
		}
	}
	if (next < network->numNodes) { displayNetwork(FULL_NOTIFICATIONS, network); fatalError("Graph given to topologicalOrder contains a cycle"); }

	deleteVector(indegree);
	deleteQueue(&LIST);
}

/*
forwardStarOrder is a comparison function; a pointer to this function can be passed
to qsort or other sorting routines.  In forward star order, a link precedes another
if its tail node has a lower index.  This function implements a tiebreaking rule based
on the head node.
*/
int forwardStarOrder(const void *arc1, const void *arc2) {
	arc_type first = *(arc_type *)arc1;
	arc_type second = *(arc_type *)arc2;
	if (first.tail < second.tail) {
		return -1;
	} else if (first.tail > second.tail) {
		return 1;
	} else if (first.head < second.head) {
		return -1;
	} else if (first.head > second.head) {
		return 1;
	} else {
		return 0;
	}
}

/*
makeStronglyConnectedNetwork creates artificial links, if needed, to ensure that every link
in the network is reachable from every other.  It does this by calling 'search' in both the
forward and reverse directions from an arbitrarily chosen node (the one with the highest index),
call this node *.  Artificial links are created to and from * and any nodes which cannot be found from
these searches.  As a result, the network becomes strongly connected: any two nodes are at least
reachable through the path i -> * -> j.  The freeFlowTime (and various other parameters) on these
links are set equal to the symbolic constant ARTIFICIAL, which should be set high enough that
the links will not be used.
*/
void makeStronglyConnectedNetwork(network_type *network) {
	long i, j;
	declareVector(long, order, network->numNodes);
	declareVector(long, backnode, network->numNodes);
	declareVector(long, forwardnode, network->numNodes);

   /* Run forward and reverse searches to see what links must be created */
	search(network->numNodes - 1, order, backnode, network, FIFO, FORWARD);
	search(network->numNodes - 1, order, forwardnode, network, FIFO, REVERSE);

	long newArcs = 0;
	for (i = 0; i < network->numNodes; i++) {
		if (backnode[i] == NO_PATH_EXISTS) newArcs++;
		if (forwardnode[i] == NO_PATH_EXISTS) newArcs++;
	}
	if (newArcs == 0) { /* Nothing to do, network is already strongly connected, so clean up/return */
		deleteVector(order);
		deleteVector(backnode);
		deleteVector(forwardnode);
		return;
	}

   /* Create new arc vector, with all the old ones plus the new artificial links */
	displayMessage(FULL_NOTIFICATIONS, "Warning: Creating %d artifical arcs to ensure strong connectivity.\n", newArcs);
	declareVector(arc_type, newArcVector, network->numArcs + newArcs);
	for (i = 0; i < network->numArcs; i++) {
		newArcVector[i] = network->arcs[i];
	}
	for (j = 0; j < network->numNodes; j++) {
		if (backnode[j] == NO_PATH_EXISTS) {
			displayMessage(DEBUG, "Creating (%d,%d)\n", network->numNodes, j + 1);
			newArcVector[i].tail = network->numNodes - 1;
			newArcVector[i].head = j;
			newArcVector[i].alpha = 0;
			newArcVector[i].beta = 1;
			newArcVector[i].flow = 0;
			newArcVector[i].capacity = ARTIFICIAL;
			newArcVector[i].length = ARTIFICIAL;
			newArcVector[i].toll = ARTIFICIAL;
			newArcVector[i].freeFlowTime = ARTIFICIAL;
			i++;
		}
		if (forwardnode[j] == NO_PATH_EXISTS) {
			displayMessage(DEBUG, "Creating (%d,%d)\n", j + 1, network->numNodes);
			newArcVector[i].tail = j;
			newArcVector[i].head = network->numNodes - 1;
			newArcVector[i].alpha = 0;
			newArcVector[i].beta = 1;
			newArcVector[i].flow = 0;
			newArcVector[i].capacity = ARTIFICIAL;
			newArcVector[i].length = ARTIFICIAL;
			newArcVector[i].toll = ARTIFICIAL;
			newArcVector[i].freeFlowTime = ARTIFICIAL;
			i++;
		}
	}
	network->numArcs += newArcs;
	deleteVector(network->arcs);
	network->arcs = newArcVector;

   /* Regenerate forward/reverse star lists */
	for (j = 0; j < network->numNodes; j++) {
		clearArcList(&(network->nodes[j].forwardStar));
		clearArcList(&(network->nodes[j].reverseStar));
	}
	finalizeNetwork(network);

   /* Clean up and return */
	deleteVector(order);
	deleteVector(backnode);
	deleteVector(forwardnode);
}

/*
quicksortDestination: Uses Darel Finley's QuickSort implementation to sort
an array of 'nodes' according to the array of 'costs' (typically distances from
the origin)
*/
#define MAX_QUICKSORT_LEVELS 1000
void quicksortDestinations(long *nodes, double *costs, int elements) {
	long piv;
	int beg[MAX_QUICKSORT_LEVELS], end[MAX_QUICKSORT_LEVELS], i, L, R ;
	for (i = 0; i < elements; i++) nodes[i] = i;
	i = 0;
	beg[0]=1; end[0]=elements+1;
	while (i>=0) {
		L=beg[i]; R=end[i]-1;
		if (L<R) {
			piv=nodes[L]; if (i==MAX_QUICKSORT_LEVELS-1) { fatalError("quicksortdestinations: TOO MANY LEVELS"); }
			while (L<R) {
				while (costs[nodes[R]] >= costs[piv] && L<R) R--; if (L<R) nodes[L++]=nodes[R];
				while (costs[nodes[L]] <= costs[piv] && L<R) L++; if (L<R) nodes[R--]=nodes[L]; }
			nodes[L]=piv; beg[i+1]=L+1; end[i+1]=end[i]; end[i++]=L;
		} else {
			i--;
		}
	}
}

/*
ptr2arc converts a pointer to an arc to the index number for that arc; to do this, the network struct
needs to be passed along with the arc pointer.
*/
int ptr2arc(network_type *network, arc_type *arcptr) {
	return (int) (arcptr - network->arcs);
}

////////////////////////////////////
// Custom data structures for TAP //
////////////////////////////////////

/*
displayNetwork prints network data in human-readable format.  minVerbosity is used to control
whether anything needs to be printed.
*/

void displayNetwork(int minVerbosity, network_type *network) {
	long i;
	displayMessage(minVerbosity, "Network has %d nodes and %d arcs\n", network->numNodes, network->numArcs);
	displayMessage(minVerbosity, "Arc data: ID, tail, head, flow, cost, der\n");
	for (i = 0; i < network->numArcs; i++) displayMessage(minVerbosity, "%ld (%ld,%ld) %f %f %f\n", i, network->arcs[i].tail + 1, network->arcs[i].head + 1, network->arcs[i].flow, network->arcs[i].cost, network->arcs[i].der);
}

/*
deleteNetwork deallocates any memory assigned to a network struct.
*/
void deleteNetwork(network_type *network) {
   int i;
   for (i = 0; i < network->numNodes; i++) {
      clearArcList(&(network->nodes[i].forwardStar));
      clearArcList(&(network->nodes[i].reverseStar));
   }

   deleteMatrix(network->OD, network->numZones);
   deleteVector(network->nodes);
   deleteVector(network->arcs);
   deleteScalar(network);
}

/*
The functions below implement doubly-linked lists of arcs.
You probably don't need to poke around here too much.
*/
arcList *createArcList() {
	declareScalar(arcList, newdll);
	initializeArcList(newdll);
	return newdll;
}

void initializeArcList(arcList *list) {
	list->head = NULL;
	list->tail = NULL;
	list->size = 0;
}

arcListElt *insertArcList(arcList *list, arc_type *value, arcListElt *after) {

	declareScalar(arcListElt, newNode);

	newNode->arc = value;
	if (after != NULL) {
		newNode->prev = after;
		newNode->next = after->next;
		if (list->tail != after) newNode->next->prev = newNode; else list->tail = newNode;
		after->next = newNode;
	} else {
		newNode->prev = NULL;
		newNode->next = list->head;
		if (list->tail != after) newNode->next->prev = newNode; else list->tail = newNode;
		list->head = newNode;
	}
	list->size++;
	return newNode;
}

void clearArcList(arcList *list) {
	while (list->head != NULL)
		deleteArcListElt(list, list->tail);
}

void deleteArcList(arcList *list) {
	clearArcList(list);
	deleteScalar(list);
}

void deleteArcListElt(arcList *list, arcListElt *elt) {
	if (list->tail != elt) {
		if (list->head != elt) elt->prev->next = elt->next; else list->head = elt->next;
		elt->next->prev = elt->prev;
	} else {
		list->tail = elt->prev;
		if (list->head != elt) elt->prev->next = elt->next; else list->head = elt->next;
	}
	list->size--;
	deleteScalar(elt);
}

void displayArcList(arcList *list) {
	arcListElt *curnode = list->head;
	printf("Start of the list: %p\n", (void *)list->head);
	while (curnode != NULL) {
		printf("%p (%ld,%ld) %p %p\n", (void *)curnode, curnode->arc->tail, curnode->arc->head,
													(void *)curnode->prev, (void *)curnode->next);
		curnode = (*curnode).next;
	}
	printf("End of the list: %p\n", (void *)list->tail);
}

/*
The functions below implement paths as special instances of arc doubly-linked lists.
You probably don't need to poke around here either.
*/

path_type *createPath() {
	declareScalar(path_type, path);
	path->arcs = createArcList();
	path->cost = 0;
	path->der = 0;
	return path;
}

void deletePath(path_type *path) {
	deleteArcList(path->arcs);
	deleteScalar(path);
}

void displayPath(path_type *path) {
	printf("Path cost and derivative %f %f\n", path->cost, path->der);
	displayArcList(path->arcs);
}

void displayPathCompact(int minVerbosity, path_type *path) {
	arcListElt *curArc = path->arcs->head;
	displayMessage(minVerbosity, "[");
	if (curArc != NULL) displayMessage(minVerbosity, "%ld", curArc->arc->tail);
	while (curArc != NULL) {
		displayMessage(minVerbosity, ",%ld", curArc->arc->head);
		curArc = curArc->next;
	}
	displayMessage(minVerbosity, "]");
}

bool comparePaths(path_type *path1, path_type *path2) {
	arcListElt *arc1 = path1->arcs->head, *arc2 = path2->arcs->head;
	for (arc1 = path1->arcs->head, arc2 = path2->arcs->head; !(arc1 == NULL && arc2 == NULL); arc1 = arc1->next, arc2 = arc2->next) {
		if (arc1 == NULL || arc2 == NULL) return FALSE; // Unequal path lengths
		if (arc1->arc != arc2->arc) return FALSE; // Non-matching arcs
	}
	return TRUE;
}

/*
pathSets are collections of paths.
*/

pathSet *createPathSet() {
	declareScalar(pathSet, newdll);
	initializePathSet(newdll);
	return newdll;
}

void initializePathSet(pathSet *list) {
	list->head = NULL;
	list->tail = NULL;
	list->numPaths = 0;
}

bool pathsEqual(arcList *path1, arcList *path2) {
	arcListElt *curarc1, *curarc2;
	curarc1 = path1->head;
	curarc2 = path2->head;
	while (curarc1 != NULL && curarc2 != NULL) {
		if (curarc1->arc->head != curarc2->arc->head) return FALSE;
		if (curarc1->arc->tail != curarc2->arc->tail) return FALSE;
		curarc1 = curarc1->next;
		curarc2 = curarc2->next;
	}
	return TRUE;
}

pathSetElt *insertPathSet(pathSet *list, path_type *value, pathSetElt *after) {
	declareScalar(pathSetElt, newNode);
	newNode->path = value;
	if (after != NULL) {
		newNode->prev = after;
		newNode->next = after->next;
		if (list->tail != after) newNode->next->prev = newNode; else list->tail = newNode;
		after->next = newNode;
	} else {
		newNode->prev = NULL;
		newNode->next = list->head;
		if (list->tail != after) newNode->next->prev = newNode; else list->tail = newNode;
		list->head = newNode;
	}
	list->numPaths++;
	return newNode;
}

void clearPathSet(pathSet *list) {
	while (list->head != NULL)
		deletePathSetElt(list, list->tail);
}

void deletePathSet(pathSet *list) {
	clearPathSet(list);
	deleteScalar(list);
}


void deletePathSetElt(pathSet *list, pathSetElt *elt) {
#ifdef DEBUG
	clearArcList(elt->arcs);
	deleteScalar(elt->arcs);
	if (list->tail != elt) {
		if (list->head != elt) elt->prev->next = elt->next; else list->head = elt->next;
		elt->next->prev = elt->prev;
	} else {
		list->tail = elt->prev;
		if (list->head != elt) elt->prev->next = elt->next; else list->head = elt->next;
	}
	deleteScalar(elt);
	list->numPaths--;
#endif
}


void displayPathSet(pathSet *list) {
	printf("BEGINNING PATH SET DISPLAY\n");
	pathSetElt *curnode = list->head;
	while (curnode != NULL) {
		displayPath(curnode->path);
		curnode = curnode->next;
	}
	printf("END OF PATH SET DISPLAY\n");
	waitForKey();
}

