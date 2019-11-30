#ifndef NETWORKS_H
#define NETWORKS_H

#include <math.h>
#include <string.h>
#include "datastructures.h"

#define NO_PATH_EXISTS -1
#define ARTIFICIAL 99999 /* Value used for costs, etc. on artificial links */


typedef enum {
	FORWARD,
	REVERSE
} direction_type;

typedef struct {
	long   	tail;
	long   	head;
	double 	flow;
	double 	cost;
	double 	der;
	double 	freeFlowTime;
	double 	capacity;
	double 	alpha;
	double 	beta;
	double 	toll;
	double 	speedLimit;
	double 	length;
	double	fixedCost; // Reflects toll and distance
	int		linkType;
} arc_type;



#define arcListElt struct AL
arcListElt {
	arc_type 	*arc;
	arcListElt	*prev;
	arcListElt	*next;
};

typedef struct {
	arcListElt	*head;
	arcListElt 	*tail;
	int			size;
} arcList;

typedef struct {
	arcList		*arcs;
	double		cost;
	double		der;
} path_type;


#define pathSetElt struct PSE
pathSetElt	{
	path_type 	*path;
	pathSetElt	*prev;
	pathSetElt	*next;
};

typedef struct {
	pathSetElt *head;
	pathSetElt *tail;
	long			numPaths;
} pathSet;

typedef struct {
	arcList	forwardStar;
	arcList	reverseStar;
} node_type;

typedef struct {
	double demand;
	double cost;
} od_type;

typedef struct {
	node_type*	nodes;
	arc_type* 	arcs;
	od_type** 	OD;
	long 			numNodes;
	long 			numArcs;
	long 			numZones;
	long 			firstThroughNode;
	double 		totalODFlow;
	double 		tollFactor;
	double 		distanceFactor;
	double 		beckmann;
	double 		beckmannLB;
} network_type;

void BellmanFord(long origin, double *label, long *backnode, network_type *network, queueDiscipline q);
void arcBellmanFord(long origin, double *label, arc_type **backarc, network_type *network, queueDiscipline q);
void heapDijkstra(long origin, double *label, long *backnode, network_type *network);
void finalizeNetwork(network_type *network);
void makeStronglyConnectedNetwork(network_type *network);
void quicksortDestinations(long *nodes, double *costs, int elements);
void search(long origin, long* order, long *backnode, network_type *network, queueDiscipline q, direction_type d);
void topologicalOrder(network_type *network, long* sequence, long* order);
void findPrimaryLink(network_type *network, arc_type *arc);
int forwardStarOrder(const void *arc1, const void *arc2);
int ptr2arc(network_type *network, arc_type *arcptr);


/////////////////////////////////
// Custom linked lists for TAP //
/////////////////////////////////

void displayNetwork(int minVerbosity, network_type *network);
void deleteNetwork(network_type *network);

arcList *createArcList();
void initializeArcList(arcList *list);
arcListElt *insertArcList(arcList *list, arc_type *value, arcListElt *after);
void clearArcList(arcList *list);
void deleteArcList(arcList *list);
void deleteArcListElt(arcList *list, arcListElt *elt);
void displayArcList(arcList *list);
path_type *createPath();
void displayPath(path_type *path);
void displayPathCompact(int minVerbosity, path_type *path);
bool comparePaths(path_type *path1, path_type *path2);

pathSet *createPathSet();
void initializePathSet(pathSet *list);
bool pathsEqual(arcList *path1, arcList *path2);
pathSetElt *insertPathSet(pathSet *list, path_type *value, pathSetElt *after);
void clearPathSet(pathSet *list);
void deletePathSet(pathSet *list);
void deletePathSetElt(pathSet *list, pathSetElt *elt);
void displayPathSet(pathSet *list);



#endif

