#ifndef CONVEXCOMBINATIONS_H
#define CONVEXCOMBINATIONS_H

#include <limits.h>
#include <stdio.h>
#include <time.h>
#include "fileio.h"
#include "networks.h"
#include "tap.h"
#include "utils.h"

typedef enum {
   BISECTION,
   NEWTON
} lineSearch_type;

typedef struct {
   gap_type gapFunction;
   double   convergenceGap;
   double   maxTime;
   long     maxIterations;
   bool     warmStart;
   double   lineSearchIterations;
   lineSearch_type lineSearch;
} FrankWolfeParameters_type;

typedef struct {
   gap_type gapFunction;
   double   convergenceGap;
   double   maxTime;
   long     maxIterations;
   bool     warmStart;
} MSAparameters_type;

FrankWolfeParameters_type initializeFrankWolfeParameters();
void FrankWolfe(network_type *network, FrankWolfeParameters_type *parameters);

double FrankWolfeLineSearch(network_type *network, FrankWolfeParameters_type *parameters, double *oldFlows, double *SPFlows);
double FrankWolfeBisection(network_type *network, FrankWolfeParameters_type *parameters, double *oldFlows, double *SPFlows);
double FrankWolfeNewton(network_type *network, FrankWolfeParameters_type *parameters, double *oldFlows, double *SPFlows);

MSAparameters_type initializeMSAparameters();
void MSA(network_type *network, MSAparameters_type *parameters);

void getSPFlows(network_type *network, double *SPFlows, network_type *SPtree, double *tempFlows);
void initializeFlows(network_type *network, network_type *SPtree, double *newArcFlows);
void loadOrigin(long origin, network_type *network, network_type *SPtree, double *newArcFlows);

#endif
