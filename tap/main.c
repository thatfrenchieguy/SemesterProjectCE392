/*
   This program performs user-equilibrium assignment.  The network file and OD matrix
   are specified as the second and third arguments of the readOBANetwork function.
   The program outputs the resulting link flows and costs to the file specified in the
   second argument of the writeLinkOutputFile function.

   The program uses the OBA file format, as documented on Hillel Bar-Gera's website:
      https://github.com/bstabler/TransportationNetworks
   This site also includes a number of test networks of different size.  If you are unfamiliar
   with this file format, please

	Three algorithms are provided: MSA, Frank-Wolfe, and Algorithm B.  The current main
	function reads a network then runs these three algorithms in succession.  The parameters
	for each algorithm are set in a struct passed to the equilibrium algorithm function.
	Default values are provided in the initialization function.

	BE SURE TO SET AT LEAST ONE OF THE THREE CONVERGENCE CRITERIA (BASED ON GAP, MAXIMUM
	TIME, OR MAXIMUM NUMBER OF ITERATIONS) OR ELSE THE PROGRAM WILL RUN FOREVER.  The program
	terminates as soon as any one of these convergence criteria is satisfied.

	FOR ALGORITHM B
      Bparameters.gapFunction -- how you want to check gap convergence.  Allowable options
                                 are AEC, RELATIVE_GAP_1, and RELATIVE_GAP_2 (There are
                                 slightly different ways of defining relative gap;
                                 RELATIVE_GAP_1 is the one I taught in class.)
                                 Default value: RELATIVE_GAP_1
      Bparameters.convergenceGap -- Program terminates when the gap function falls below
                                    this value.  If no value is set, the program will
                                    continue running until another termination criterion is
                                    met.
      Bparameters.maxTime -- Maximum time in seconds before the program terminates.
                             If no value is set, the program will continue running
                             until one of the other termination criteria is met.
      Bparameters.maxIterations -- Maximum number of iterations before the program terminates.
                                   If no value is set, the program will continue running
                                   until one of the other termination criteria is met.
      Bparameters.innerIterations -- Number of bush flow adjustments to make before updating
                                     bush links.  Default value = 20
      Bparameters.minLinkFlow -- When bushes are updated, any link flow value smaller than
                                 this is set to zero.  Set to a small positive value to avoid
                                 problems due to round-off errors.  Default value = 1e-14
      Bparameters.minDerivative -- Lower bound on derivative calculations to avoid division
                                   by zero in Newton's method.  Default value = 1e-6
      Bparameters.newtonStep -- Step size for use in Newton's method.  Default value = 1
      Bparameters.nodeGapTolerance -- Minimum difference between longest and shortest path
                                      travel times to a node for shifting flow.  Default
                                      value = 1e-14

   FOR FRANK-WOLFE
      FWparameters.gapFunction -- same as for Algorithm B
      FWparameters.convergenceGap -- same as for Algorithm B
      FWparameters.maxTime -- same as for Algorithm B
      FWparameters.maxIterations -- same as for Algorithm B
      FWparameters.warmStart -- Equals TRUE iff you already have a feasible solution in the
                                network structure that you want to use as an initial solution.
                                Default value is FALSE, and FW begins by finding an initial solution
                                based on all-or-nothing assignment to free-flow shortest paths.
      FWparameters.lineSearch -- The method used to find the step size lambda.  Allowable choices
                                 are BISECTION and NEWTON.  Default value is BISECTION.
      FWparameters.lineSearchIterations -- The number of iterations to perform when determining
                                           the step size lambda (for either BISECTION or NEWTON).
                                           Default value is 11, which is more appropriate for
                                           BISECTION.  (NEWTON generally requires fewer.)

   FOR MSA:
      MSAparameters.gapFunction -- same as for Algorithm B
      MSAparameters.convergenceGap -- same as for Algorithm B
      MSAparameters.maxTime -- same as for Algorithm B
      MSAparameters.maxIterations -- same as for Algorithm B
      MSAparameters.warmStart -- same as for Frank-Wolfe

   A few notes on reading the source code:
      bush.c and bush.h contain code for Algorithm B.
      convexcombination.c/h contain code for MSA and Frank-Wolfe.
      fileio.c/h contain code for reading network files and formatting output.
      networks.c/h contain code for basic network algorithms (shortest path, topological order) and
         working with network data structures.  Look here for definitions of the network struct.
      tap.c/h contain code common to all traffic assignment algorithms -- BPR functions, gap
         calculations, TSTT/SPTT, etc.  Look here for linkCost.
      datastructures.c/h contain code implementing standard data structures -- linked lists, queues,
         and binary heaps -- and memory (de)allocation.  Look here for declareScalar/declareVector/declareMatrix,
         newScalar, deleteVector, etc.
      utils.c/h contain code for assorted tasks -- displaying messages, memory (de)allocation, timing, and so on.
         Look here for displayMessage, warning, fatalError,
*/
#include "main.h"

int main() {
   verbosity = FULL_NOTIFICATIONS; /* verbosity is a global variable affecting how much text is output; see utils.c */
   #ifdef DEBUG_MODE
      debugFile = openFile("full_log.txt", "w");
   #endif

   network_type *network = newScalar(network_type);
   algorithmBParameters_type Bparameters = initializeAlgorithmBParameters();
   FrankWolfeParameters_type FWparameters = initializeFrankWolfeParameters();
   MSAparameters_type MSAparameters = initializeMSAparameters();

   readOBANetwork(network, "SiouxFalls_net.txt", "SiouxFalls_trips.txt");
   makeStronglyConnectedNetwork(network); /* Check connectivity */

   displayMessage(FULL_NOTIFICATIONS, "Starting Frank-Wolfe...\n");
   FWparameters.maxIterations = 11;
   FWparameters.gapFunction = RELATIVE_GAP_1;
   FrankWolfe(network, &FWparameters);

   //displayMessage(FULL_NOTIFICATIONS, "Starting Algorithm B...\n");
   //Bparameters.convergenceGap = 1e-13;
   //Bparameters.gapFunction = AEC;
   //AlgorithmB(network, &Bparameters);

   displayMessage(FULL_NOTIFICATIONS, "Equilibrium done, press a key to see link flows.\n");
   waitForKey();
   displayNetwork(LOW_NOTIFICATIONS, network);

   deleteNetwork(network);

   #ifdef DEBUG_MODE
      fclose(debugFile);
   #endif

   return EXIT_SUCCESS;
}
