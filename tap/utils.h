#ifndef UTILS_H
#define UTILS_H

#include <string.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define DEBUG_MODE /* Uncomment this line to echo output to log file. */

#define IS_MISSING -1
#define STRING_SIZE 9999

#define PAUSE_ON_ERROR FALSE
#define PAUSE_ON_WARNING FALSE


/*
Standard units: feet, seconds
Multiplying a quantity by these values will convert it to standard units
Dividing a quantity by these values will convert it from standard units
*/
#define HOURS      3600.0
#define MINUTES    60.0
#define SECONDS    1.0
#define MILES      5280.0
#define KILOMETERS 3280.839895
#define METERS     3.280839895
#define FEET       1.0
#define INCHES     0.083333333

#define min(x,y)	( ((x)<(y)) ? (x) : (y) )
#define max(x,y)	( ((x)>(y)) ? (x) : (y) )
#define swap(a,b) SWAP(&a, &b, sizeof(a))
#define round2int(x)  (int)((x) < 0 ? ((x) - 0.5) : ((x) + 0.5))
#define round2long(x)  (long)((x) < 0 ? ((x) - 0.5) : ((x) + 0.5))

#ifndef __cplusplus
typedef enum {
	FALSE,
	TRUE
} bool;
#endif

#ifdef DEBUG_MODE
	char debugFileName[STRING_SIZE];
	FILE *debugFile;
#endif


enum { /* Verbosity levels for status messages */
	NOTHING,
	LOW_NOTIFICATIONS,
	MEDIUM_NOTIFICATIONS,
	FULL_NOTIFICATIONS,
	DEBUG,
	FULL_DEBUG
};


int verbosity;


#ifdef MEMCHECK
extern long memcheck_numScalars, memcheck_numVectors, memcheck_numMatrices;
#endif

void waitForKey();
void SWAP(void *a, void *b, int size);
double updateElapsedTime(clock_t startTime, double *elapsedTime);

/*********************
 ** Status messages **
 *********************/

void displayMessage(int minVerbosity, char *format, ...);
void fatalError(char *format, ...);
void warning(int minVerbosity, char *format, ...);

#endif
