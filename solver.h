#ifndef ITP_SOLVER
#define ITP_SOLVER
#include <gsl/gsl_complex.h>
#include <fftw3.h>

#include "state.h"
#include "parameters.h"
#include "operator.h"

int loadParameters(parameters* p, char* filename);
stateset* initializeSolver(parameters* p, double (*init_func) (double, double, double, double, int));
int solveSchrodingerEquation(parameters* p, char* outputFilename);
void advanceImaginaryTime(parameters* p, stateset* s, operators* ops, fftw_plan* plans);
void orthonormalize(stateset *s, parameters* p);

#endif