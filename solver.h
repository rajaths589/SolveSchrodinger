#ifndef ITP_SOLVER
#define ITP_SOLVER
#include <gsl/gsl_complex.h>

#include "state.h"
#include "wavefninit.h"

typedef struct parameters {
	double timestep;
	double deltaT;
	
	double xres;	// distance between two discrete points along x-axis, i.e. delta_x
	double yres;	// distance between two discrete points along y-axis, i.e. delta_y
	
	double xMax;	// boundary assumed from (-xMax, xMax)
	double yMax;	// boundary assumed from (-yMax, yMax)
	
	int xsteps;		// number of x points on x-axis
	int ysteps;		// number of y points on y-axis

	double tol;		// a factor in the range 0.1 ... 1 determined emperically : refer to papers
	double chi;		// factor by which timestep is changed	

	int max_iter;
	int current_iter;
	int n;			// number of wavefunctions to solve
	int stop;

	init_types inittype;
	operator_types optype;
} parameters;

int loadParameters(char* filename, parameters* p);
stateset* initializeSolver(parameters* p, double (*init_func) (double, double, double, double, int));
int solveSchrodingerEquation(parameters* p, char* outputFilename);
void advanceImaginaryTime(parameters* p, stateset* s, operators* ops, fftw_plan* plans);
inline double getX(parameters* p, int i);
inline double gety(parameters* p, int j);
void orthonormalize(stateset *s);

#endif