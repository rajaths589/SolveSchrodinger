#ifndef ITP_OPS
#define ITP_OPS

#include "state.h"
#include "solver.h"
#include <fftw3.h>

typedef double (*op) (double, double);
typedef void (*ev_op) (stateset* , parameters, op, op, fftw_plan* fp);
typedef void (*H_op) (stateset*, op, op, fftw_plan* fp);

typedef enum
{
	HARMONIC
} operator_types;

typedef struct operator
{
	ev_op evolution;
	H_op hamiltonian;
	op potential;
	op kinetic;
} operators;

operators* init_ops(parameters* p);
double harmonic_potential(double x, double y);
double kinetic_2d(double kx, double ky);
void hamiltonian_4(stateset* s, parameters* p, op potential, op kinetic, fftw_plan* fp);
void evolution_operator_4(stateset* s, parameters* p, op potential, op kinetic, fftw_plan* fp);

#endif