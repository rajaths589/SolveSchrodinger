#ifndef ITP_STATESET
#define ITP_STATESET

#include <gsl/gsl_complex_math.h>

typedef struct state {
	gsl_matrix_complex* eigenfn;
	double eigenval;
	double normalization_energy;
} state;


typedef struct stateset {
	state** eigenspectrum;
	state** trial_eigenspectrum;
	int n;	// number of eigenstates
	int correct;	// refers to which set of states has the more correct wavevectors

	double rms_nEnergy;
	double rms_nEnergyDelta;
	double t_rms_nEnergy;
	double t_rms_nEnergyDelta;
	double rms_energyExpectation;
} stateset;

void normalize(state* eigenfn);
gsl_complex dotproduct(state* stateA, state* stateB);

stateset* createcopy(stateset* set);

#endif