#ifndef ITP_STATESET
#define ITP_STATESET

#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>

typedef struct state {
	gsl_matrix_complex* eigenfn;
	double eigenval;
	double normalization_energy;
} state;


typedef struct stateset {
	state** eigenspectrum;
	state** trial_eigenspectrum;
	int n;	// number of eigenstates	

	double rms_nEnergy;
	double rms_nEnergyDelta;
	double t_rms_nEnergy;
	double t_rms_nEnergyDelta;
	double rms_energyExpectation;
} stateset;

gsl_complex dotproduct(state* stateA, state* stateB);

stateset* createcopy(stateset* set);

#endif