#include "state.h"
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

void normalize(state* eigenfn)
{

}

//computes <A|B>
gsl_complex dotproduct(state* stateA, state* stateB)
{
	gsl_complex result = gsl_complex_rect(0.0, 0.0);

	int i,j;

	for(i=0; i<stateA->eigenfn.size1; i++)
	{
		for(j=0; j<stateA->eigenfn.size2; j++)
		{
			result = gsl_complex_add(result, (gsl_complex_mul(gsl_complex_conjugate(gsl_matrix_complex_get(stateA->eigenfn,i,j)),gsl_matrix_complex_get(stateB->eigenfn,i,j))));
		}
	}

	return result;
}

void normalize(state* state)
{
	gsl_complex norm = dotproduct(state, state);

	norm = gsl_complex_sqrt(norm);
	gsl_matrix_complex_scale(state->eigenfn, norm);
}

void copy(stateset* set)
{
	
	if(set->correct == 0)
	{		
		for(k=0;k< p->n;k++)
		{
			gsl_matrix_complex_memcpy(set->trial_eigenspectrum[k]->eigenfn, set->eigenspectrum[k]->eigenfn);
		}
	}
	else
	{
		for(k=0;k< p->n;k++)
		{
			gsl_matrix_complex_memcpy(set->eigenspectrum[k]->eigenfn, set->trial_eigenspectrum[k]->eigenfn);
		}	
	}	
}