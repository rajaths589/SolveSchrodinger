#include <stdio.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_vector.h>
#include <wavefninit.h>

#include "solver.h"

inline double getX(parameters* p, int i)
{
	return -(p->xMax) + (i*p->xres);
}

inline double getY(parameters* p, int j)
{
	return -(p->yMax) + (j*p->yres);
}

int loadParameters(parameters* p, char* filename)
{
	FILE* configure_fp = fopen(filename, "r+");

	if(fopen == NULL)
		return 1;

	char line[21];
	int lsize = 20;

	p->deltaT = 0.0f;

	getline(&line, &lsize, configure_fp);
	sscanf(line, "%f", &(p->timestep));

	getline(&line, &lsize, configure_fp);
	sscanf(line, "%f %f", &(p->xMax), &(p->yMax));

	getline(&line, &lsize, configure_fp);
	sscanf(line, "%d %d", &(p->xsteps), &(p->ysteps));

	getline(&line, &lsize, configure_fp);
	sscanf(line, "%f", &(p->tol));

	getline(&line, &lsize, configure_fp);
	sscanf(line, "%f", &(p->chi));

	getline(&line, &lsize, configure_fp);
	sscanf(line, "%d", &(p->max_iter));

	getline(&line, &lsize, configure_fp);
	sscanf(line, "%d", &(p->n));

	p->xres = (2*p->xMax) / ((p->xsteps)-1);
	p->yres = (2*p->yMax) / ((p->ysteps)-1);

	p->current_iter = 0;
	p->stop = 0;

	//Define parameter decoding logic to decipher the init&op conditions.
	p->inittype = RANDOM;
	p->optype = HARMONIC;

	fclose(fp);
	return 0;
}

stateset* initializeSolver(parameters* p, double (*init_func) (double, double, double, double, int))
{
	stateset* wavefns = (stateset*) malloc(sizeof(stateset));
	wavefns->correct = 0;
	wavefns->eigenspectrum = (state**) malloc(sizeof(state*)* p->n);
	wavefns->trial_eigenspectrum = (state**) malloc(sizeof(state*)* p->n);

	wavefns->rms_nEnergy = 0.0f;
	wavefns->rms_nEnergyDelta = 0.0f;
	wavefns->t_rms_nEnergy = 0.0f;
	wavefns->t_rms_nEnergyDelta = 0.0f;
	wavefns->rms_energyExpectation = 0.0f;

	int i,j,k;

	for(i=0; i<n; i++)		
	{
		wavefns->eigenspectrum[i] = (state*) malloc(sizeof(state));
		wavefns->eigenspectrum[i]->eigenfn = gsl_matrix_complex_alloc(p->xsteps, p->ysteps);
		wavefns->eigenspectrum[i]->eigenval = 0.0f;
		wavefns->eigenspectrum[i]->normalization_energy = 0.0f;

		wavefns->trial_eigenspectrum[i] = (state*) malloc(sizeof(state));
		wavefns->trial_eigenspectrum[i]->eigenfn = gsl_matrix_complex_alloc(p->xsteps, p->ysteps);
		wavefns->trial_eigenspectrum[i]->eigenval = 0.0f;
		wavefns->trial_eigenspectrum[i]->normalization_energy = 0.0f;
	}

	for(k=0;k< p->n;k++)
	{
		for(i=0; i< p->xsteps; i++)
		{
			for(j=0; j < p->ysteps; j++)
			{
				gsl_matrix_complex_set(wavefns->eigenspectrum[k]->eigenfn, i, j, gsl_complex_rect(init_func(getX(p,i),init_func(getY(p,j),p->xMax, p->yMax,k)), 0.0));
				gsl_matrix_complex_set(wavefns->trial_eigenspectrum[k]->eigenfn, i, j, gsl_complex_rect(init_func(getX(p,i),init_func(getY(p,j),p->xMax, p->yMax,k)), 0.0));
			}
		}
	}

	return stateset;
}

void orthonormalize(stateset *s)
{
	double rms_nEnergy = 0.0f;
	double rms_nEnergyDelta = 0.0f;
	double t_rms_nEnergy = 0.0f;
	double t_rms_nEnergyDelta = 0.0f;

	gsl_matrix_complex* overlap = gsl_matrix_complex_alloc(s->n, s->n);
	int i,j;

	for(i=0; i<s->n;i++)
	{
		for(j=0;j<s->n;j++)
		{
			gsl_matrix_complex_set(overlap, i, j, dotproduct(s->eigenspectrum[i], s->eigenspectrum[j]);
		}
	}	

	gsl_matrix_complex *evec = gsl_matrix_complex_alloc(s->n,s->n);
	gsl_vector_complex *eval = gsl_vector_complex_alloc(s->n);
	gsl_eigen_hermv_workspace* w = gsl_eigen_hermv_workspace(s->n);

	gsl_eigen_hermv(overlap, eval, evec, w);
	gsl_eigen_hermv_free(w);

	gsl_matrix_complex** newvectors = (gsl_matrix_complex**) malloc(s->n * sizeof(gsl_matrix_complex*));

	for(i=0; i<s->n; i++)
	{
		newvectors[i] = gsl_matrix_complex_alloc(s->eigenspectrum[i]->eigenfn->size1, s->eigenspectrum[i]->eigenfn->size2);
		gsl_matrix_complex_set_zero(newvectors[i]);
	}

	gsl_complex eigenval;
	gsl_vector_complex_view eigenvec_i;
	gsl_matrix_complex tempMatrix = gsl_matrix_complex_alloc(s->eigenspectrum[0]->eigenfn->size1, s->eigenspectrum[0]->eigenfn->size2);

	for(i=0;i<s->n;i++)
	{
		eigenval = gsl_vector_complex_get(eval, i);
		eigenvec_i = gsl_matrix_complex_column(evec, i);

		rms_nEnergy += pow((log(GSL_REAL(eigenval))/(-2.0*p->timestep)),2);
		rms_nEnergyDelta += pow(s->eigenspectrum[i]->normalization_energy-(log(GSL_REAL(eigenval))/(-2.0*p->timestep)),2);
		s->eigenspectrum[i]->normalization_energy = eigenval;

		for(j=0; j<s->n;j++)
		{
			gsl_matrix_complex_memcpy(tempMatrix, s->eigenspectrum[j]->eigenfn);
			gsl_matrix_complex_scale(tempMatrix, gsl_vector_complex_get(eigenvec_i,j));
			gsl_matrix_complex_add(newvectors[i], tempMatrix);			
		}
		gsl_matrix_complex_scale(newvectors[i], gsl_complex_inverse(gsl_complex_sqrt_real(eigenval)));
	}

	for(i=0; i<s->n; i++)
	{
		gsl_matrix_complex_memcpy(s->eigenspectrum[i]->eigenfn, newvectors[i]);		
		gsl_matrix_complex_set_zero(newvectors[i]);
	}

	rms_nEnergy /= s->n;
	rms_nEnergy = sqrt(rms_nEnergy);
	s->rms_nEnergy = rms_nEnergy;
	rms_nEnergyDelta /= s->n;
	rms_nEnergyDelta = sqrt(rms_nEnergyDelta);
	s->rms_nEnergyDelta = rms_nEnergyDelta;

	for(i=0; i<s->n;i++)
	{
		for(j=0;j<s->n;j++)
		{
			gsl_matrix_complex_set(overlap, i, j, dotproduct(s->trial_eigenspectrum[i], s->trial_eigenspectrum[j]);
		}
	}

	w = gsl_eigen_hermv_workspace(s->n);
	gsl_eigen_hermv(overlap, eval, evec, w);
	gsl_eigen_hermv_free(w);

	for(i=0;i<s->n;i++)
	{
		eigenval = gsl_vector_complex_get(eval, i);
		eigenvec_i = gsl_matrix_complex_column(evec, i);

		t_rms_nEnergy += pow((log(GSL_REAL(eigenval))/(-2.0*p->timestep)),2);
		t_rms_nEnergyDelta += pow(s->trial_eigenspectrum[i]->normalization_energy-(log(GSL_REAL(eigenval))/(-2.0*p->chi*p->timestep)),2);
		s->trial_eigenspectrum[i]->normalization_energy = eigenval;

		for(j=0; j<s->n;j++)
		{
			gsl_matrix_complex_memcpy(tempMatrix, s->trial_eigenspectrum[j]->eigenfn);
			gsl_matrix_complex_scale(tempMatrix, gsl_vector_complex_get(eigenvec_i,j));
			gsl_matrix_complex_add(newvectors[i], tempMatrix);			
		}
		gsl_matrix_complex_scale(newvectors[i], gsl_complex_inverse(gsl_complex_sqrt_real(eigenval)));
	}	

	for(i=0; i<s->n; i++)
	{
		gsl_matrix_complex_memcpy(s->trial_eigenspectrum[i]->eigenfn, newvectors[i]);
		gsl_matrix_complex_free(newvectors[i]);		
	}

	t_rms_nEnergy /= s->n;
	t_rms_nEnergy = sqrt(t_rms_nEnergy);
	t_rms_nEnergyDelta /= s->n;
	t_rms_nEnergyDelta = sqrt(t_rms_nEnergyDelta);
	s->t_rms_nEnergy = t_rms_nEnergy;
	s->t_rms_nEnergyDelta = t_rms_nEnergyDelta;

	free(newvectors);
	gsl_matrix_complex_free(tempMatrix);
	gsl_vector_complex_free(eigenvec_i);

}

void advanceImaginaryTime(parameters* p, stateset* s, operators* ops, fftw_plan* plans)
{
	p->current_iter += 1;	

	copy(s);
	
	ops->evolution(*s, p, ops->potential, ops->kinetic, plans);
	orthonormalize(*s);

	if((s->rms_nEnergy < s->t_rms_nEnergy) || (s->rms_nEnergyDelta < s->t_rms_nEnergyDelta))
	{
		s->chi = sqrt(s->chi);
	}
	else
	{
		temp = s->eigenspectrum;
		s->eigenspectrum = s->trial_eigenspectrum;
		s->trial_eigenspectrum = temp;
		s->timestep *= s->chi;
		s->chi = pow(s->chi,2);
		s->rms_nEnergy = s->t_rms_nEnergy;
		s->rms_nEnergyDelta = s->t_rms_nEnergyDelta;
	}

	ops->hamiltonian(s,p,ops->potential, ops->kinetic, plans);
	if(s->rms_nEnergy <= s->tol*s->rms_energyExpectation)
	{
		s->stop = 1;
	}
}


static fftw_plan* create_fftwplan(int xsize, int ysize, gsl_matrix_complex* m)
{
	fftw_complex* f = (fftw_complex*)m->data;
	fftw_plan* fp = (fftw_plan*) malloc(2*sizeof(fftw_plan));
	fp[0] = fftw_plan_dft_2d(xsize, ysize, f, f, FFTW_FORWARD,  FFTW_MESAURE);
	fp[1] = fftw_plan_dft_2d(xsize, ysize, f, f, FFTW_BACKWARD,  FFTW_MESAURE);

	return fp;
}

int solveSchrodingerEquation(parameters* p, char* outputFilename)
{	
	initialization initn = get_initalizationfn(p);
	stateset* s = initializeSolver(p, initn);
	operators* ops = init_ops(p);
	fftw_plan* fpw = create_fftwplan(p->xsteps, p->ysteps, s->eigenspectrum[0]->eigenfn);
	while(p->current_iter <= p->max_iter && p->stop==0)
	{
		advanceImaginaryTime(p, s, ops, fpw);
	}
	
	//Write a better print function
	int i;
	for(i=0;i<s->n;i++)
	{
		printf("%d\n",s->eigenspectrum[i]->eigenval);
	}

	free(ops);
	fftw_destroy_plan(fpw[0]);
	fftw_destroy_plan(fpw[1]);
	free(fpw);
}

int main()
{
	parameters* p = (parameters*) malloc(sizeof(parameters));
	loadParameters(p, "config.txt");
	solveSchrodingerEquation(p,"solution.txt");
}