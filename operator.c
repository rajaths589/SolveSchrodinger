#include "operator.h"
#include "math.h"
#include "solver.h"
#include <fftw3.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

double harmonic_potential(double x, double y)
{
	return ((x*x) + (y*y))/2;
}

double kinetic_2d(double kx, double ky)
{
	return ((kx*k_x) + (ky*ky)/2);
}

static void transform_to_momentumspace(stateset* s, fftw_plan* fp)
{
	fftw_complex* in;
	int k;
	for(k=0; k<s->n; k++)
	{
		in = (fftw_complex*) (s->eigenspectrum[k]->eigenfn->data);
		fftw_execute_dft(fp[0], in , in);
		gsl_matrix_complex_scale(s->eigenspectrum[k]->eigenfn, gsl_complex_rect((1.0/(s->eigenspectrum[k]->eigenfn->size1*s->eigenspectrum[k]->eigenfn->size2)),0));

		in = (fftw_complex*) (s->trial_eigenspectrum[k]->eigenfn->data);
		fftw_execute_dft(fp[0], in , in);
		gsl_matrix_complex_scale(s->trial_eigenspectrum[k]->eigenfn, gsl_complex_rect((1.0/(s->trial_eigenspectrum[k]->eigenfn->size1*s->trial_eigenspectrum[k]->eigenfn->size2)),0));
	}
}

//by default eigenvectors are in position representation. Do not use this function if they are not already transformed into momentum space.
static void transform_to_positionspace(stateset* s, fftw_plan* fp)
{
	fftw_complex* in;
	int k;
	for(k=0; k<s->n; k++)
	{
		in = (fftw_complex*) (s->eigenspectrum[k]->eigenfn->data);
		fftw_execute_dft(fp[1], in , in);		

		in = (fftw_complex*) (s->trial_eigenspectrum[k]->eigenfn->data);
		fftw_execute_dft(fp[1], in , in);		
	}
}

static void transform_to_momentumspace_i(stateset* s, fftw_plan* fp, int i)
{
	fftw_complex* in;
	in = (fftw_complex*) (s->eigenspectrum[i]->eigenfn->data);
	fftw_execute_dft(fp[0], in , in);
	gsl_matrix_complex_scale(s->eigenspectrum[i]->eigenfn, gsl_complex_rect((1.0/(s->eigenspectrum[k]->eigenfn->size1*s->eigenspectrum[k]->eigenfn->size2)),0));
}

static void transform_to_positionspace_i(stateset* s, fftw_plan* fp, int i)
{
	fftw_complex* in;
	in = (fftw_complex*) (s->eigenspectrum[i]->eigenfn->data);
	fftw_execute_dft(fp[1], in , in);	
}

static double reduceto2pi(double k)
{
	if(k>2*M_PI)
		return k-2*M_PI;

	return k;
}

// |delV|^2 is calculated
static double finitedifference_gradient(op potential, double x, double y, double xres. double yres)
{
	double a = (potential(x+(xres/2), y) - potential(x-(xres/2), y))/xres;
	double b = (potential(x, y+(yres/2)) - potential(x, y-(yres/2))/yres;

	return sqrt(a*a + b*b);
}

void evolution_operator_4(stateset* s, parameters* p, op potential, op kinetic, fftw_plan* fp)
{
	int i, k, j;


	//multiply e^(-exp*V(x,y)/6)
	for(k=0;k<s->n;k++)
	{
		for(i=0;i<p->xsteps;i++)
		{
			for(j=0;j<y->steps;j++)
			{
				gsl_matrix_complex_set(s->eigenspectrum[k]->eigenfn, i, j, gsl_complex_mul_real(gsl_matrix_complex_get(s->eigenspectrum[k]->eigenfn,i,j),exp((-(p->timestep/6.0))*potential(getX(i),getY(j)))));
				gsl_matrix_complex_set(s->trial_eigenspectrum[k]->eigenfn, i, j, gsl_complex_mul_real(gsl_matrix_complex_get(s->trial_eigenspectrum[k]->eigenfn,i,j),exp((-(p->timestep*p->chi/6.0))*potential(getX(i),getY(j)))));
			}
		}
	}

	//fourier transform
	transform_to_momentumspace(s,fp);

	int kx, ky;
	for(k=0;k<s->n;k++)
	{		
		for(i=0; i<p->xsteps; i++)
		{
			for(j=0; j<p->ysteps; j++)
			{
				kx = i*2*M_PI/p->xsteps;
				ky = j*2*M_PI/p->ysteps;

				kx = reduceto2pi(kx);
				ky = reduceto2pi(ky);

				gsl_matrix_complex_set(s->eigenspectrum[k]->eigenfn, i, j, gsl_complex_mul_real(gsl_matrix_complex_get(s->eigenspectrum[k]->eigenfn,i,j),exp((-(p->timestep/2.0))*kinetic(kx, ky));
				gsl_matrix_complex_set(s->trial_eigenspectrum[k]->eigenfn, i, j, gsl_complex_mul_real(gsl_matrix_complex_get(s->trial_eigenspectrum[k]->eigenfn,i,j),exp((-(p->timestep*p->chi/2.0))*kinetic(kx, ky))));
			}
		}
	}

	transform_to_positionspace(s, fp);

	double pot, x, y;
	for(k=0;k<s->n;k++)
	{
		for(i=0;i<p->xsteps;i++)
		{
			for(j=0;j<y->steps;j++)
			{
				x = getX(i);
				y = getY(j);
				pot = potential(x,y) + (finitedifference_gradient(potential, x, y, p->xres, p->yres)*(p->timestep*p->timestep)/48.0);
				gsl_matrix_complex_set(s->eigenspectrum[k]->eigenfn, i, j, gsl_complex_mul_real(gsl_matrix_complex_get(s->eigenspectrum[k]->eigenfn,i,j),exp((-(p->timestep*2.0/3.0))*pot)));
				
				pot = potential(x,y) + (finitedifference_gradient(potential, x, y, p->xres, p->yres)*(p->timestep*p->timestep*p->chi*p->chi)/48.0);
				gsl_matrix_complex_set(s->trial_eigenspectrum[k]->eigenfn, i, j, gsl_complex_mul_real(gsl_matrix_complex_get(s->trial_eigenspectrum[k]->eigenfn,i,j),exp((-(p->timestep*2.0*p->chi/3.0))*pot)));
			}
		}
	}

	transform_to_momentumspace(s,fp);
	for(k=0;k<s->n;k++)
	{		
		for(i=0; i<p->xsteps; i++)
		{
			for(j=0; j<p->ysteps; j++)
			{
				kx = i*2*M_PI/p->xsteps;
				ky = j*2*M_PI/p->ysteps;

				kx = reduceto2pi(kx);
				ky = reduceto2pi(ky);

				gsl_matrix_complex_set(s->eigenspectrum[k]->eigenfn, i, j, gsl_complex_mul_real(gsl_matrix_complex_get(s->eigenspectrum[k]->eigenfn,i,j),exp((-(p->timestep/2.0))*kinetic(kx, ky));
				gsl_matrix_complex_set(s->trial_eigenspectrum[k]->eigenfn, i, j, gsl_complex_mul_real(gsl_matrix_complex_get(s->trial_eigenspectrum[k]->eigenfn,i,j),exp((-(p->timestep*p->chi/2.0))*kinetic(kx, ky))));
			}
		}
	}

	transform_to_positionspace(s,fp);
	for(k=0;k<s->n;k++)
	{
		for(i=0;i<p->xsteps;i++)
		{
			for(j=0;j<y->steps;j++)
			{
				gsl_matrix_complex_set(s->eigenspectrum[k]->eigenfn, i, j, gsl_complex_mul_real(gsl_matrix_complex_get(s->eigenspectrum[k]->eigenfn,i,j),exp((-(p->timestep/6.0))*potential(getX(i),getY(j)))));
				gsl_matrix_complex_set(s->trial_eigenspectrum[k]->eigenfn, i, j, gsl_complex_mul_real(gsl_matrix_complex_get(s->trial_eigenspectrum[k]->eigenfn,i,j),exp((-(p->timestep*p->chi/6.0))*potential(getX(i),getY(j)))));
			}
		}
	}
}

void hamiltonian_4(stateset* s, parameters* p, op potential, op kinetic, fftw_plan* fp)
{
	double delH = 0.0;
	gsl_complex energy;
	gsl_complex t;
	//gsl_matrix_complex* temp = gsl_matrix_complex_alloc(p->xsteps, p->ysteps);
	int i,j,k;

	for(k=0; k<s->n; k++)
	{		
		energy = gsl_complex_rect(0,0);
		for(i=0; i<p->xsteps; i++)
		{
			for(j=0; j<p->ysteps; j++)
			{
				t = gsl_matrix_complex_get(s->eigenspectrum[k]->eigenfn, i, j);
				t = gsl_complex_mul(t, gsl_complex_rect(potential(getX(p,i),getY(p,j)),0.0));
				t = gsl_complex_mul(t, gsl_complex_conjugate(gsl_complex_get(s->eigenspectrum[k]->eigenfn, i, j)));
				energy = gsl_complex_add(t,energy);
			}
		}

		transform_to_momentumspace_i(s, fp, k);
	
		for(i=0; i<p->xsteps; i++)
		{
			for(j=0; j<p->ysteps; j++)
			{
				kx = i*2*M_PI/p->xsteps;
				ky = j*2*M_PI/p->ysteps;

				kx = reduceto2pi(kx);
				ky = reduceto2pi(ky);

				t = gsl_matrix_complex_get(s->eigenspectrum[k]->eigenfn, i, j);
				t = gsl_complex_mul(t, gsl_complex_rect(kinetic(kx, ky),0.0));
				t = gsl_complex_mul(t, gsl_complex_conjugate(gsl_complex_get(s->eigenspectrum[k]->eigenfn, i, j)));
				energy = gsl_complex_add(t,energy);
			}
		}

		transform_to_positionspace_i(s, fp, k);

		s->eigenspectrum[k]->eigenval = GSL_REAL(energy);
		delH += pow(s->eigenspectrum[k]->eigenval-s->eigenspectrum[k]->normalization_energy, 2);
	}

	delH = delH / s->n;
	delH = sqrt(delH);

	s->rms_energyExpectation = delH;
}

operators* init_ops(parameters* p)
{
	operators* op = (operators*) malloc(sizof(operators));
	
	//conditions can be added.
	//I'm seting it to the implemented functions
	op->potential = harmonic_potential;
	op->kinetic = kinetic_2d;
	op->hamiltonian = hamiltonian_4;
	op->evolution = evolution_operator_4;

	return op;
}
