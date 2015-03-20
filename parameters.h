#ifndef PARAMETERS_DEFN_ITP
#define PARAMETERS_DEFN_ITP

typedef enum
{
	RANDOM	
} init_types;

typedef enum
{
	HARMONIC
} operator_types;

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

inline double getX(parameters* p, int i);
inline double gety(parameters* p, int j);

#endif