#ifndef WAVEFN_INIT
#define WAVEFN_INIT
#include "solver.h"

typedef enum
{
	RANDOM	
} init_types;

typedef double (initialization*)(double x, double y, double xMax, double yMax, int n);
initialization get_initalizationfn(parameters* p);

#endif