#ifndef WAVEFN_INIT
#define WAVEFN_INIT
#include "parameters.h"

typedef double (*initialization)(double , double , double , double , int);
initialization get_initalizationfn(parameters* p);

#endif