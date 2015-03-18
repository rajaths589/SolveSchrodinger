#include "wavefninit.h"
#include <time.h>
#include <stdlib.h>

double random_init(double x, double y, double xMax, double yMax, int n)
{
	srand(time(NULL));
	return (double)rand()/15.0;
}

initialization get_initalizationfn(parameters* p)
{
	//implement decoding logic later
	return random_init;
}