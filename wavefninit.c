#include "wavefninit.h"
#include <time.h>
#include <stdlib.h>

double random_init(double x, double y, double xMax, double yMax, int n)
{
	int r = rand()%100;
	return r/100.00;
}

initialization get_initalizationfn(parameters* p)
{
	//implement decoding logic later
	srand(time(NULL));
	return random_init;
}
