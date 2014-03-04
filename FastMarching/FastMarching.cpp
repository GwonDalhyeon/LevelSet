#include "MarchingForward.h"
#include <limits>

void main()
{
	double* testFunction = new  double[9];
	testFunction[0] = 0.0;testFunction[1]=0.0;testFunction[2] = 0.0;testFunction[3]=2.0;
	testFunction[4] = 3.0;testFunction[5]=1.0;testFunction[6] = numeric_limits<double>::infinity();testFunction[7]=numeric_limits<double>::infinity();
	testFunction[8] = numeric_limits<double>::infinity();//testFunction[9]=574.5;testFunction[10] = -3.2;testFunction[11]=-9.5;
	struct LevelSet* initial =new LevelSet(3, 1.0,testFunction);
	initial->pointStatus[3]=0;
	initial->pointStatus[4]=0;
	initial->pointStatus[5]=0;
	
}