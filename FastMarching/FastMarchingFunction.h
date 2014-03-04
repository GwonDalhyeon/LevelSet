#include "FastMarchingLevelset.h"
#include <cmath>
#include <algorithm>
#include <stdio.h>

double minmod(double a, double b)
{
	if (a*b<0)
	{
		return 0;
	}
	else if (abs(a)>=abs(b))
	{
		return b;
	}
	else
	{
		return a;
	}
}

double LevelSet ::derivation2x(int i,int j, double* A)
{
	double result;

	if (i>0 && i<(xNum-1))
	{
		result = A[i-1 +j*xNum]-2.0*A[i+j*xNum]+A[i+1+j*xNum];
	}
	else if (i==0)
	{
		result = A[i+j*xNum]-2.0*A[i+1+j*xNum]+A[i+2+j*xNum];
	}
	else if (i==(xNum-1))
	{
		result = A[i+j*xNum]-2.0*A[i-1+j*xNum]+A[i-2+j*xNum];
	}
	else
	{
		result = 0;
	}

	result = result/(dx*dx);
	return result;

}

double LevelSet ::derivation2y(int i,int j, double* A)
{
	double result;

	if (j>0 && j<(yNum-1))
	{
		result = A[i+(j-1)*xNum]-2.0*A[i+j*xNum]+A[i+(j+1)*xNum];
	}
	else if (j==0)
	{
		result = A[i+j*xNum]-2.0*A[i+(j+1)*xNum]+A[i+(j+2)*xNum];
	}
	else if (j==(yNum-1))
	{
		result = A[i+j*xNum]-2.0*A[i+(j-1)*xNum]+A[i+(j-2)*xNum];
	}
	else
	{
		result = 0;
	}

	result = result/(dy*dy);
	return result;
}


double LevelSet ::derivation1x_plus(int i,int j, double* A)
{
	double result;

	// calculate D+xPSI_ij
	if (i<xNum-1)
	{
		result = (A[i+1+j*xNum] -A[i+j*xNum])/dx;// - dx/2.0*minmod(derivation2x(i,j,A),derivation2x(i+1,j,A));
	}
	else
	{
		result = (A[i+j*xNum] -A[i-1+j*xNum])/dx;// - dx/2.0*minmod(derivation2x(i,j,A),derivation2x(i+1,j,A));

	}
	return result;
}


double LevelSet:: derivation1y_plus(int i,int j, double* A)
{
	double result;

	// calculate D+xPSI_ij
	if (j<yNum-1)
	{
		result = (A[i+(j+1)*xNum]-A[i+j*xNum])/dy;// - dy/2.0*minmod(derivation2y(i,j,A),derivation2y(i,j+1,A));
	}
	else
	{
		result = (A[i+j*xNum]-A[i+(j-1)*xNum])/dy;// - dy/2.0*minmod(derivation2y(i,j,A),derivation2y(i,j+1,A));

	}
	return result;
}


double LevelSet :: derivation1x_minus(int i,int j, double* A)
{
	double result;

	// calculate D+xPSI_ij
	if (i>0)
	{
		result = (A[i+j*xNum]-A[i-1+j*xNum])/dx;// + dx/2.0*minmod(derivation2x(i,j,A),derivation2x(i-1,j,A));
	}
	else
	{
		result = (A[i+1+j*xNum]-A[i+j*xNum])/dx;// +dx/2.0*minmod(derivation2x(i,j,A),derivation2x(i-1,j,A));
	}
	return result;
}

double LevelSet :: derivation1y_minus(int i,int j, double* A)
{
	double result;

	// calculate D+xPSI_ij
	if (j>0)
	{
		result = (A[i+j*xNum]-A[i+(j-1)*xNum])/dy;// + dy/2.0*minmod(derivation2y(i,j,A),derivation2y(i,j-1,A));
	}
	else
	{
		result = (A[i+(j+1)*xNum]-A[i+j*xNum])/dy;// + dy/2.0*minmod(derivation2y(i,j,A),derivation2y(i,j-1,A));
	}

	return result;
}

int extractMinNarrowBand(struct LevelSet* inputLevelSet)
{
	int result;
	int narrowBandCount = count (inputLevelSet->pointStatus, inputLevelSet->pointStatus+inputLevelSet->xNum*inputLevelSet->yNum, 0);

	int* tempIndex = new int(narrowBandCount);
	double* tempPhi = new double(narrowBandCount);

	int tempPlace=0;
	for (int i = 0; i < inputLevelSet->xNum*inputLevelSet->yNum; i++)
	{
		if (inputLevelSet->pointStatus[i] ==0)
		{
			
			*(tempIndex + tempPlace) = i;
			*(tempPhi + tempPlace) = inputLevelSet->Phi[i];
			tempPlace++;
		}
	}

	struct HeapsortTree* tempHeap = HeapInitial(tempPhi,narrowBandCount);
	struct HeapsortTree* sortedNarrowBand = HeapExtract(tempHeap, narrowBandCount);

	return result = tempIndex[sortedNarrowBand->indexArray[0]];
}




