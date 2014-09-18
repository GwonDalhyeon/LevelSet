#include <iostream>
#include <cmath>
#include <fstream>
#include "memoryuse.h"
using namespace std;


struct Phi2D
{
	double x, y;
	double phi;
	double tempPhiValue;
	double originPhi;

	double xVelocity;
	double yVelocity;
	double velocityField;

	double k1,k2,k3;

	int xIndex, yIndex;
	Phi2D* PhiLeft;
	Phi2D* PhiRight;
	Phi2D* PhiTop;
	Phi2D* PhiBottom;

	Phi2D* PhiNext;
	Phi2D* PhiBefore;
	Phi2D();
};



Phi2D::Phi2D()
{
	PhiLeft   = NULL;
	PhiRight  = NULL;
	PhiTop    = NULL;
	PhiBottom = NULL;

	PhiNext   = NULL;
	PhiBefore = NULL;
}


void deleteAllPhi(Phi2D* inputPhi)
{
	if (inputPhi != NULL)
	{
		deleteAllPhi(inputPhi ->PhiNext);
		delete inputPhi;
	}
}


struct ZeroLevelSetPoint
{
	double x,y;

	ZeroLevelSetPoint* PointBefore;
	ZeroLevelSetPoint* PointNext;
	
	ZeroLevelSetPoint();
};

ZeroLevelSetPoint::ZeroLevelSetPoint()
{
	PointBefore = NULL;
	PointNext   = NULL;
}

void deleteAllZeroLevelSetPoint(ZeroLevelSetPoint* inputPoint)
{
	if (inputPoint != NULL)
	{
		deleteAllZeroLevelSetPoint(inputPoint->PointNext);
		delete inputPoint;
	}
}

