#include <iostream>
#include <cmath>
#include <fstream>
#include "memoryuse.h"
using namespace std;

struct Phi1D
{
	double x;
	double phi;

	struct Phi1D* PhiLeft;
	struct Phi1D* PhiRight;

	Phi1D();
};

struct Phi2D
{
	double x, y;
	double phi;
	double tempPhiValue;
	double originPhi;

	double k1,k2,k3;

	int xIndex, yIndex;
	struct Phi2D* PhiLeft;
	struct Phi2D* PhiRight;
	struct Phi2D* PhiTop;
	struct Phi2D* PhiBottom;

	struct Phi2D* PhiNext;
	struct Phi2D* PhiBefore;
	Phi2D();
};

Phi1D::Phi1D()
{
	PhiLeft = NULL;
	PhiRight = NULL;
}

Phi2D::Phi2D()
{
	PhiLeft   = NULL;
	PhiRight  = NULL;
	PhiTop    = NULL;
	PhiBottom = NULL;

	PhiNext   = NULL;
	PhiBefore = NULL;
}

void deleteAllPhi(Phi1D* phiHead)
{
	if (phiHead!=NULL)
	{
		deleteAllPhi(phiHead->PhiRight);
		delete phiHead;
	}
}

void deleteAllPhi(Phi2D* phiHead)
{
	if (phiHead != NULL)
	{
		deleteAllPhi(phiHead ->PhiNext);
		delete phiHead;
	}
}