#include "shape.h"

void Shape2D::reinitialTVDRK3()
{
	Phi2D* tempPhi = PhiHead;
	double tempDxPlusPhi;
	double tempDxMinusPhi;
	double tempDyPlusPhi;
	double tempDyMinusPhi;

	while (tempPhi!=NULL)
	{
		tempPhi->originPhi = tempPhi->phi;

		tempDxPlusPhi  = wenoApproximationXPlus(tempPhi);
		tempDxMinusPhi = wenoApproximationXMinus(tempPhi);
		tempDyPlusPhi  = wenoApproximationYPlus(tempPhi);
		tempDyMinusPhi = wenoApproximationYMinus(tempPhi);

		tempPhi->k1 = -sign2(tempPhi->phi)*deltaT*reinitialGodonov(tempDxPlusPhi,tempDxMinusPhi,tempDyPlusPhi,tempDyMinusPhi, tempPhi->phi);
		tempPhi = tempPhi->PhiNext;
	}
	
	tempPhi = PhiHead;
	while (tempPhi!=NULL)
	{
		tempPhi->phi = tempPhi->originPhi + tempPhi->k1;
		tempPhi = tempPhi->PhiNext;
	}

	tempPhi = PhiHead;
	while (tempPhi!=NULL)
	{
		tempDxPlusPhi  = wenoApproximationXPlus(tempPhi);
		tempDxMinusPhi = wenoApproximationXMinus(tempPhi);
		tempDyPlusPhi  = wenoApproximationYPlus(tempPhi);
		tempDyMinusPhi = wenoApproximationYMinus(tempPhi);

		tempPhi->k2 = -sign2(tempPhi->phi)*deltaT*reinitialGodonov(tempDxPlusPhi,tempDxMinusPhi,tempDyPlusPhi,tempDyMinusPhi, tempPhi->phi);
		tempPhi = tempPhi->PhiNext;
	}
	
	tempPhi = PhiHead;
	while (tempPhi!=NULL)
	{
		tempPhi->phi = 3.0/4.0*tempPhi->originPhi + 1.0/4.0*tempPhi->phi + 1.0/4.0*tempPhi->k2;
		tempPhi = tempPhi->PhiNext;
	}

	tempPhi = PhiHead;
	while (tempPhi!=NULL)
	{
		tempDxPlusPhi  = wenoApproximationXPlus(tempPhi);
		tempDxMinusPhi = wenoApproximationXMinus(tempPhi);
		tempDyPlusPhi  = wenoApproximationYPlus(tempPhi);
		tempDyMinusPhi = wenoApproximationYMinus(tempPhi);

		tempPhi->k3 = -sign2(tempPhi->phi)*deltaT*reinitialGodonov(tempDxPlusPhi,tempDxMinusPhi,tempDyPlusPhi,tempDyMinusPhi, tempPhi->phi);
		tempPhi = tempPhi->PhiNext;
	}
	
	tempPhi = PhiHead;
	while (tempPhi!=NULL)
	{
		tempPhi->phi = 1.0/3.0*tempPhi->originPhi + 2.0/3.0*tempPhi->phi + 2.0/3.0*tempPhi->k3;

		tempPhi = tempPhi->PhiNext;
	}
}


void Shape2D::propagatingTVDRK3()
{
	Phi2D* tempPhi = PhiHead;
	double tempDxPlusPhi;
	double tempDxMinusPhi;
	double tempDyPlusPhi;
	double tempDyMinusPhi;

	while (tempPhi!=NULL)
	{
		tempPhi->originPhi = tempPhi->phi;

		tempDxPlusPhi  = wenoApproximationXPlus(tempPhi);
		tempDxMinusPhi = wenoApproximationXMinus(tempPhi);
		tempDyPlusPhi  = wenoApproximationYPlus(tempPhi);
		tempDyMinusPhi = wenoApproximationYMinus(tempPhi);

		tempPhi->k1 = -deltaT*propagatingGodonov(tempDxPlusPhi,tempDxMinusPhi,tempDyPlusPhi,tempDyMinusPhi, tempPhi->phi);
		tempPhi = tempPhi->PhiNext;
	}

	tempPhi = PhiHead;
	while (tempPhi!=NULL)
	{
		tempPhi->phi = tempPhi->originPhi + tempPhi->k1;
		tempPhi = tempPhi->PhiNext;
	}

	tempPhi = PhiHead;
	while (tempPhi!=NULL)
	{
		tempDxPlusPhi  = wenoApproximationXPlus(tempPhi);
		tempDxMinusPhi = wenoApproximationXMinus(tempPhi);
		tempDyPlusPhi  = wenoApproximationYPlus(tempPhi);
		tempDyMinusPhi = wenoApproximationYMinus(tempPhi);

		tempPhi->k2 = -deltaT*propagatingGodonov(tempDxPlusPhi,tempDxMinusPhi,tempDyPlusPhi,tempDyMinusPhi, tempPhi->phi);
		tempPhi = tempPhi->PhiNext;
	}

	tempPhi = PhiHead;
	while (tempPhi!=NULL)
	{
		tempPhi->phi = 3.0/4.0*tempPhi->originPhi + 1.0/4.0*tempPhi->phi + 1.0/4.0*tempPhi->k2;
		tempPhi = tempPhi->PhiNext;
	}

	tempPhi = PhiHead;
	while (tempPhi!=NULL)
	{
		tempDxPlusPhi  = wenoApproximationXPlus(tempPhi);
		tempDxMinusPhi = wenoApproximationXMinus(tempPhi);
		tempDyPlusPhi  = wenoApproximationYPlus(tempPhi);
		tempDyMinusPhi = wenoApproximationYMinus(tempPhi);

		tempPhi->k3 = -deltaT*propagatingGodonov(tempDxPlusPhi,tempDxMinusPhi,tempDyPlusPhi,tempDyMinusPhi, tempPhi->phi);
		tempPhi = tempPhi->PhiNext;
	}

	tempPhi = PhiHead;
	while (tempPhi!=NULL)
	{
		tempPhi->phi = 1.0/3.0*tempPhi->originPhi + 2.0/3.0*tempPhi->phi + 2.0/3.0*tempPhi->k3;
		tempPhi = tempPhi->PhiNext;
	}
}