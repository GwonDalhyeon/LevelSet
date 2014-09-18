#include "flux.h"

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



void Shape2D::example3LaxFriedrichs()
{
	Phi2D* tempPhi = PhiHead;
	double tempDxPlusPhi;
	double tempDxMinusPhi;
	double tempDyPlusPhi;
	double tempDyMinusPhi;
	//double velocity;
	while (tempPhi!=NULL)
	{
		tempPhi->originPhi = tempPhi->phi;

		tempDxPlusPhi  = wenoApproximationXPlus(tempPhi);
		tempDxMinusPhi = wenoApproximationXMinus(tempPhi);
		tempDyPlusPhi  = wenoApproximationYPlus(tempPhi);
		tempDyMinusPhi = wenoApproximationYMinus(tempPhi);

		tempPhi->k1 = -deltaT*example3LaxFriedrichsFlux(tempDxPlusPhi,tempDxMinusPhi,tempDyPlusPhi,tempDyMinusPhi);
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

		tempPhi->k2 = -deltaT*example3LaxFriedrichsFlux(tempDxPlusPhi,tempDxMinusPhi,tempDyPlusPhi,tempDyMinusPhi);
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

		tempPhi->k3 = -deltaT*example3LaxFriedrichsFlux(tempDxPlusPhi,tempDxMinusPhi,tempDyPlusPhi,tempDyMinusPhi);
		tempPhi = tempPhi->PhiNext;
	}

	tempPhi = PhiHead;
	while (tempPhi!=NULL)
	{
		tempPhi->phi = 1.0/3.0*tempPhi->originPhi + 2.0/3.0*tempPhi->phi + 2.0/3.0*tempPhi->k3;
		tempPhi = tempPhi->PhiNext;
	}
}


void Shape2D::singleVortex()
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

		tempPhi->k1 = -tempPhi->xVelocity*deltaT*minmod(tempDxPlusPhi,tempDxMinusPhi) -tempPhi->yVelocity*deltaT*minmod(tempDyPlusPhi,tempDyMinusPhi);
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

		tempPhi->k2 = -tempPhi->xVelocity*deltaT*minmod(tempDxPlusPhi,tempDxMinusPhi) -tempPhi->yVelocity*deltaT*minmod(tempDyPlusPhi,tempDyMinusPhi);
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

		tempPhi->k3 = -tempPhi->xVelocity*deltaT*minmod(tempDxPlusPhi,tempDxMinusPhi) -tempPhi->yVelocity*deltaT*minmod(tempDyPlusPhi,tempDyMinusPhi);
		tempPhi = tempPhi->PhiNext;
	}

	tempPhi = PhiHead;
	while (tempPhi!=NULL)
	{
		tempPhi->phi = 1.0/3.0*tempPhi->originPhi + 2.0/3.0*tempPhi->phi + 2.0/3.0*tempPhi->k3;
		tempPhi = tempPhi->PhiNext;
	}
}




void Shape2D::singleVortexParticle()
{
	Particle2D* tempParticle = ParticleHead;
	double xCoord, yCoord;
	//double originPhiValue;

	while (tempParticle!=NULL)
	{
		xCoord = tempParticle->x0;
		yCoord = tempParticle->y0;

		if (xCoord<X0 || xCoord >X1 || yCoord<Y0 || yCoord >Y1)
		{
			tempParticle = tempParticle->ParticleNext;
			deleteParticle(tempParticle->ParticleBefore);
			continue;
		}
		//findCellContainingParticle(tempParticle, CellHead);
		//interpolationParticleVelocity(tempParticle);

		tempParticle->phi1 = tempParticle->phi0;
		//tempParticle->phi1 = originPhiValue;

		tempParticle->kx1 = tempParticle->xVelocity*deltaT;
		tempParticle->ky1 = tempParticle->yVelocity*deltaT;

		tempParticle->x1 = tempParticle->x0 + tempParticle->kx1;
		tempParticle->y1 = tempParticle->y0 + tempParticle->ky1;

		tempParticle->x0 = tempParticle->x1;
		tempParticle->y0 = tempParticle->y1;

		findCellContainingParticle(tempParticle, CellHead);
		interpolationParticleVelocity(tempParticle);

		tempParticle->kx2 = tempParticle->xVelocity*deltaT;
		tempParticle->ky2 = tempParticle->yVelocity*deltaT;

		tempParticle->x1 = 3.0/4.0*xCoord + 1.0/4.0*tempParticle->x1 + 1.0/4.0*tempParticle->kx2;
		tempParticle->y1 = 3.0/4.0*yCoord + 1.0/4.0*tempParticle->y1 + 1.0/4.0*tempParticle->ky2;

		tempParticle->x0 = tempParticle->x1;
		tempParticle->y0 = tempParticle->y1;

		findCellContainingParticle(tempParticle, CellHead);
		interpolationParticleVelocity(tempParticle);

		tempParticle->kx3 = tempParticle->xVelocity*deltaT;
		tempParticle->ky3 = tempParticle->yVelocity*deltaT;

		tempParticle->x0 = 1.0/3.0*xCoord + 2.0/3.0*tempParticle->x1 + 2.0/3.0*tempParticle->kx3;
		tempParticle->y0 = 1.0/3.0*yCoord + 2.0/3.0*tempParticle->y1 + 2.0/3.0*tempParticle->ky3;

		findCellContainingParticle(tempParticle, CellHead);
		interpolationParticleVelocity(tempParticle);


		interpolationParticlePhi(tempParticle);


		//if ((abs(tempParticle->phi0) > tempParticle->radius) && (tempParticle->phi0 * tempParticle->phi1 < 0) )
		if ((tempParticle->phi0 * tempParticle->phi1 < 0) )
		{
			tempParticle->escapedFlag = true;
		}
		else
		{
			tempParticle->escapedFlag = false;
		}

		tempParticle = tempParticle->ParticleNext;
	}
}



void Shape2D::reductionError()
{
	Phi2D* tempPhi = PhiHead;
	Particle2D* tempParticle;

	double tempPhiPlus, tempPhiMinus;
	double phiP1, phiP2;

	double tempPhiX0, tempPhiX1, tempPhiY0, tempPhiY1, tempParticleX, tempParticleY;

	while (tempPhi!=NULL)
	{
		//if (abs(tempPhi->phi)<bandMax)
		//{

		tempPhiPlus  = tempPhi->phi;
		tempPhiMinus = tempPhi->phi;
		tempParticle = ParticleHead;

		tempPhiX0 = tempPhi->x;
		tempPhiY0 = tempPhi->y;

		while (tempParticle!=NULL)
		{
			if (tempParticle->escapedFlag)
			{
				if (tempParticle->phi1>0)
				{
					phiP1 = sign2(tempParticle->phi1)*(tempParticle->radius - sqrt( (tempParticle->x0 - tempPhiX0)*(tempParticle->x0 - tempPhiX0) + (tempParticle->y0 - tempPhiY0)*(tempParticle->y0 - tempPhiY0) ));
					tempPhiPlus = max(tempPhiPlus, phiP1);
				}
				else
				{
					phiP2 = sign2(tempParticle->phi1)*(tempParticle->radius - sqrt( (tempParticle->x0 - tempPhiX0)*(tempParticle->x0 - tempPhiX0) + (tempParticle->y0 - tempPhiY0)*(tempParticle->y0 - tempPhiY0) ));
					tempPhiMinus = min(tempPhiMinus, phiP2);
				}
			}

			//particleRadius(tempParticle);
			tempParticle = tempParticle->ParticleNext;
		}

		if (abs(tempPhiPlus)<=abs(tempPhiMinus))
		{
			tempPhi->phi = tempPhiPlus;
		}
		else
		{
			tempPhi->phi = tempPhiMinus;
		}
		
			//tempPhi->phi = (tempPhiPlus + tempPhiMinus)/2.0;
		

		//}


		tempPhi = tempPhi->PhiNext;
	}
}





