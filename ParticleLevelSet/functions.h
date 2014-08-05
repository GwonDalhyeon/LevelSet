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
	double velocity;
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

	while (tempParticle!=NULL)
	{
		interpolationParticleVelocity(tempParticle);

		tempParticle->kx1 = -tempParticle->xVelocity*deltaT;
		tempParticle->ky1 = -tempParticle->yVelocity*deltaT;

		tempParticle->x1 = tempParticle->x0 + tempParticle->kx1;
		tempParticle->y1 = tempParticle->y0 + tempParticle->ky1;


		tempParticle->kx2 = -tempParticle->xVelocity*deltaT;
		tempParticle->ky2 = -tempParticle->yVelocity*deltaT;

		tempParticle->x1 = tempParticle->x1 + tempParticle->kx2;
		tempParticle->y1 = tempParticle->y1 + tempParticle->ky2;

		tempParticle->kx3 = -tempParticle->xVelocity*deltaT;
		tempParticle->ky3 = -tempParticle->yVelocity*deltaT;

		tempParticle->x0 = 1.0/3.0*tempParticle->x0 + 2.0/3.0*tempParticle->x1 + 2.0/3.0*tempParticle->kx3;
		tempParticle->y0 = 1.0/3.0*tempParticle->y0 + 2.0/3.0*tempParticle->y1 + 2.0/3.0*tempParticle->ky3;

		tempParticle = tempParticle->ParticleNext;
	}
}


