#include "shapeCell.h"



void Shape2D::savingParticle(int timeIndex)
{
	Particle2D* tempParticle=ParticleHead;

	ofstream particleData;
	particleData.open("D:\Data/particle"+to_string(timeIndex)+".txt");

	while (tempParticle!=NULL)
	{
		particleData<<tempParticle->x0<<" "<<tempParticle->y0;
		particleData<<endl;

		tempParticle = tempParticle->ParticleNext;
	}
	particleData.close();
}


void Shape2D::deleteParticle(Particle2D* inputParticle)
{
	if (inputParticle->ParticleBefore == NULL && inputParticle->ParticleNext != NULL)
	{
		ParticleHead = inputParticle->ParticleNext;
		delete inputParticle;
	}
	else if (inputParticle->ParticleBefore != NULL && inputParticle->ParticleNext == NULL)
	{
		ParticleTail = inputParticle->ParticleBefore;
		delete inputParticle;
	}
	else if (inputParticle->ParticleBefore != NULL && inputParticle->ParticleNext != NULL)
	{
		inputParticle->ParticleBefore->ParticleNext = inputParticle->ParticleNext;
		inputParticle->ParticleNext->ParticleBefore = inputParticle->ParticleBefore;
		delete inputParticle;
	}
}


void Shape2D::findCellContainingParticle(Particle2D* inputParticle)
{
	Cell2D* tempCell = CellHead;
	while (tempCell != NULL)
	{
		if ((int)tempCell->CellLevel == cellLevel)
		{
			if (tempCell->x0<=inputParticle->x0 && tempCell->x1>=inputParticle->x0 && tempCell->y0<=inputParticle->x0 && tempCell->y1>=inputParticle->y0)
			{
				inputParticle->containdCell = tempCell;
				return;
			}
		}
		tempCell = tempCell->CellNext;
	}

}


void Shape2D::interpolationParticlePhi(Particle2D* inputParticle)
{
	double xCoord = inputParticle->x0;
	double yCoord = inputParticle->y0;
	if (xCoord<X0 || xCoord >X1 || yCoord<Y0 || yCoord >Y1)
	{
		deleteParticle(inputParticle);
		return;
	}

	findCellContainingParticle(inputParticle);

	Phi2D* tempLeftBottomPhi = inputParticle->containdCell->PhiLeftBottom;
	Phi2D* tempRightBottomPhi= inputParticle->containdCell->PhiRightBottom;
	Phi2D* tempLeftTopPhi    = inputParticle->containdCell->PhiLeftTop;
	Phi2D* tempRightTopPhi   = inputParticle->containdCell->PhiRightTop;
	 

	double tempX0 = tempLeftBottomPhi->x;
	double tempX1 = tempRightTopPhi->x;
	double tempY0 = tempLeftBottomPhi->y;
	double tempY1 = tempRightTopPhi->y;
	

	double phixx = dxxPhi(tempLeftBottomPhi);
	double phiyy = dyyPhi(tempLeftBottomPhi);

	if (abs(phixx) > abs(dxxPhi(tempRightBottomPhi)))
	{
		phixx = dxxPhi(tempRightBottomPhi);
	}
	if (abs(phixx) > abs(dxxPhi(tempLeftTopPhi)))
	{
		phixx = dxxPhi(tempLeftTopPhi);
	}
	if (abs(phixx) > abs(dxxPhi(tempRightTopPhi)))
	{
		phixx = dxxPhi(tempRightTopPhi);
	}

	if (abs(phiyy) > abs(dyyPhi(tempRightBottomPhi)))
	{
		phiyy = dyyPhi(tempRightBottomPhi);
	}
	if (abs(phiyy) > abs(dyyPhi(tempLeftTopPhi)))
	{
		phiyy = dyyPhi(tempLeftTopPhi);
	}
	if (abs(phiyy) > abs(dyyPhi(tempRightTopPhi)))
	{
		phiyy = dyyPhi(tempRightTopPhi);
	}

	double a1 = tempLeftBottomPhi->phi*(tempX0-xCoord)/(tempX1-tempX0)*(yCoord-tempY0)/(tempY1-tempY0);
	double a2 = tempLeftTopPhi->phi*(tempX0-xCoord)/(tempX1-tempX0)*(tempY1-yCoord)/(tempY1-tempY0);
	double a3 = tempRightBottomPhi->phi*(tempX1-xCoord)/(tempX1-tempX0)*(yCoord-tempY0)/(tempY1-tempY0);
	double a4 = tempRightTopPhi->phi*(tempX1-xCoord)/(tempX1-tempX0)*(tempY1-yCoord)/(tempY1-tempY0);

	inputParticle->phi0 = a1 + a2 + a3 + a4 - phixx*(tempX1-xCoord)*(xCoord-tempX0)/2.0 - phiyy*(tempY1-yCoord)*(yCoord-tempY0)/2.0 ;
}



void Shape2D::interpolationParticleVelocity(Particle2D* inputParticle)
{
	double xCoord = inputParticle->x0;
	double yCoord = inputParticle->y0;
	//if (xCoord<X0 || xCoord >X1 || yCoord<Y0 || yCoord >Y1)
	//{
	//	deleteParticle(inputParticle);
	//	return;
	//}


	findCellContainingParticle(inputParticle);
	
	Phi2D* tempLeftBottomPhi  = inputParticle->containdCell->PhiLeftBottom;
	Phi2D* tempRightBottomPhi = inputParticle->containdCell->PhiRightBottom;
	Phi2D* tempLeftTopPhi     = inputParticle->containdCell->PhiLeftTop;
	Phi2D* tempRightTopPhi    = inputParticle->containdCell->PhiRightTop;

	double tempX0 = tempLeftBottomPhi->x;
	double tempX1 = tempRightTopPhi->x;
	double tempY0 = tempLeftBottomPhi->y;
	double tempY1 = tempRightTopPhi->y;

	double a1 = tempLeftBottomPhi->xVelocity*(xCoord-tempX0)/(tempX1-tempX0)*(yCoord-tempY0)/(tempY1-tempY0);
	double a2 = tempLeftTopPhi->xVelocity*(xCoord-tempX0)/(tempX1-tempX0)*(tempY1-yCoord)/(tempY1-tempY0);
	double a3 = tempRightBottomPhi->xVelocity*(tempX1-xCoord)/(tempX1-tempX0)*(yCoord-tempY0)/(tempY1-tempY0);
	double a4 = tempRightTopPhi->xVelocity*(tempX1-xCoord)/(tempX1-tempX0)*(tempY1-yCoord)/(tempY1-tempY0);

	inputParticle->xVelocity = a1 + a2 + a3 + a4;

	double b1 = tempLeftBottomPhi->yVelocity*(xCoord-tempX0)/(tempX1-tempX0)*(yCoord-tempY0)/(tempY1-tempY0);
	double b2 = tempLeftTopPhi->yVelocity*(xCoord-tempX0)/(tempX1-tempX0)*(tempY1-yCoord)/(tempY1-tempY0);
	double b3 = tempRightBottomPhi->yVelocity*(tempX1-xCoord)/(tempX1-tempX0)*(yCoord-tempY0)/(tempY1-tempY0);
	double b4 = tempRightTopPhi->yVelocity*(tempX1-xCoord)/(tempX1-tempX0)*(tempY1-yCoord)/(tempY1-tempY0);

	inputParticle->yVelocity = b1 + b2 + b3 + b4;

	//inputParticle->xVelocity =  2.0*sin(PI*(xCoord))*sin(PI*(xCoord))*sin(PI*(yCoord))*cos(PI*(yCoord));
	//inputParticle->yVelocity = -2.0*sin(PI*(xCoord))*cos(PI*(xCoord))*sin(PI*(yCoord))*sin(PI*(yCoord));
}

void Shape2D::initializationParticle()
{
	double particleBandWidth = 3.0*max(deltaX, deltaY);
	int particlePerDim = 4;

	Cell2D* tempCell = CellHead;
	while (tempCell != NULL)
	{
		if ((int)tempCell->CellLevel == cellLevel)
		{
			if (abs(tempCell->PhiLeftBottom->phi)<particleBandWidth || abs(tempCell->PhiLeftTop->phi)<particleBandWidth || abs(tempCell->PhiRightBottom->phi)<particleBandWidth || abs(tempCell->PhiRightTop->phi)<particleBandWidth)
			{
				tempCell->particlePlacedFlag = 1;
				sprinkleParticle(tempCell);
			}
		}
		tempCell = tempCell->CellNext;
	}
}

void Shape2D::sprinkleParticle(Cell2D* inputCell)
{
	double tempX0 = inputCell->x0;
	double tempX1 = inputCell->x1;
	double tempY0 = inputCell->y0;
	double tempY1 = inputCell->y1;

	double randX;
	double randY;
	
	int particleNum = 0;
	
	

	while (particleNum < 16)
	{
		randX = (double) rand()/RAND_MAX*(tempX1- tempX0)+tempX0;
		randY = (double) rand()/RAND_MAX*(tempY1- tempY0)+tempY0;
		
		if (ParticleHead == NULL)
		{
			ParticleHead = new Particle2D();
			ParticleHead->x0 = randX;
			ParticleHead->y0 = randY;
			ParticleHead->ParticleNext = ParticleTail;
			ParticleTail = ParticleHead;
			ParticleHead->containdCell = inputCell;
			interpolationParticleVelocity(ParticleHead);
			interpolationParticlePhi(ParticleHead);
		}
		else
		{
			
			Particle2D* newParticle = new Particle2D();

			newParticle->x0 = randX;
			newParticle->y0 = randY;

			newParticle->containdCell = inputCell;

			interpolationParticleVelocity(newParticle);
			interpolationParticlePhi(newParticle);

			Particle2D* oldParticleTail = ParticleTail;
			oldParticleTail->ParticleNext = newParticle;
			newParticle->ParticleBefore = oldParticleTail;
			ParticleTail = newParticle;
		}

		particleNum = particleNum + 1;
	}

}




