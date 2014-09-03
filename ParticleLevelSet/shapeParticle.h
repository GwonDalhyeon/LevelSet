#include "shapeCell.h"



void Shape2D::savingParticle(int timeIndex)
{
	Particle2D* tempParticle=ParticleHead;

	ofstream particleData;
	particleData.open("D:\Data/particle"+to_string(timeIndex)+".txt");

	while (tempParticle!=NULL)
	{
		particleData<<tempParticle->x0<<" "<<tempParticle->y0<< " " << tempParticle->phi0;
		particleData<<endl;
		particleRadius(tempParticle);
		tempParticle = tempParticle->ParticleNext;
	}
	particleData.close();



	//tempParticle=ParticleHead;

	//ofstream particleData;
	//particleData.open("D:\Data/negativeParticle"+to_string(timeIndex)+".txt");

	//while (tempParticle!=NULL && tempParticle->phi0<0)
	//{
	//	particleData<<tempParticle->x0<<" "<<tempParticle->y0;
	//	particleData<<endl;

	//	tempParticle = tempParticle->ParticleNext;
	//}
	//particleData.close();
}


void Shape2D::deleteParticle(Particle2D* inputParticle)
{
	if (inputParticle->ParticleBefore == NULL)
	{
		ParticleHead = inputParticle->ParticleNext;
		inputParticle->ParticleNext->ParticleBefore = NULL;
		delete inputParticle;
	}
	else if (inputParticle->ParticleNext == NULL)
	{
		ParticleTail = inputParticle->ParticleBefore;
		inputParticle->ParticleBefore->ParticleNext = NULL;
		delete inputParticle;
	}
	else if (inputParticle->ParticleBefore != NULL && inputParticle->ParticleNext != NULL)
	{
		Particle2D* tempParticle = inputParticle;
		inputParticle->ParticleBefore->ParticleNext = inputParticle->ParticleNext;
		inputParticle->ParticleNext->ParticleBefore = inputParticle->ParticleBefore;
		delete inputParticle;
	}
}


void Shape2D::findCellContainingParticle(Particle2D* inputParticle, Cell2D* inputCell)
{

	if (inputCell->x0<=inputParticle->x0 && inputCell->x1>=inputParticle->x0 && inputCell->y0<=inputParticle->y0 && inputCell->y1>=inputParticle->y0)
	{
		if ((int)inputCell->CellLevel == cellLevel)
		{
			inputParticle->containedCell = inputCell;
			return;
		}
		else
		{
			findCellContainingParticle(inputParticle, inputCell->CellChildLeftBottom);
			findCellContainingParticle(inputParticle, inputCell->CellChildLeftTop);
			findCellContainingParticle(inputParticle, inputCell->CellChildRightBottom);
			findCellContainingParticle(inputParticle, inputCell->CellChildRightTop);
		}
	}


}


void Shape2D::interpolationParticlePhi(Particle2D* inputParticle)
{
	double xCoord = inputParticle->x0;
	double yCoord = inputParticle->y0;


	Phi2D* tempLeftBottomPhi = inputParticle->containedCell->PhiLeftBottom;
	Phi2D* tempRightBottomPhi= inputParticle->containedCell->PhiRightBottom;
	Phi2D* tempLeftTopPhi    = inputParticle->containedCell->PhiLeftTop;
	Phi2D* tempRightTopPhi   = inputParticle->containedCell->PhiRightTop;


	double tempX0 = tempLeftBottomPhi->x;
	double tempX1 = tempRightTopPhi->x;
	double tempY0 = tempLeftBottomPhi->y;
	double tempY1 = tempRightTopPhi->y;

	double phixx1 = min(abs(dxxPhi(tempLeftBottomPhi)),abs(dxxPhi(tempRightBottomPhi)));
	double phixx2 = min(abs(dxxPhi(tempLeftTopPhi)),abs(dxxPhi(tempRightTopPhi)) );
	double phixx  = min(phixx1, phixx2);

	double phiyy1 = min(abs(dyyPhi(tempLeftBottomPhi)),abs(dyyPhi(tempRightBottomPhi)));
	double phiyy2 = min(abs(dyyPhi(tempLeftTopPhi)),abs(dyyPhi(tempRightTopPhi)) );
	double phiyy  = min(phiyy1, phiyy2);

	double a1 = tempLeftBottomPhi->phi*(tempX1-xCoord)/(tempX1-tempX0)*(tempY1-yCoord)/(tempY1-tempY0);
	double a2 = tempLeftTopPhi->phi*(tempX1-xCoord)/(tempX1-tempX0)*(yCoord-tempY0)/(tempY1-tempY0);
	double a3 = tempRightBottomPhi->phi*(xCoord-tempX0)/(tempX1-tempX0)*(tempY1-yCoord)/(tempY1-tempY0);
	double a4 = tempRightTopPhi->phi*(xCoord-tempX0)/(tempX1-tempX0)*(yCoord-tempY0)/(tempY1-tempY0);

	inputParticle->phi0 = a1 + a2 + a3 + a4 - phixx*(tempX1-xCoord)*(xCoord-tempX0)/2.0 - phiyy*(tempY1-yCoord)*(yCoord-tempY0)/2.0 ;
}



void Shape2D::interpolationParticleVelocity(Particle2D* inputParticle)
{
	double xCoord = inputParticle->x0;
	double yCoord = inputParticle->y0;
	
	Phi2D* tempLeftBottomPhi  = inputParticle->containedCell->PhiLeftBottom;
	Phi2D* tempRightBottomPhi = inputParticle->containedCell->PhiRightBottom;
	Phi2D* tempLeftTopPhi     = inputParticle->containedCell->PhiLeftTop;
	Phi2D* tempRightTopPhi    = inputParticle->containedCell->PhiRightTop;

	double tempX0 = tempLeftBottomPhi->x;
	double tempX1 = tempRightTopPhi->x;
	double tempY0 = tempLeftBottomPhi->y;
	double tempY1 = tempRightTopPhi->y;

	double a1 = tempLeftBottomPhi->xVelocity*(tempX1-xCoord)/(tempX1-tempX0)*(tempY1-yCoord)/(tempY1-tempY0);
	double a2 = tempLeftTopPhi->xVelocity*(tempX1-xCoord)/(tempX1-tempX0)*(yCoord-tempY0)/(tempY1-tempY0);
	double a3 = tempRightBottomPhi->xVelocity*(xCoord-tempX0)/(tempX1-tempX0)*(tempY1-yCoord)/(tempY1-tempY0);
	double a4 = tempRightTopPhi->xVelocity*(xCoord-tempX0)/(tempX1-tempX0)*(yCoord-tempY0)/(tempY1-tempY0);

	inputParticle->xVelocity = a1 + a2 + a3 + a4;

	double b1 = tempLeftBottomPhi->yVelocity*(tempX1-xCoord)/(tempX1-tempX0)*(tempY1-yCoord)/(tempY1-tempY0);
	double b2 = tempLeftTopPhi->yVelocity*(tempX1-xCoord)/(tempX1-tempX0)*(yCoord-tempY0)/(tempY1-tempY0);
	double b3 = tempRightBottomPhi->yVelocity*(xCoord-tempX0)/(tempX1-tempX0)*(tempY1-yCoord)/(tempY1-tempY0);
	double b4 = tempRightTopPhi->yVelocity*(xCoord-tempX0)/(tempX1-tempX0)*(yCoord-tempY0)/(tempY1-tempY0);

	inputParticle->yVelocity = b1 + b2 + b3 + b4;

	//inputParticle->xVelocity =  2.0*sin(PI*(xCoord))*sin(PI*(xCoord))*sin(PI*(yCoord))*cos(PI*(yCoord));
	//inputParticle->yVelocity = -2.0*sin(PI*(xCoord))*cos(PI*(xCoord))*sin(PI*(yCoord))*sin(PI*(yCoord));

	double tempVel1 = 2.0*sin(PI*(xCoord))*sin(PI*(xCoord))*sin(PI*(yCoord))*cos(PI*(yCoord));
	double tempVel2 = -2.0*sin(PI*(xCoord))*cos(PI*(xCoord))*sin(PI*(yCoord))*sin(PI*(yCoord));

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
			ParticleHead->containedCell = inputCell;

			interpolationParticleVelocity(ParticleHead);
			interpolationParticlePhi(ParticleHead);
			attractingParticle(ParticleHead);
		}
		else
		{

			Particle2D* newParticle = new Particle2D();

			newParticle->x0 = randX;
			newParticle->y0 = randY;

			newParticle->containedCell = inputCell;

			interpolationParticleVelocity(newParticle);
			interpolationParticlePhi(newParticle);
			attractingParticle(newParticle);

			if (newParticle !=NULL)
			{
				Particle2D* oldParticleTail = ParticleTail;
				oldParticleTail->ParticleNext = newParticle;
				newParticle->ParticleBefore = oldParticleTail;
				ParticleTail = newParticle;
			}

		}

		particleNum = particleNum + 1;
	}

}


void Shape2D::phiNormalVector(Phi2D* inputPhi, double& phix, double& phiy)
{
	//double phix, phiy;
	phix = dxPhi(inputPhi);
	phiy = dyPhi(inputPhi);

	phix = phix/sqrt(phix*phix + phiy*phiy);
	phiy = phiy/sqrt(phix*phix + phiy*phiy);
}

void Shape2D::particleNormalVector(Particle2D* inputParticle)
{
	double phix1, phix2, phix3, phix4;
	double phiy1, phiy2, phiy3, phiy4;
	double xCoord = inputParticle->x0;
	double yCoord = inputParticle->y0;

	double tempX0 = inputParticle->containedCell->PhiLeftBottom->x;
	double tempX1 = inputParticle->containedCell->PhiRightTop->x;
	double tempY0 = inputParticle->containedCell->PhiLeftBottom->y;
	double tempY1 = inputParticle->containedCell->PhiRightTop->y;

	phiNormalVector(inputParticle->containedCell->PhiLeftBottom, phix1, phiy1);
	phiNormalVector(inputParticle->containedCell->PhiLeftTop, phix2, phiy2);
	phiNormalVector(inputParticle->containedCell->PhiRightBottom, phix3, phiy3);
	phiNormalVector(inputParticle->containedCell->PhiRightTop, phix4, phiy4);

	double a1 = phix1*(tempX1-xCoord)/(tempX1-tempX0)*(tempY1-yCoord)/(tempY1-tempY0);
	double a2 = phix2*(tempX1-xCoord)/(tempX1-tempX0)*(yCoord-tempY0)/(tempY1-tempY0);
	double a3 = phix3*(xCoord-tempX0)/(tempX1-tempX0)*(tempY1-yCoord)/(tempY1-tempY0);
	double a4 = phix4*(xCoord-tempX0)/(tempX1-tempX0)*(yCoord-tempY0)/(tempY1-tempY0);

	inputParticle->normalX = a1 + a2 + a3 + a4;

	double b1 = phiy1*(tempX1-xCoord)/(tempX1-tempX0)*(tempY1-yCoord)/(tempY1-tempY0);
	double b2 = phiy2*(tempX1-xCoord)/(tempX1-tempX0)*(yCoord-tempY0)/(tempY1-tempY0);
	double b3 = phiy3*(xCoord-tempX0)/(tempX1-tempX0)*(tempY1-yCoord)/(tempY1-tempY0);
	double b4 = phiy4*(xCoord-tempX0)/(tempX1-tempX0)*(yCoord-tempY0)/(tempY1-tempY0);

	inputParticle->normalY = b1 + b2 + b3 + b4;
}

void Shape2D::attractingParticle(Particle2D* inputParticle)
{
	double goalPhi = (double) sign2(inputParticle->phi0)* rand()/RAND_MAX*(bandMax-bandMin)+bandMin;
	int iterationNum = 0;
	double lambda = 1.0;
	//double phix, phiy;

	while (iterationNum < 15)
	{
		interpolationParticlePhi(inputParticle);
		particleRadius(inputParticle);


		if (abs(inputParticle->phi0)>bandMin && abs(inputParticle->phi0)<bandMax )
		{
			break;
		}

		particleNormalVector(inputParticle);

		inputParticle->x0 = inputParticle->x0 + lambda*(goalPhi-inputParticle->phi0)*inputParticle->normalX;
		inputParticle->y0 = inputParticle->y0 + lambda*(goalPhi-inputParticle->phi0)*inputParticle->normalY;

		//interpolationParticlePhi(inputParticle);

		lambda = lambda/2.0;
		iterationNum += 1;
	}

	//deleteParticle(inputParticle);
	return;
}

void Shape2D::particleRadius(Particle2D* inputParticle)
{
	if (abs(inputParticle->phi0)>radiusMax)
	{
		inputParticle->radius = radiusMax;
	}
	else if (abs(inputParticle->phi0)<radiusMin)
	{
		inputParticle->radius = radiusMin;
	}
	else
	{
		inputParticle->radius = abs(inputParticle->phi0);
	}
}