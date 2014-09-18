#include "phi.h"
#include <vector>

struct Cell2D;
struct Particle2D;

struct Cell2D
{
	Phi2D* PhiLeftBottom;
	Phi2D* PhiLeftTop;
	Phi2D* PhiRightBottom;
	Phi2D* PhiRightTop;

	double x0,x1,y0,y1;
	double cellLength;
	double CellLevel;

	Cell2D* CellBefore;
	Cell2D* CellNext;
	Cell2D* CellParent;
	Cell2D* CellChildLeftBottom;
	Cell2D* CellChildRightBottom;
	Cell2D* CellChildLeftTop;
	Cell2D* CellChildRightTop;

	Particle2D* containParticle;
	vector<Particle2D*> particleVector;
	vector<double> containParticleValue;
	

	int numContainParticle;
	int particlePlacedFlag;

	Cell2D();
};

Cell2D::Cell2D()
{
	PhiLeftBottom  = NULL;
	PhiLeftTop     = NULL;
	PhiRightBottom = NULL;
	PhiRightTop    = NULL;



	CellBefore = NULL;
	CellNext   = NULL;
	CellParent = NULL;
	CellChildLeftBottom  = NULL;
	CellChildRightBottom = NULL;
	CellChildLeftTop     = NULL;
	CellChildRightTop    = NULL;

	//containParticle = NULL;
	
	//containParticleValue = NULL;

	numContainParticle = 0;
	particlePlacedFlag = 0;
}


void deleteAllCell(Cell2D* inputCell)
{
	if (inputCell != NULL)
	{
		deleteAllCell(inputCell->CellNext);
		delete inputCell;
	}
}





struct Particle2D
{
	double x0, y0;
	double x1, y1;
	double radius;

	double normalX, normalY;

	double phi0;
	double phi1;
	double tempPhiValue;
	double originPhi;

	bool escapedFlag;

	double xVelocity;
	double yVelocity;
	//double velocityField;

	double kx1,kx2,kx3;
	double ky1,ky2,ky3;

	int xIndex, yIndex;
	//struct Particle2D* PhiLeft;
	//struct Particle2D* PhiRight;
	//struct Particle2D* PhiTop;
	//struct Particle2D* PhiBottom;

	Cell2D* containedCell;

	Particle2D* ParticleNext;
	Particle2D* ParticleBefore;
	Particle2D();
};

Particle2D::Particle2D()
{
	//PhiLeft   = NULL;
	//PhiRight  = NULL;
	//PhiTop    = NULL;
	//PhiBottom = NULL;

	ParticleNext   = NULL;
	ParticleBefore = NULL;
}


void deleteAllParticle(Particle2D* inputParticle)
{
	if (inputParticle != NULL)
	{
		deleteAllParticle(inputParticle ->ParticleNext);
		delete inputParticle;
	}
}