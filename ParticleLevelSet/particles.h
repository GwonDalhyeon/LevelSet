#include "phi.h"

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