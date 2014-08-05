#include "operator.h"


class Shape2D
{
public:

	double deltaT;
	double deltaX;
	double deltaY;
	double T, X0, X1, Y0, Y1;
	int xPointNum;
	int yPointNum;
	int numPoint;

	Phi2D* PhiHead;
	Phi2D* PhiTail;

	Cell2D* CellHead;
	Cell2D* CellTail;
	double cellLevel;

	Particle2D* ParticleHead;
	Particle2D* ParticleTail;



	double** savedPhi;
	double** cellParticle;

	//Phi based contructor.
	Shape2D(double inputX0, double inputX1, int inputXPointNum, double inputY0, double inputY1, int inputYPointNum, double inputT, double inputdeltaT);
	Shape2D(double inputX0, double inputX1, double inputDeltaX, double inputY0, double inputY1, double inputDeltaY, double inputT, double inputdeltaT);
	
	// Cell based constructor.
	Shape2D(double inputX0, double inputX1, double inputY0, double inputY1, double inputLevel, double inputT, double inputdeltaT);
	~Shape2D();

	// About Phi
	Phi2D* addPhi(double inputX, double inputY, Phi2D* leftPhi, Phi2D* rightPhi, Phi2D* lowerPhi, Phi2D* UpperPhi);
	void connectingPhi(Phi2D* inputPhi, int inputXIndex, int inputYIndex);
	void connectingPhi();
	Phi2D* findPoint(Phi2D* phiHead, int inputXIndex, int inputYIndex);
	bool findPointBool(Phi2D* phiHead, int inputXIndex, int inputYIndex);
	void savingPhi(int timeIndex);
	void indexingPhi(Phi2D* inputPhi);

	// About Cell
	Cell2D* addCell(Phi2D* leftBottomPhi, Phi2D* rightBottomPhi, Phi2D* leftTopPhi, Phi2D* RightTopPhi, double inputLevel);
	void combiningCellPhi(Cell2D* inputCell, double inputLevel);


	// About Particle
	void initializationParticle();
	void savingParticle(int timeIndex);
	void interpolationParticlePhi(Particle2D* inputParticle);
	void interpolationParticleVelocity(Particle2D* inputParticle);
	void deleteParticle(Particle2D* inputParticle);
	void sprinkleParticle(Cell2D* inputCell);
	void findCellContainingParticle(Particle2D* inputParicle);

	void reinitialTVDRK3();
	void propagatingTVDRK3();


	double example3LaxFriedrichsFlux(double uPlus, double uMinus, double vPlus, double vMinus);
	void   example3LaxFriedrichs();

	void zalesakDisk();
	void singleVortex();
	void singleVortexParticle();
private:

};

Shape2D::Shape2D(double inputX0, double inputX1, int inputXPointNum, double inputY0, double inputY1, int inputYPointNum, double inputT, double inputdeltaT)
{
	deltaT    = inputdeltaT;
	deltaX    = (inputX1 - inputX0)/(double)inputXPointNum;
	deltaY    = (inputY1 - inputY0)/(double)inputYPointNum;

	T  = inputT;
	X0 =  inputX0;
	X1 = inputX1;
	Y0 = inputY0;
	Y1 = inputY1;

	xPointNum = inputXPointNum+1;
	yPointNum = inputYPointNum+1;
	numPoint = 0;
	
	savedPhi = new double* [xPointNum];
	for (int i = 0; i < xPointNum; i++)
	{
		savedPhi[i] = new double[yPointNum];
	}

	cellParticle = new double* [xPointNum-1];
	for (int i = 0; i < xPointNum-1; i++)
	{
		cellParticle[i] = new double[yPointNum-1];
	}

	for (int i = 0; i < xPointNum-1; i++)
	{
		for (int j = 0; j < yPointNum-1; j++)
		{
			cellParticle[i][j] = 0;
		}
	}

	PhiHead    = new Phi2D();
	PhiHead->x = inputX0;
	PhiHead->y = inputY0;
	PhiHead->xIndex = 0;
	PhiHead->yIndex = 0;
	PhiTail    = NULL;

	CellHead = NULL;
	CellTail = NULL;

	ParticleHead = NULL;
	ParticleTail = NULL;
}


Shape2D::Shape2D(double inputX0, double inputX1, double inputDeltaX, double inputY0, double inputY1, double inputDeltaY, double inputT, double inputdeltaT)
{
	deltaT    = inputdeltaT;
	deltaX    = inputDeltaX;
	deltaY    = inputDeltaY;

	T  = inputT;
	X0 = inputX0;
	X1 = inputX1;
	Y0 = inputY0;
	Y1 = inputY1;

	xPointNum = (int)floor((inputX1 - inputX0)/inputDeltaX)+1;
	yPointNum = (int)floor((inputY1 - inputY0)/inputDeltaY)+1;
	numPoint = 0;

	savedPhi = new double* [xPointNum];
	for (int i = 0; i < xPointNum; i++)
	{
		savedPhi[i] = new double[yPointNum];
	}

	cellParticle = new double* [xPointNum-1];
	for (int i = 0; i < xPointNum-1; i++)
	{
		cellParticle[i] = new double[yPointNum-1];
	}

	for (int i = 0; i < xPointNum-1; i++)
	{
		for (int j = 0; j < yPointNum-1; j++)
		{
			cellParticle[i][j] = 0;
		}
	}



	PhiHead    = new Phi2D();
	PhiHead->x = inputX0;
	PhiHead->y = inputY0;
	PhiHead->xIndex = 0;
	PhiHead->yIndex = 0;
	PhiTail    = NULL;

	CellHead = NULL;
	CellTail = NULL;

	ParticleHead = NULL;
	ParticleTail = NULL;
}

Shape2D::Shape2D(double inputX0, double inputX1, double inputY0, double inputY1, double inputLevel, double inputT, double inputdeltaT)
{
	deltaT    = inputdeltaT;
	deltaX    = (inputX1-inputX0)/pow(2.0,inputLevel);
	deltaY    = (inputY1-inputY0)/pow(2.0,inputLevel);

	T  = inputT;
	X0 = inputX0;
	X1 = inputX1;
	Y0 = inputY0;
	Y1 = inputY1;

	xPointNum = (int)pow(2.0,inputLevel)+1;
	yPointNum = (int)pow(2.0,inputLevel)+1;
	numPoint = 0;

	savedPhi = new double* [xPointNum];
	for (int i = 0; i < xPointNum; i++)
	{
		savedPhi[i] = new double[yPointNum];
	}

	cellParticle = new double* [xPointNum-1];
	for (int i = 0; i < xPointNum-1; i++)
	{
		cellParticle[i] = new double[yPointNum-1];
	}

	for (int i = 0; i < xPointNum-1; i++)
	{
		for (int j = 0; j < yPointNum-1; j++)
		{
			cellParticle[i][j] = 0;
		}
	}



	PhiHead    = new Phi2D();
	PhiHead->x = inputX0;
	PhiHead->y = inputY0;
	PhiHead->xIndex = 0;
	PhiHead->yIndex = 0;
	PhiTail    = NULL;

	CellHead = NULL;
	CellTail = NULL;
	cellLevel = inputLevel;

	ParticleHead = NULL;
	ParticleTail = NULL;
}


Shape2D::~Shape2D()
{
	if (PhiHead != NULL)
	{
		deleteAllPhi(PhiTail);
	}

	if (ParticleHead != NULL)
	{
		deleteAllParticle(ParticleTail);
	}

	if (CellHead != NULL)
	{
		deleteAllCell(CellHead);
	}

	for (int i = 0; i < xPointNum; i++)
	{
		delete[] savedPhi[i];
	}
	delete[] savedPhi;

	for (int i = 0; i < xPointNum-1; i++)
	{
		delete[] cellParticle[i];
	}
	delete[] cellParticle;
}
