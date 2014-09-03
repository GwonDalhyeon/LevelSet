#include "functions.h"


int main()
{

	cout<< "Connecting phi on domain."<<endl;
	//Shape2D testShape2 = Shape2D(0,1,100,0,1,100,1.0,1.0/500.0);
	//testShape2.connectingPhi(testShape2.PhiHead,0,0);

	Shape2D testShape2 = Shape2D(0,1,0,1,6,1.0,1.0/500.0);
	testShape2.combiningCellPhi(NULL,0);
	testShape2.connectingPhi();


	cout<< "Initializing phi."<<endl;

	/////////////////////////////////////
	///////  Sing vortex
	////////////////////////////////////
	Phi2D* tempPhi = testShape2.PhiHead;
	double tempX, tempY;
	while (tempPhi!=NULL)
	{
		tempX = tempPhi->x;
		tempY = tempPhi->y;
		tempPhi->xVelocity = 2.0*sin(PI*(tempX))*sin(PI*(tempX))*sin(PI*(tempY))*cos(PI*(tempY));
		tempPhi->yVelocity = -2.0*sin(PI*(tempX))*cos(PI*(tempX))*sin(PI*(tempY))*sin(PI*(tempY));
		tempPhi->phi = sqrt( (tempX - 0.5)*(tempX - 0.5)+(tempY - 0.75)*(tempY - 0.75)) - 0.15;
		tempPhi = tempPhi->PhiNext;	
	}



	cout<< "Sprinkle particle."<<endl;
	testShape2.initializationParticle();

	cout<< "Start iterating."<<endl;
	testShape2.savingPhi(0);
	testShape2.savingParticle(0);
	for (int i = 1; i <= 1500; i++)
	{
		cout<< "Time step : " << i <<endl;
		testShape2.singleVortex();
		testShape2.singleVortexParticle();
		//testShape2.reductionError();

		//int j=0;
		//while (j<10)
		//{
		//	testShape2.reinitialTVDRK3();
		//	j=j+1;
		//}

		testShape2.savingPhi(i);
		testShape2.savingParticle(i);
		if (i%10==0)
		{
			int j=0;
			while (j<10)
			{
				testShape2.reinitialTVDRK3();
				j=j+1;
			}
			
		}
	}

	return 0;
}