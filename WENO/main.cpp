#include "flux.h"


int main()
{
	//Shape1D testShape = Shape1D(0.0, 1.0, 5, 1.0, 1e-3);

	//testShape.connectingPhi(testShape.Phi0,0);

	//cout<<testShape.Phi0->PhiRight->PhiRight->PhiLeft->x<<"\n";
	
	cout<< "Connecting phi on domain."<<endl;
	Shape2D testShape2 = Shape2D(-1,1,100,-1,1,100,1.0,1.0/100.0);
	testShape2.connectingPhi(testShape2.PhiHead,0,0);
	
	cout<< "Initializing phi"<<endl;
	// Sphere
	//Phi2D* tempPhi = testShape2.PhiHead;
	//double tempX, tempY;
	//double tempValue;
	//while (tempPhi!=NULL)
	//{
	//	tempX = tempPhi->x;
	//	tempY = tempPhi->y;
	//	tempValue = sqrt( (tempX - 0.5)*(tempX - 0.5)+(tempY - 0.5)*(tempY - 0.5)) - 0.3;
	//	tempPhi->phi = sqrt( (tempX - 0.5)*(tempX - 0.5)+(tempY - 0.5)*(tempY - 0.5)) - 0.3;
	//	tempPhi = tempPhi->PhiNext;
	//}
	

	 //A square.
	Phi2D* tempPhi = testShape2.PhiHead;
	double tempX, tempY;
	while (tempPhi!=NULL)
	{
		tempX = tempPhi->x;
		tempY = tempPhi->y;
		if (abs(tempX) < 0.5 && abs(tempY)<0.5)
			{
				tempPhi->phi = -1.0;
			}
			else
			{
				tempPhi->phi = 1.0;

			}
		tempPhi = tempPhi->PhiNext;
	}
			



	
	cout<< "Start iterating"<<endl;
	testShape2.savingPhi(0);
	for (int i = 1; i < 100; i++)
	{
		testShape2.reinitialTVDRK3();
		testShape2.savingPhi(i);
	}

	return 0;
}