#include "flux.h"


int main()
{
	//Shape1D testShape = Shape1D(0.0, 1.0, 5, 1.0, 1e-3);

	//testShape.connectingPhi(testShape.Phi0,0);

	//cout<<testShape.Phi0->PhiRight->PhiRight->PhiLeft->x<<"\n";

	Shape2D testShape2 = Shape2D(0,1,100,0,1,100,1.0,1e-3);
	testShape2.connectingPhi(testShape2.PhiHead,0,0);
	cout<<testShape2.numPoint<<endl;
	//deleteAllPhi(testShape.Phi0);
	return 0;
}