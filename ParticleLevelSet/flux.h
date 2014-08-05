#include "shapeParticle.h"

// Thr Lax-Friedrichs flux, for example3 on
// "WENO scheme for Hamilton-Jacobi equations" - Jiang and Peng
double Shape2D::example3LaxFriedrichsFlux(double uPlus, double uMinus, double vPlus, double vMinus)
{
	
	double tempAlpha1 = 0.0, tempAlpha2 = 0.0;
	double tempBeta1 = 0.0, tempBeta2 = 0.0;
	double minU = min(uPlus, uMinus);
	double maxU = max(uPlus, uMinus);
	double minV = min(vPlus, vMinus);
	double maxV = max(vPlus, vMinus);

	double diffOfU = maxU - minU;
	double diffOfV = maxV - minV;
	double i = 0.0, j = 0.0;

	while (maxU > minU + deltaX*i)
	{
		j = 0.0;
		while (maxV > minV + deltaY*j )
		{
			tempAlpha2 = abs(cos(minU + deltaX*i + minV + deltaY*j));
			tempBeta2  = abs(cos(minU + deltaX*i + minV + deltaY*j));;
			if (tempAlpha1<tempAlpha2)
			{
				tempAlpha1 = tempAlpha2;
			}

			if (tempBeta1<tempBeta2)
			{
				tempBeta1 = tempBeta2;
			}

			j = j + 1.0;
		}
		i = i +1.0;
	}

	double return1 = sin((uPlus + uMinus)/2.0 + (vPlus +vMinus)/2.0);

	return return1 -tempAlpha1*(uPlus + uMinus)/2.0  - tempBeta1*(vPlus -vMinus)/2.0;
}




