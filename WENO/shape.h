#include "operator.h"

class Shape1D
{
public:
	double deltaT;
	double deltaX;
	double T, X0, X1;
	int xPointNum;

	Phi1D* Phi0;
	Phi1D* Phi1;

	Shape1D();
	Shape1D(double inputX0, double inputX1, int inputXPointNum, double inputT, double inputdeltaT);
	Shape1D(double inputX0, double inputX1, double inputDeltaX, double inputT, double inputdeltaT);
	~Shape1D();

	void connectingPhi(Phi1D* inputPhi, int inputXIndex);
	Phi1D* addPhi(double inputX, Phi1D* leftPhi, Phi1D* rightPhi);

private:

};

Shape1D::Shape1D()
{
	Phi0 = NULL;
	Phi1 = NULL;
}

Shape1D::Shape1D(double inputX0, double inputX1, int inputXPointNum,double inputT, double inputdeltaT)
{
	deltaT    = inputdeltaT;
	deltaX    = (inputX1 - inputX0)/(double)inputXPointNum;

	T  = inputT;
	X0 = inputX0;
	X1 = inputX1;

	xPointNum = inputXPointNum;

	Phi0 = new Phi1D();
	Phi0->x = inputX0;
	Phi1 = new Phi1D();
	Phi1->x = inputX1;
}

Shape1D::Shape1D(double inputX0, double inputX1, double inputDeltaX, double inputT, double inputdeltaT)
{
	deltaT    = inputdeltaT;
	deltaX    = inputDeltaX;

	T  = inputT;
	X0 = inputX0;
	X1 = inputX1;

	xPointNum = floor((inputX1 - inputX0)/inputDeltaX);

	Phi0 = new Phi1D();
	Phi0->x = inputX0;
	Phi1 = new Phi1D();
	Phi1->x = inputX1;
}

Shape1D::~Shape1D()
{
	if (Phi0 != NULL)
	{
		deleteAllPhi(Phi0);
	}
}

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

	Shape2D(double inputX0, double inputX1, int inputXPointNum, double inputY0, double inputY1, int inputYPointNum, double inputT, double inputdeltaT);
	Shape2D(double inputX0, double inputX1, double inputDeltaX, double inputY0, double inputY1, double inputDeltaY, double inputT, double inputdeltaT);
	~Shape2D();

	Phi2D* addPhi(double inputX, double inputY, Phi2D* leftPhi, Phi2D* rightPhi, Phi2D* lowerPhi, Phi2D* UpperPhi);
	void connectingPhi(Phi2D* inputPhi, int inputXIndex, int inputYIndex);
	Phi2D* findPoint(Phi2D* phiHead, int inputXIndex, int inputYIndex);
	bool findPointBool(Phi2D* phiHead, int inputXIndex, int inputYIndex);
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

	xPointNum = inputXPointNum;
	yPointNum = inputYPointNum;

	PhiHead    = new Phi2D();
	PhiHead->x = inputX0;
	PhiHead->y = inputY0;
	PhiTail    = NULL;
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

	xPointNum = floor((inputX1 - inputX0)/inputDeltaX);

	PhiHead    = new Phi2D();
	PhiHead->x = inputX0;
	PhiHead->y = inputY0;
	PhiTail    = NULL;

}


Shape2D::~Shape2D()
{
	if (PhiHead != NULL)
	{
		deleteAllPhi(PhiTail);
	}
}

Phi1D* Shape1D::addPhi(double inputX, Phi1D* leftPhi, Phi1D* rightPhi)
{
	Phi1D* returnPhi = new Phi1D();

	if (inputX <X0  || inputX>X1)
	{
		return NULL;
	}

	if (rightPhi != NULL)
	{
		rightPhi->PhiLeft = returnPhi;
	}

	if (leftPhi != NULL)
	{
		leftPhi->PhiRight = returnPhi;
	}

	returnPhi->x = inputX;

	return returnPhi;
}

void Shape1D::connectingPhi(Phi1D* inputPhi, int inputXIndex)
{
	Phi1D* newRightPhi = addPhi((double)(inputXIndex+1)*deltaX, inputPhi, NULL) ;

	if (newRightPhi != NULL)
	{
		newRightPhi->PhiLeft = inputPhi;
		connectingPhi(inputPhi->PhiRight, inputXIndex+1);
	}
}


Phi2D* Shape2D::findPoint(Phi2D* inputPhi, int inputXIndex, int inputYIndex)
{
	Phi2D* tempPhi = inputPhi;

	while (tempPhi!=NULL)
	{
		if (abs(tempPhi->x - (double)inputXIndex*deltaX)<DBL_EPSILON && abs(tempPhi->y - (double)inputYIndex*deltaY)<DBL_EPSILON)
		{
			return tempPhi;
		}
		tempPhi = tempPhi->PhiNext;
	}
	return NULL;
}

bool Shape2D::findPointBool(Phi2D* inputPhi, int inputXIndex, int inputYIndex)
{
	Phi2D* tempPhi = inputPhi;

	while (tempPhi!=NULL)
	{
		if (abs(tempPhi->x - (double)inputXIndex*deltaX)<DBL_EPSILON && abs(tempPhi->y - (double)inputYIndex*deltaY)<DBL_EPSILON)
		{
			return true;
		}
		tempPhi = tempPhi->PhiNext;
	}
	return false;
}


Phi2D* Shape2D::addPhi(double inputX, double inputY, Phi2D* leftPhi, Phi2D* rightPhi, Phi2D* bottomPhi, Phi2D* topPhi)
{
	Phi2D* returnPhi = new Phi2D();

	if (inputX <X0  || inputX>X1 || inputY<Y0 || inputY>Y1)
	{
		return NULL;
	}

	if (leftPhi != NULL && leftPhi->PhiRight==NULL)
	{
		leftPhi->PhiRight  = returnPhi;
		returnPhi->PhiLeft = leftPhi;
	}

	if (rightPhi != NULL && rightPhi->PhiLeft==NULL)
	{
		rightPhi->PhiLeft   = returnPhi;
		returnPhi->PhiRight = rightPhi;
	}

	if (bottomPhi != NULL && bottomPhi->PhiTop==NULL)
	{
		bottomPhi->PhiTop = returnPhi;
		returnPhi->PhiBottom = bottomPhi;
	}

	if (topPhi != NULL && topPhi->PhiBottom==NULL)
	{
		topPhi->PhiBottom = returnPhi;
		returnPhi->PhiTop = topPhi;
	}
	returnPhi->x = inputX;
	returnPhi->y = inputY;
	//numPoint = numPoint+1;

	//size_t currentSize = getPeakRSS();
	//	size_t peakSize = getCurrentRSS();
	//	cout <<numPoint+1 << "\n";
	//	cout << currentSize << "   " << peakSize << endl;
	return returnPhi;
}



void Shape2D::connectingPhi(Phi2D* inputPhi, int inputXIndex, int inputYIndex)
{
	if (findPointBool(PhiHead, inputXIndex+1, inputYIndex))
	{
		inputPhi->PhiRight = findPoint(PhiHead,inputXIndex+1,inputYIndex);
	}
	else if (inputPhi->PhiRight == NULL && !findPointBool(PhiHead, inputXIndex+1, inputYIndex))
	{
		Phi2D* newRightPhi   = NULL;
		newRightPhi   = addPhi((double)(inputXIndex+1)*deltaX, (double)inputYIndex*deltaY, inputPhi, NULL, NULL, NULL);

		if (newRightPhi != NULL)
		{
			Phi2D* oldPhiTail;
			oldPhiTail = PhiTail;

			if (oldPhiTail != NULL)
			{
				newRightPhi->PhiBefore = oldPhiTail;
				oldPhiTail->PhiNext = newRightPhi;
				PhiTail = newRightPhi;
			}
			else
			{
				PhiTail = newRightPhi;
				PhiHead->PhiNext = newRightPhi;
				PhiTail->PhiBefore = PhiHead;
			}

			connectingPhi(newRightPhi, inputXIndex+1, inputYIndex);
		}
	}
	

	if (findPointBool(PhiHead, inputXIndex, inputYIndex+1))
	{
		inputPhi->PhiTop = findPoint(PhiHead,inputXIndex,inputYIndex+1);
	}
	else if (inputPhi->PhiTop == NULL && !findPointBool(PhiHead, inputXIndex, inputYIndex+1))
	{
		Phi2D* newTopPhi     = NULL;
		newTopPhi     = addPhi((double)(inputXIndex)*deltaX, (double)(inputYIndex+1)*deltaY, NULL, NULL, inputPhi, NULL);

		if (newTopPhi != NULL)
		{
			Phi2D* oldPhiTail;
			oldPhiTail = PhiTail;

			if (oldPhiTail != NULL)
			{
				newTopPhi->PhiBefore = oldPhiTail;
				oldPhiTail->PhiNext = newTopPhi;
				PhiTail = newTopPhi;
			}
			else
			{
				PhiTail = newTopPhi;
				PhiHead->PhiNext = newTopPhi;
				PhiTail->PhiBefore = PhiHead;
			}

			connectingPhi(newTopPhi, inputXIndex, inputYIndex+1);
		}
	}
	


	if (findPointBool(PhiHead, inputXIndex-1, inputYIndex))
	{
		inputPhi->PhiLeft = findPoint(PhiHead,inputXIndex-1,inputYIndex);
	}
	else if (inputPhi->PhiLeft == NULL && !findPointBool(PhiHead, inputXIndex-1, inputYIndex))
	{
		Phi2D* newLeftPhi    = NULL;
		newLeftPhi    = addPhi((double)(inputXIndex-1)*deltaX, (double)inputYIndex*deltaY, NULL, inputPhi, NULL, NULL);

		if (newLeftPhi != NULL)
		{
			Phi2D* oldPhiTail;
			oldPhiTail = PhiTail;

			if (oldPhiTail != NULL)
			{
				newLeftPhi->PhiBefore = oldPhiTail;
				oldPhiTail->PhiNext = newLeftPhi;
				PhiTail = newLeftPhi;
			}
			else
			{
				PhiTail = newLeftPhi;
				PhiHead->PhiNext = newLeftPhi;
				PhiTail->PhiBefore = PhiHead;
			}

			connectingPhi(newLeftPhi, inputXIndex-1, inputYIndex);
		}
	}
	


	if (findPointBool(PhiHead, inputXIndex, inputYIndex-1))
	{
		inputPhi->PhiBottom = findPoint(PhiHead,inputXIndex-1,inputYIndex);
	}
	else if (inputPhi->PhiBottom == NULL && !findPointBool(PhiHead, inputXIndex, inputYIndex-1))
	{
		Phi2D* newBottomPhi  = NULL;
		newBottomPhi  = addPhi((double)(inputXIndex)*deltaX, (double)(inputYIndex-1)*deltaY, NULL, NULL, NULL, inputPhi);
		if (newBottomPhi != NULL)
		{

			Phi2D* oldPhiTail;
			oldPhiTail = PhiTail;

			if (oldPhiTail != NULL)
			{
				newBottomPhi->PhiBefore = oldPhiTail;
				oldPhiTail->PhiNext = newBottomPhi;
				PhiTail = newBottomPhi;
			}
			else
			{
				PhiTail = newBottomPhi;
				PhiHead->PhiNext = newBottomPhi;
				PhiTail->PhiBefore = PhiHead;
			}

			connectingPhi(newBottomPhi, inputXIndex, inputYIndex-1);
		}

	}
	

}




//void Shape2D::connectingPhi(Phi2D* inputPhi, int inputXIndex, int inputYIndex)
//{
//
//	Phi2D* newRightPhi   = NULL;
//	Phi2D* newLeftPhi    = NULL;
//	Phi2D* newTopPhi     = NULL;
//	Phi2D* newBottomPhi  = NULL;
//
//	if (inputPhi->PhiRight == NULL)
//	{
//		newRightPhi   = addPhi((double)(inputXIndex+1)*deltaX, (double)inputYIndex*deltaY, inputPhi, NULL, NULL, NULL);
//	}
//	if (inputPhi->PhiLeft == NULL)
//	{
//		newLeftPhi    = addPhi((double)(inputXIndex-1)*deltaX, (double)inputYIndex*deltaY, NULL, inputPhi, NULL, NULL);
//	}
//	if (inputPhi->PhiTop == NULL)
//	{
//		newTopPhi     = addPhi((double)(inputXIndex)*deltaX, (double)(inputYIndex+1)*deltaY, NULL, NULL, inputPhi, NULL);
//	}
//	if (inputPhi->PhiBottom == NULL)
//	{
//		newBottomPhi  = addPhi((double)(inputXIndex)*deltaX, (double)(inputYIndex-1)*deltaY, NULL, NULL, NULL, inputPhi);
//	}
//
//
//	if (newRightPhi != NULL && !findPointBool(PhiHead, inputXIndex+1, inputYIndex))
//	{
//		
//			Phi2D* oldPhiTail;
//			oldPhiTail = PhiTail;
//
//			if (oldPhiTail != NULL)
//			{
//				newRightPhi->PhiBefore = oldPhiTail;
//				oldPhiTail->PhiNext = newRightPhi;
//				PhiTail = newRightPhi;
//			}
//			else
//			{
//				PhiTail = newRightPhi;
//				PhiHead->PhiNext = newRightPhi;
//				PhiTail->PhiBefore = PhiHead;
//			}
//
//			connectingPhi(newRightPhi, inputXIndex+1, inputYIndex);
//	}
//
//	if (newTopPhi != NULL && !findPointBool(PhiHead, inputXIndex, inputYIndex+1))
//	{
//			Phi2D* oldPhiTail;
//			oldPhiTail = PhiTail;
//
//			if (oldPhiTail != NULL)
//			{
//				newTopPhi->PhiBefore = oldPhiTail;
//				oldPhiTail->PhiNext = newTopPhi;
//				PhiTail = newTopPhi;
//			}
//			else
//			{
//				PhiTail = newTopPhi;
//				PhiHead->PhiNext = newTopPhi;
//				PhiTail->PhiBefore = PhiHead;
//			}
//
//			connectingPhi(newTopPhi, inputXIndex, inputYIndex+1);
//	}
//
//	if (newLeftPhi != NULL && !findPointBool(PhiHead, inputXIndex-1, inputYIndex))
//	{
//			Phi2D* oldPhiTail;
//			oldPhiTail = PhiTail;
//
//			if (oldPhiTail != NULL)
//			{
//				newLeftPhi->PhiBefore = oldPhiTail;
//				oldPhiTail->PhiNext = newLeftPhi;
//				PhiTail = newLeftPhi;
//			}
//			else
//			{
//				PhiTail = newLeftPhi;
//				PhiHead->PhiNext = newLeftPhi;
//				PhiTail->PhiBefore = PhiHead;
//			}
//
//			connectingPhi(newLeftPhi, inputXIndex-1, inputYIndex);
//	}
//
//	if (newBottomPhi != NULL && !findPointBool(PhiHead, inputXIndex, inputYIndex-1))
//	{
//			Phi2D* oldPhiTail;
//			oldPhiTail = PhiTail;
//
//			if (oldPhiTail != NULL)
//			{
//				newBottomPhi->PhiBefore = oldPhiTail;
//				oldPhiTail->PhiNext = newBottomPhi;
//				PhiTail = newBottomPhi;
//			}
//			else
//			{
//				PhiTail = newBottomPhi;
//				PhiHead->PhiNext = newBottomPhi;
//				PhiTail->PhiBefore = PhiHead;
//			}
//
//			connectingPhi(newBottomPhi, inputXIndex, inputYIndex-1);
//		}
//}