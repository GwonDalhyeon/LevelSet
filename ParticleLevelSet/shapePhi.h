#include "shape.h"


Phi2D* Shape2D::findPoint(Phi2D* inputPhi, int inputXIndex, int inputYIndex)
{
	Phi2D* tempPhi = inputPhi;

	while (tempPhi!=NULL)
	{
		if (abs(tempPhi->x - (double)(inputXIndex*deltaX+X0))<DBL_EPSILON && abs(tempPhi->y - (double)(inputYIndex*deltaY+Y0))<DBL_EPSILON)
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
		if (abs(tempPhi->x - (double)(inputXIndex*deltaX+X0))<DBL_EPSILON && abs(tempPhi->y - (double)(inputYIndex*deltaY+Y0))<DBL_EPSILON)
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

	//returnPhi->xIndex = (int)floor( (inputX+X0)/deltaX+DBL_EPSILON);
	//returnPhi->yIndex = (int)floor( (inputY+Y0)/deltaY+DBL_EPSILON);
	numPoint = numPoint+1;

	//size_t currentSize = getPeakRSS();
	//size_t peakSize = getCurrentRSS();
	//cout <<numPoint+1 << "\n";
	//cout << currentSize << "   " << peakSize << endl;
	return returnPhi;
}



//void Shape2D::connectingPhi(Phi2D* inputPhi, int inputXIndex, int inputYIndex)
//{
//	if (findPointBool(PhiHead, inputXIndex+1, inputYIndex))
//	{
//		inputPhi->PhiRight = findPoint(PhiHead,inputXIndex+1,inputYIndex);
//	}
//	else if (inputPhi->PhiRight == NULL && !findPointBool(PhiHead, inputXIndex+1, inputYIndex))
//	{
//		Phi2D* newRightPhi   = NULL;
//		newRightPhi   = addPhi((double)(inputXIndex+1)*deltaX+X0, (double)inputYIndex*deltaY+Y0, inputPhi, NULL, NULL, NULL);
//
//		if (newRightPhi != NULL)
//		{
//			newRightPhi->xIndex = inputXIndex+1;
//			newRightPhi->yIndex = inputYIndex;
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
//		}
//	}
//
//
//	if (findPointBool(PhiHead, inputXIndex, inputYIndex+1))
//	{
//		inputPhi->PhiTop = findPoint(PhiHead,inputXIndex,inputYIndex+1);
//	}
//	else if (inputPhi->PhiTop == NULL && !findPointBool(PhiHead, inputXIndex, inputYIndex+1))
//	{
//		Phi2D* newTopPhi     = NULL;
//		newTopPhi     = addPhi((double)(inputXIndex)*deltaX+X0, (double)(inputYIndex+1)*deltaY+Y0, NULL, NULL, inputPhi, NULL);
//
//		if (newTopPhi != NULL)
//		{
//			newTopPhi->xIndex = inputXIndex;
//			newTopPhi->yIndex = inputYIndex+1;
//
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
//		}
//	}
//
//
//
//	if (findPointBool(PhiHead, inputXIndex-1, inputYIndex))
//	{
//		inputPhi->PhiLeft = findPoint(PhiHead,inputXIndex-1,inputYIndex);
//	}
//	else if (inputPhi->PhiLeft == NULL && !findPointBool(PhiHead, inputXIndex-1, inputYIndex))
//	{
//		Phi2D* newLeftPhi    = NULL;
//		newLeftPhi    = addPhi((double)(inputXIndex-1)*deltaX+X0, (double)inputYIndex*deltaY+Y0, NULL, inputPhi, NULL, NULL);
//
//		if (newLeftPhi != NULL)
//		{
//			newLeftPhi->xIndex = inputXIndex-1;
//			newLeftPhi->yIndex = inputYIndex;
//
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
//		}
//	}
//
//
//
//	if (findPointBool(PhiHead, inputXIndex, inputYIndex-1))
//	{
//		inputPhi->PhiBottom = findPoint(PhiHead,inputXIndex,inputYIndex-1);
//	}
//	else if (inputPhi->PhiBottom == NULL && !findPointBool(PhiHead, inputXIndex, inputYIndex-1))
//	{
//		Phi2D* newBottomPhi  = NULL;
//		newBottomPhi  = addPhi((double)(inputXIndex)*deltaX+X0, (double)(inputYIndex-1)*deltaY+X0, NULL, NULL, NULL, inputPhi);
//		if (newBottomPhi != NULL)
//		{
//			newBottomPhi->xIndex = inputXIndex;
//			newBottomPhi->yIndex = inputYIndex-1;
//
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
//	}
//}


void Shape2D::connectingPhi()
{
	Phi2D* tempPhi = PhiHead;
	while (tempPhi != NULL)
	{
		tempPhi->PhiRight  = findPoint(PhiHead,tempPhi->xIndex+1,tempPhi->yIndex);
		tempPhi->PhiTop    = findPoint(PhiHead,tempPhi->xIndex,tempPhi->yIndex+1);
		tempPhi->PhiLeft   = findPoint(PhiHead,tempPhi->xIndex-1,tempPhi->yIndex);
		tempPhi->PhiBottom = findPoint(PhiHead,tempPhi->xIndex,tempPhi->yIndex-1);

		tempPhi = tempPhi->PhiNext;
	}

}

void Shape2D::savingPhi(int timeIndex)
{
	Phi2D* tempPhi = PhiHead;
	while (tempPhi!=NULL)
	{
		savedPhi[tempPhi->xIndex][tempPhi->yIndex]=tempPhi->phi;
		tempPhi = tempPhi->PhiNext;
	}

	ofstream phiData;
	phiData.open("D:\Data/phi"+to_string(timeIndex)+".txt");

	for (int i = 0; i < xPointNum; i++)
	{
		for (int j = 0; j < yPointNum; j++)
		{
			phiData<<savedPhi[i][j]<<" ";
		}
		phiData<<endl;
	}
	phiData.close();
}

void Shape2D::indexingPhi(Phi2D* inputPhi)
{
	inputPhi->xIndex = floor((inputPhi->x-X0)/deltaX+DBL_EPSILON);
	inputPhi->yIndex = floor((inputPhi->y-Y0)/deltaY+DBL_EPSILON);
}


void countPositvePhiAroundCell(Cell2D* inputCell, int positvePoint,int bottomZero, int rightZero, int leftZero, int topZero)
{
	if (inputCell->PhiLeftBottom->phi * inputCell->PhiRightBottom->phi<0)
	{
		bottomZero = 0;
	}
	if (inputCell->PhiRightBottom->phi * inputCell->PhiRightTop->phi<0)
	{
		rightZero = 0;
	}
	if (inputCell->PhiLeftTop->phi * inputCell->PhiRightTop->phi<0)
	{
		topZero = 0;
	}
	if (inputCell->PhiLeftBottom->phi * inputCell->PhiLeftTop->phi<0)
	{
		leftZero = 0;
	}

	if (inputCell->PhiLeftBottom->phi > 0)
	{
		positvePoint = positvePoint + 1;

	}
	if (inputCell->PhiRightBottom->phi > 0)
	{
		positvePoint = positvePoint + 1;
	}
	if (inputCell->PhiLeftTop->phi > 0)
	{
		positvePoint = positvePoint + 1;
	}
	if (inputCell->PhiLeftBottom->phi > 0)
	{
		positvePoint = positvePoint + 1;
	}
	return;
}


void Shape2D::findZeroLevelSetPoint(Shape2D* inputShape2D)
{
	Cell2D* tempCell = CellHead;

	//double tempPhiLeftBottom;
	//double tempPhiLeftTop;
	//double tempPhiRightBottom;
	//double tempPhiRightTop;

	//double zeroBottomSign;
	//double zeroRightSign;
	//double zeroLeftSign;
	//double zeroTopSign;

	//double zeroCenterDiagSign;
	//double zeroCenteOffDiagSign;
	int numPositvePoint, bottomZero, rightZero, leftZero, topZero;
	
	while (tempCell != NULL)
	{
		if ((int) tempCell->CellLevel == (int) cellLevel)
		{
			//tempPhiLeftBottom  = tempCell->PhiLeftBottom->phi;
			//tempPhiLeftTop     = tempCell->PhiLeftTop->phi;
			//tempPhiRightBottom = tempCell->PhiRightBottom->phi;
			//tempPhiRightTop    = tempCell->PhiRightTop->phi;

			//zeroBottomSign = tempPhiLeftBottom*tempPhiRightBottom;
			//zeroRightSign  = tempPhiRightBottom*tempPhiRightTop;
			//zeroLeftSign   = tempPhiLeftBottom*tempPhiLeftTop;
			//zeroTopSign    = tempPhiLeftTop*tempPhiRightTop;
			//zeroCenterDiagSign = tempPhiLeftBottom*tempPhiRightTop;
			//zeroCenteOffDiagSign = tempPhiLeftTop*tempPhiRightBottom;
			
			numPositvePoint = 0;
			bottomZero = 0;
			rightZero = 0;
			leftZero = 0;
			topZero = 0;

			countPositvePhiAroundCell(tempCell, numPositvePoint, bottomZero, rightZero, leftZero, topZero);

			if (numPositvePoint > 0 && numPositvePoint<4)
			{
				addZeroLevelSetPoint(tempCell, numPositvePoint);
			}
			
		}


		tempCell = tempCell->CellNext;
	}
}


void Shape2D::addZeroLevelSetPoint(Cell2D* inputCell, int positivePhi)
{
	double tempPhiLeftBottom;
	double tempPhiLeftTop;
	double tempPhiRightBottom;
	double tempPhiRightTop;

	double zeroBottomSign;
	double zeroRightSign;
	double zeroLeftSign;
	double zeroTopSign;
	double zeroCenterDiagSign;
	double zeroCenteOffDiagSign;

	tempPhiLeftBottom  = inputCell->PhiLeftBottom->phi;
	tempPhiLeftTop     = inputCell->PhiLeftTop->phi;
	tempPhiRightBottom = inputCell->PhiRightBottom->phi;
	tempPhiRightTop    = inputCell->PhiRightTop->phi;

	zeroBottomSign = tempPhiLeftBottom*tempPhiRightBottom;
	zeroRightSign  = tempPhiRightBottom*tempPhiRightTop;
	zeroLeftSign   = tempPhiLeftBottom*tempPhiLeftTop;
	zeroTopSign    = tempPhiLeftTop*tempPhiRightTop;
	zeroCenterDiagSign = tempPhiLeftBottom*tempPhiRightTop;
	zeroCenteOffDiagSign = tempPhiLeftTop*tempPhiRightBottom;

	if (positivePhi == 1)
	{

	}


	if (inputCell->PhiLeftBottom->phi * inputCell->PhiRightBottom->phi<0)
	{
		
	}
	if (inputCell->PhiRightBottom->phi * inputCell->PhiRightTop->phi<0)
	{
		
	}
	if (inputCell->PhiLeftTop->phi * inputCell->PhiRightTop->phi<0)
	{
		
	}
	if (inputCell->PhiLeftBottom->phi * inputCell->PhiLeftTop->phi<0)
	{
		
	}
}