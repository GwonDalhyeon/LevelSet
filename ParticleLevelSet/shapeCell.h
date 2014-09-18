#include "shapePhi.h"



Cell2D* Shape2D::addCell(Phi2D* leftBottomPhi, Phi2D* rightBottomPhi, Phi2D* leftTopPhi, Phi2D* RightTopPhi, double inputLevel)
{
	Cell2D* returnCell = new Cell2D();

	returnCell->CellLevel = inputLevel;

	returnCell->PhiLeftBottom = leftBottomPhi;
	returnCell->PhiRightBottom = rightBottomPhi;
	returnCell->PhiLeftTop = leftTopPhi;
	returnCell->PhiRightTop = RightTopPhi;

	returnCell->x0 = leftBottomPhi->x;
	returnCell->x1 = RightTopPhi->x;
	returnCell->y0 = leftBottomPhi->y;
	returnCell->y1 = RightTopPhi->y;
	

	returnCell->numContainParticle = 0;

	return returnCell;
}



void Shape2D::combiningCellPhi(Cell2D* inputCell, double inputLevel)
{
	if (inputCell == NULL)
	{
		inputCell = new Cell2D();
		inputCell->CellLevel = 0.0;
		inputCell->x0 = X0;
		inputCell->x1 = X1;
		inputCell->y0 = Y0;
		inputCell->y1 = Y1;

		inputCell->PhiLeftBottom = PhiHead;
		Phi2D* newRightBottomPhi = new Phi2D();
		Phi2D* newLeftTopPhi     = new Phi2D();
		Phi2D* newRightTopPhi    = new Phi2D();

		inputCell->PhiLeftBottom  = PhiHead;
		inputCell->PhiRightBottom = newRightBottomPhi;
		inputCell->PhiLeftTop     = newLeftTopPhi;
		inputCell->PhiRightTop    = newRightTopPhi;

		CellHead = inputCell;
		CellTail = inputCell;

		newRightBottomPhi->x = X1;
		newRightBottomPhi->y = Y0;
		newLeftTopPhi->x     = X0;
		newLeftTopPhi->y     = Y1;
		newRightTopPhi->x    = X1;
		newRightTopPhi->y    = Y1;

		PhiHead->PhiNext = newRightBottomPhi;
		newRightBottomPhi->PhiBefore = PhiHead;
		newRightBottomPhi->PhiNext   = newRightTopPhi;
		newRightTopPhi->PhiBefore    = newRightBottomPhi;
		newRightTopPhi->PhiNext      = newLeftTopPhi;
		newLeftTopPhi->PhiBefore     = newRightTopPhi;
		PhiTail = newLeftTopPhi;
		
		//PhiHead->PhiRight          = newRightBottomPhi;
		//PhiHead->PhiTop            = newRightTopPhi;
		//newRightBottomPhi->PhiLeft = PhiHead;
		//newRightBottomPhi->PhiTop  = newRightTopPhi;
		//newRightTopPhi->PhiBottom  = newRightBottomPhi;
		//newRightTopPhi->PhiLeft    = newLeftTopPhi;
		//newLeftTopPhi->PhiRight    = newRightTopPhi;
		//newLeftTopPhi->PhiBottom   = PhiHead;

		indexingPhi(newRightBottomPhi);
		indexingPhi(newRightTopPhi);
		indexingPhi(newLeftTopPhi);
	}


	if (abs(inputLevel - cellLevel)<DBL_EPSILON)
	{
		return;
	}

	


	Phi2D* newCenterPhi      = new Phi2D();
	Phi2D* newBottomCenterPhi= new Phi2D();
	Phi2D* newRightCenterPhi = new Phi2D();
	Phi2D* newTopCenterPhi   = new Phi2D();
	Phi2D* newLeftCenterPhi  = new Phi2D();

	newCenterPhi->x = (inputCell->x0 + inputCell->x1)/2.0;
	newCenterPhi->y = (inputCell->y0 + inputCell->y1)/2.0;

	newBottomCenterPhi->x = (inputCell->x0 + inputCell->x1)/2.0;
	newBottomCenterPhi->y = inputCell->y0;

	newLeftCenterPhi->x = inputCell->x0;
	newLeftCenterPhi->y = (inputCell->y0 + inputCell->y1)/2.0;


	newRightCenterPhi->x = inputCell->x1;
	newRightCenterPhi->y = (inputCell->y0 + inputCell->y1)/2.0;

	newTopCenterPhi->x = (inputCell->x0 + inputCell->x1)/2.0;
	newTopCenterPhi->y =inputCell->y1;
	
	indexingPhi(newCenterPhi);
	indexingPhi(newBottomCenterPhi);
	indexingPhi(newRightCenterPhi);
	indexingPhi(newTopCenterPhi);
	indexingPhi(newLeftCenterPhi);


	Phi2D* oldPhiTail = PhiTail;
	oldPhiTail->PhiNext = newCenterPhi;

	newCenterPhi->PhiBefore = oldPhiTail;
	newCenterPhi->PhiNext = newBottomCenterPhi;

	newBottomCenterPhi->PhiBefore = newCenterPhi;
	newBottomCenterPhi->PhiNext = newLeftCenterPhi;

	newLeftCenterPhi->PhiBefore = newBottomCenterPhi;
	newLeftCenterPhi->PhiNext = newRightCenterPhi;


	newRightCenterPhi->PhiBefore = newLeftCenterPhi;
	newRightCenterPhi->PhiNext = newTopCenterPhi;

	newTopCenterPhi->PhiBefore = newRightCenterPhi;
	PhiTail = newTopCenterPhi;




	Cell2D* newLeftBottomCell  = addCell(inputCell->PhiLeftBottom, newBottomCenterPhi, newLeftCenterPhi,newCenterPhi, inputLevel + 1.0);
	Cell2D* newRightBottomcell = addCell(newBottomCenterPhi, inputCell->PhiRightBottom, newCenterPhi, newRightCenterPhi, inputLevel + 1.0);
	Cell2D* newLeftTopCell     = addCell(newLeftCenterPhi, newCenterPhi, inputCell->PhiLeftTop, newTopCenterPhi, inputLevel + 1.0);
	Cell2D* newRightTopCell    = addCell(newCenterPhi, newRightCenterPhi, newTopCenterPhi, inputCell->PhiRightTop, inputLevel + 1.0);



	newLeftBottomCell->CellParent  = inputCell;
	newRightBottomcell->CellParent = inputCell;
	newLeftTopCell->CellParent     = inputCell;
	newRightTopCell->CellParent    = inputCell;

	inputCell->CellChildLeftBottom = newLeftBottomCell;
	inputCell->CellChildRightBottom = newRightBottomcell;
	inputCell->CellChildLeftTop = newLeftTopCell;
	inputCell->CellChildRightTop = newRightTopCell;



	Cell2D* oldCellTail = CellTail;
	oldCellTail->CellNext = newLeftBottomCell;

	newLeftBottomCell->CellBefore  = oldCellTail;
	newLeftBottomCell->CellNext  = newRightBottomcell;

	newRightBottomcell->CellBefore  = newLeftBottomCell;
	newRightBottomcell->CellNext  = newLeftTopCell;
	
	newLeftTopCell->CellBefore  = newRightBottomcell;
	newLeftTopCell->CellNext  = newRightTopCell;
	
	newRightTopCell->CellBefore  = newLeftTopCell;
	CellTail = newRightTopCell;



	combiningCellPhi(newLeftBottomCell, newLeftBottomCell->CellLevel);
	combiningCellPhi(newRightBottomcell, newRightBottomcell->CellLevel);
	combiningCellPhi(newLeftTopCell, newLeftTopCell->CellLevel);
	combiningCellPhi(newRightTopCell, newRightTopCell->CellLevel);
}