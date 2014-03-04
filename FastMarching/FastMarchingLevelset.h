#include <iostream>

using namespace std;

struct LevelSet
{
	int* pointIndex;
	int xNum, yNum;
	double* xCood;
	double* yCood;
	double* Phi;
	double dx, dy;
	int* pointStatus;
	int alive, narrowBand, farAway;
	int aliveCount, narrowBandCount, farAwayCount;
	LevelSet();
	LevelSet(int gridNum, double dX);
	LevelSet(int gridNum, double dX, double* functionPhi);
	LevelSet(int gridxNum, int gridyNum, double dX, double dY);
	LevelSet(int gridxNum, int gridyNum, double dX, double dY, double* functionPhi);

	LevelSet operator=(LevelSet orignal_LevelSet)
	{
		xNum = orignal_LevelSet.xNum;
		yNum = orignal_LevelSet.yNum;
		dx = orignal_LevelSet.dx;
		dy = orignal_LevelSet.dy;
		alive =orignal_LevelSet.alive;
		narrowBand = orignal_LevelSet.narrowBand;
		farAway = orignal_LevelSet.farAway;


		pointIndex = new int [xNum*yNum];
		pointStatus = new int [xNum*yNum];
		xCood = new double [xNum*yNum];
		yCood = new double [xNum*yNum];
		Phi = new double [xNum*yNum];

		for (int i = 0; i < xNum; i++)
		{
			for (int j = 0; j < yNum; j++)
			{
				pointIndex[i+j*xNum]= orignal_LevelSet.pointIndex[i+j*xNum];
				xCood[i+j*xNum] = orignal_LevelSet.xCood[i+j*xNum];
				yCood[i+j*xNum] = orignal_LevelSet.yCood[i+j*xNum];
				pointStatus[i+j*xNum] = orignal_LevelSet.pointStatus[i+j*xNum];
			}
		}
		return *this;
	}

	double derivation2x(int i,int j, double* A);
	double derivation2y(int i,int j, double* A);
	double derivation1x_plus(int i,int j, double* A);
	double derivation1y_plus(int i,int j, double* A);
	double derivation1x_minus(int i,int j, double* A);
	double derivation1y_minus(int i,int j, double* A);
};

LevelSet :: LevelSet()
{
	dx = 0.0; dy = 0.0; pointStatus = NULL;
	xNum=0; yNum=0; pointIndex=NULL; xCood=NULL; yCood=NULL; Phi = NULL;
	alive=-1, narrowBand=0, farAway=1;
	aliveCount=0, narrowBandCount=0, farAwayCount=xNum*yNum;
}

LevelSet::LevelSet(int gridNum, double dX)
{
	alive=-1, narrowBand=0, farAway=1;
	dx = dX;
	dy = dX;
	xNum = gridNum; yNum = gridNum;
	pointIndex = new int [xNum*yNum];
	xCood = new double [xNum*yNum];
	yCood = new double [xNum*yNum];
	Phi = new double [xNum*yNum];
	pointStatus = new int[xNum*yNum];
	aliveCount=0, narrowBandCount=0, farAwayCount=xNum*yNum;
	for (int i = 0; i < xNum; i++)
	{
		for (int j = 0; j < yNum; j++)
		{
			pointIndex[i+j*xNum]=i+j*xNum;
			xCood[i+j*xNum] = -dx/2*((double)xNum-1) + dx* ((double)i);
			yCood[i+j*xNum] = -dx/2*((double)yNum-1) + dx* ((double)j);
			pointStatus[i+j*xNum] = farAway;
		}
	}
}

LevelSet::LevelSet(int gridNum, double dX, double* functionPhi)
{
	dx = dX;
	dy = dX;
	alive=-1, narrowBand=0, farAway=1;
	xNum = gridNum; yNum = gridNum;
	pointIndex = new int [xNum*yNum];
	xCood = new double [xNum*yNum];
	yCood = new double [xNum*yNum];
	Phi = new double [xNum*yNum];
	pointStatus = new int[xNum*yNum];
	aliveCount=0, narrowBandCount=0, farAwayCount=xNum*yNum;

	for (int i = 0; i < xNum; i++)
	{
		for (int j = 0; j < yNum; j++)
		{
			pointIndex[i+j*xNum]=i+j*xNum;
			xCood[i+j*xNum] = -dx/2*((double)xNum-1) + dx* ((double)i);
			yCood[i+j*xNum] = -dx/2*((double)yNum-1) + dx* ((double)j);
			Phi[i+j*xNum]= functionPhi[i+j*xNum];
			pointStatus[i+j*xNum] = farAway;
		}
	}
}

LevelSet::LevelSet(int gridxNum, int gridyNum, double dX, double dY)
{
	dx = dX;
	dy = dY;
	alive=-1, narrowBand=0, farAway=1;
	xNum = gridxNum; yNum = gridyNum;
	pointIndex = new int [xNum*yNum];
	xCood = new double [xNum*yNum];
	yCood = new double [xNum*yNum];
	Phi = new double [xNum*yNum];
	pointStatus = new int[xNum*yNum];
	aliveCount=0, narrowBandCount=0, farAwayCount=xNum*yNum;

	for (int i = 0; i < xNum; i++)
	{
		for (int j = 0; j < yNum; j++)
		{
			pointIndex[i+j*xNum]=i+j*xNum;
			xCood[i+j*xNum] = -dx/2*((double)xNum-1) + dx* ((double)i);
			yCood[i+j*xNum] = -dx/2*((double)yNum-1) + dx* ((double)j);
			pointStatus[i+j*xNum] = farAway;
		}
	}
}

LevelSet::LevelSet(int gridxNum, int gridyNum, double dX, double dY, double* functionPhi)
{
	dx = dX;
	dy = dY;
	alive=-1, narrowBand=0, farAway=1;
	xNum = gridxNum; yNum = gridyNum;
	pointIndex = new int [xNum*yNum];
	xCood = new double [xNum*yNum];
	yCood = new double [xNum*yNum];
	Phi = new double [xNum*yNum];
	pointStatus = new int[xNum*yNum];
	aliveCount=0, narrowBandCount=0, farAwayCount=xNum*yNum;
	for (int i = 0; i < xNum; i++)
	{
		for (int j = 0; j < yNum; j++)
		{
			pointIndex[i+j*xNum]=i+j*xNum;
			xCood[i+j*xNum] = -dx/2*((double)xNum-1) + dx* ((double)i);
			yCood[i+j*xNum] = -dx/2*((double)yNum-1) + dx* ((double)j);
			Phi[i+j*xNum]= functionPhi[i+j*xNum];
			pointStatus[i+j*xNum] = farAway;
		}
	}
}




#ifndef HeapsortAlgorithm_heapsort_h
#define HeapsortAlgorithm_heapsort_h
#include <iostream>
#include <cmath>
#include <stdio.h>
using namespace std;

struct HeapsortTree
{
	int length;
	int* indexArray;
	double* dataArray;
	bool* flagArray;

	HeapsortTree();
	HeapsortTree(int l);
	void deleteHeap(struct HeapsortTree* inputHeap);

	HeapsortTree operator=(HeapsortTree& originHeap)
	{
		length = NULL;
		indexArray = NULL;
		dataArray = NULL;
		flagArray = NULL;

		length = originHeap.length;
		indexArray = new int [length];
		dataArray = new double [length];
		flagArray = new bool [length];

		for (int i = 0; i < length; i++)
		{
			indexArray[i]=originHeap.indexArray[i];
			dataArray[i]=originHeap.dataArray[i];
			flagArray[i]=originHeap.flagArray[i];
		}
	}
	//    HeapsortTree(struct HeapsortTree *inputHeap);
	//    struct HeapsortTree *low;
	//    void makingTree(struct HeapsortTree inputheap);

};

HeapsortTree :: HeapsortTree()
{
	length = 0;
	indexArray= NULL;
	dataArray = NULL;
	flagArray = NULL;
}




HeapsortTree :: HeapsortTree(int l)
{
	length =l; //(int) sizeof(inputArray)/2;
	indexArray = new int [length];
	dataArray = new double [length];
	flagArray = new bool [length];
}



void HeapsortTree :: deleteHeap(struct HeapsortTree* inputHeap)
{
	inputHeap->dataArray=NULL;
	inputHeap->flagArray=NULL;
	inputHeap->indexArray = NULL;
}




struct HeapsortTree* checkingUpSorting(struct HeapsortTree *inputHeap,int child)
{
	int parent = floor((child-1)/2);
	if (parent<0)
	{
		return NULL;
	}


	if (inputHeap->dataArray[parent]>inputHeap->dataArray[child])
	{
		double tempData = inputHeap->dataArray[child];
		double tempIndex = inputHeap->indexArray[child];
		inputHeap->dataArray[child] = inputHeap->dataArray[parent];
		inputHeap->dataArray[parent] = tempData;
		inputHeap->indexArray[child] = inputHeap->indexArray[parent];
		inputHeap->indexArray[parent] = tempIndex;
		inputHeap->flagArray[child] = true;
		inputHeap->flagArray[parent] = false;

		checkingUpSorting(inputHeap, parent);
		return inputHeap ;
	}
	return NULL;
}



struct HeapsortTree* HeapInitial(double* inputArray, int arraySize)
{
	struct HeapsortTree* temp= new HeapsortTree(arraySize);
	temp->dataArray[0] = inputArray[0];
	temp->indexArray[0]=0;
	temp->flagArray[0]= true;

	for (int i = 0; i < temp->length; i++)
	{
		if (2*i+1<temp->length)
		{
			temp->dataArray[2*i+1] = inputArray[2*i+1];
			temp->indexArray[2*i+1] = 2*i+1;
			temp->flagArray[2*i+1] = false;

			if (temp->dataArray[i]>temp->dataArray[2*i+1])
			{
				double tempData = temp->dataArray[2*i+1];
				int tempIndex = temp->indexArray[i];
				temp->dataArray[2*i+1] = temp->dataArray[i];
				temp->dataArray[i] = tempData;
				temp->indexArray[i] = temp->indexArray[2*i+1];
				temp->indexArray[2*i+1] = tempIndex;
				temp->flagArray[2*i+1]=true;
				temp->flagArray[i] = false;

				checkingUpSorting(temp,i);
			}

		}
		if (2*i+2<temp->length) 
		{
			temp->dataArray[2*i+2] = inputArray[2*i+2];
			temp->indexArray[2*i+2] = 2*i+2;
			if (temp->dataArray[i]>temp->dataArray[2*i+2])
			{
				double tempData = temp->dataArray[2*i+2];
				int tempIndex = temp->indexArray[i];
				temp->dataArray[2*i+2] =temp-> dataArray[i];
				temp->dataArray[i] = tempData;
				temp->indexArray[i] = temp->indexArray[2*i+2];
				temp->indexArray[2*i+2] = tempIndex;
				temp->flagArray[2*i+2] = true;
				temp->flagArray[i] = false;

				checkingUpSorting(temp,i);
			}
		}
	}

	return temp;
}


struct HeapsortTree* HeapRearrange(struct HeapsortTree *inputHeap,int place)
{
	if (2*place+1>=inputHeap->length)
	{
		inputHeap->flagArray[place] = true;
		return inputHeap;
	}

	if (2*place+2 < inputHeap->length)
	{
		if (inputHeap->dataArray[2*place+1]>inputHeap->dataArray[2*place+2])
		{
			if (inputHeap->dataArray[2*place+2]<inputHeap->dataArray[place])
			{
				double tempData = inputHeap->dataArray[2*place+2];
				int tempIndex = inputHeap->indexArray[2*place+2];
				inputHeap->dataArray[2*place+2] =inputHeap-> dataArray[place];
				inputHeap->dataArray[place] = tempData;
				inputHeap->indexArray[2*place+2] = inputHeap->indexArray[2*place+2];
				inputHeap->indexArray[place] = tempIndex;
				inputHeap->flagArray[2*place+2] = false;
				inputHeap->flagArray[place] = true;

				HeapRearrange(inputHeap,2*place+2);
			}
			
		}
		else
		{
			if (inputHeap->dataArray[2*place+1]<inputHeap->dataArray[place])
			{
				double tempData = inputHeap->dataArray[2*place+1];
				int tempIndex = inputHeap->indexArray[2*place+1];
				inputHeap->dataArray[2*place+1] = inputHeap->dataArray[place];
				inputHeap->dataArray[place] = tempData;
				inputHeap->indexArray[2*place+1] = inputHeap->indexArray[place];
				inputHeap->indexArray[place] = tempIndex;
				inputHeap->flagArray[2*place+1]=false;
				inputHeap->flagArray[place] = true;

				HeapRearrange(inputHeap,2*place+1);
			}
			
		}
	}
	else
	{
		if (inputHeap->dataArray[2*place+1]<inputHeap->dataArray[place])
		{

			double tempData = inputHeap->dataArray[2*place+1];
			int tempIndex = inputHeap->indexArray[2*place+1];
			inputHeap->dataArray[2*place+1] = inputHeap->dataArray[place];
			inputHeap->dataArray[place] = tempData;
			inputHeap->indexArray[2*place+1] = inputHeap->indexArray[place];
			inputHeap->indexArray[place] = tempIndex;
			inputHeap->flagArray[2*place+1]=false;
			inputHeap->flagArray[place] = true;
			for (int i = 0; i < inputHeap->length; i++)
			{
				cout <<inputHeap->indexArray[i]<<" = "<< inputHeap->dataArray[i]<<"\n";
			}
			cout<<"\n";
			HeapRearrange(inputHeap,2*place+1);
		}
	}

	return inputHeap;
}


struct HeapsortTree* diminish(struct HeapsortTree* originHeap)
{
	struct HeapsortTree* resultHeap = new HeapsortTree(originHeap->length-1);

	resultHeap->length = originHeap->length-1;
	resultHeap->indexArray = new int [originHeap->length-1];
	resultHeap->dataArray = new double [originHeap->length-1];
	resultHeap->flagArray = new bool [originHeap->length-1];

	resultHeap->indexArray[0]=originHeap->indexArray[resultHeap->length];
	resultHeap->dataArray[0]=originHeap->dataArray[resultHeap->length];
	resultHeap->flagArray[0]=originHeap->flagArray[resultHeap->length];

	for (int i = 1; i < resultHeap->length; i++)
	{
		resultHeap->indexArray[i]=originHeap->indexArray[i];
		resultHeap->dataArray[i]=originHeap->dataArray[i];
		resultHeap->flagArray[i]=originHeap->flagArray[i];
	}

	return resultHeap;
}


struct HeapsortTree* diminish(struct HeapsortTree* originHeap, int level)
{
	struct HeapsortTree* resultHeap = new HeapsortTree(originHeap->length-1);

	resultHeap->length = level;
	resultHeap->indexArray = new int [level];
	resultHeap->dataArray = new double [level];
	resultHeap->flagArray = new bool [level];

	resultHeap->indexArray[0]=originHeap->indexArray[originHeap->length-1];
	resultHeap->dataArray[0]=originHeap->dataArray[originHeap->length-1];
	resultHeap->flagArray[0]=originHeap->flagArray[originHeap->length-1];

	for (int i = 1; i < level; i++)
	{
		resultHeap->indexArray[i]=originHeap->indexArray[i-1 + originHeap->length - level];
		resultHeap->dataArray[i]=originHeap->dataArray[i-1 + originHeap->length - level];
		resultHeap->flagArray[i]=originHeap->flagArray[i-1 + originHeap->length - level];
	}

	return resultHeap;
}



struct HeapsortTree* HeapExtract(struct HeapsortTree *inputHeap, int level)
{
	if (level==0)
	{
		return NULL;
	}

	if (level == inputHeap->length)
	{
		inputHeap->dataArray[0] = inputHeap->dataArray[0];
		inputHeap->flagArray[0] = inputHeap->flagArray[0];
		inputHeap->indexArray[0] = inputHeap->indexArray[0];

		HeapExtract(inputHeap,level-1);
	}
	else
	{
		struct HeapsortTree* tempHeap = diminish(inputHeap, level);
		struct HeapsortTree* tempResult = HeapRearrange(tempHeap,0);
		for (int i = 0; i < tempResult->length; i++)
		{
			inputHeap->dataArray[inputHeap->length - level+i] = tempResult->dataArray[i];
			inputHeap->flagArray[inputHeap->length - level+i] = tempResult->flagArray[i];
			inputHeap->indexArray[inputHeap->length - level+i] = tempResult->indexArray[i];
		}

		HeapExtract(inputHeap,level-1);
	}

	return inputHeap;
}

#endif