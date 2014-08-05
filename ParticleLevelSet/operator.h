#include "particles.h"
#include <algorithm>
#include <time.h>
#include <string>
#include <cstring>

#define CHECK_TIME_START __int64 freq, start, end; if (QueryPerformanceFrequency((_LARGE_INTEGER*)&freq)) {QueryPerformanceCounter((_LARGE_INTEGER*)&start);
// a는 float type milli second이고 b가 FALSE일때는 에러입니다
#define CHECK_TIME_END(a,b) QueryPerformanceCounter((_LARGE_INTEGER*)&end); a=(float)((double)(end - start)/freq*1000); b=TRUE; } else b=FALSE;


#define PI 3.141592653589793238

double sign(double a)
{
	if (a>0)
	{
		return 1.0;
	}
	else if (a<0)
	{
		return -1.0;
	}
	else
	{
		return 0;
	}
}

double sign2(double a)
{
	return a/sqrt(a*a + DBL_EPSILON*DBL_EPSILON);
}

double minmod(double a, double b)
{
	return (sign(a) + sign(b))*min(abs(a), abs(b))/2.0;
}

//////////////////////////////////////////////////////////////////////////////////////////////
////////////////
////////////////    2D
////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

double dxxPhi(Phi2D* inputPhi)
{
	if (inputPhi == NULL)
	{
		return 0.0;
	}
	else if ((inputPhi->PhiLeft == NULL) || (inputPhi->PhiRight == NULL))
	{
		return 0.0;
	}
	else
	{
		double dxLeft  = inputPhi->x - inputPhi->PhiLeft->x;
		double dxRight = inputPhi->PhiRight->x - inputPhi->x;
		return (inputPhi->PhiRight->phi - inputPhi->phi) / dxRight *2.0/(dxLeft + dxRight) - (inputPhi->phi - inputPhi->PhiLeft->phi) / dxLeft *2.0/(dxLeft + dxRight);
	}
}

double dxMinusPhi(Phi2D* inputPhi)
{
	if ((inputPhi->PhiLeft == NULL) || (inputPhi == NULL))
	{
		return 0.0;
	}
	else
	{
		double dxLeft  = inputPhi->x - inputPhi->PhiLeft->x;
		return (inputPhi->phi - inputPhi->PhiLeft->phi) / dxLeft + dxLeft/2.0*minmod(dxxPhi(inputPhi),dxxPhi(inputPhi->PhiLeft));
	}
}

double dxPlusPhi(Phi2D* inputPhi)
{
	if ((inputPhi == NULL) || (inputPhi->PhiRight == NULL))
	{
		return 0.0;
	}
	else
	{
		double dxRight = inputPhi->PhiRight->x - inputPhi->x;
		return (inputPhi->PhiRight->phi - inputPhi->phi) / dxRight - dxRight/2.0*minmod(dxxPhi(inputPhi),dxxPhi(inputPhi->PhiRight));
	}
}

double dxPhi(Phi2D* inputPhi)
{
	if (inputPhi == NULL)
	{
		return 0.0;
	}
	else if ((inputPhi->PhiLeft == NULL) && (inputPhi->PhiRight == NULL))
	{
		return 0.0;
	}
	else if (inputPhi->PhiLeft == NULL)
	{
		return dxPlusPhi(inputPhi);
	}
	else if (inputPhi->PhiRight == NULL)
	{
		return dxMinusPhi(inputPhi);
	}
	else
	{
		double dxLeft  = inputPhi->x - inputPhi->PhiLeft->x;
		double dxRight = inputPhi->PhiRight->x - inputPhi->x;
		return (inputPhi->phi - inputPhi->PhiLeft->phi) / dxLeft *dxRight/(dxLeft + dxRight) + (inputPhi->PhiRight->phi - inputPhi->phi) / dxRight *dxLeft/(dxLeft + dxRight);
	}

}

double dyyPhi(Phi2D* inputPhi)
{
	if (inputPhi == NULL)
	{
		return 0.0;
	}
	else if ((inputPhi->PhiTop == NULL) || (inputPhi->PhiBottom == NULL))
	{
		return 0.0;
	}
	else
	{
		double dyBottom  = inputPhi->y - inputPhi->PhiBottom->y;
		double dyTop     = inputPhi->PhiTop->y - inputPhi->y;
		return (inputPhi->PhiTop->phi - inputPhi->phi) / dyTop *2.0/(dyTop + dyBottom) - (inputPhi->phi - inputPhi->PhiBottom->phi) / dyBottom *2.0/(dyBottom + dyTop);
	}
}

double dyMinusPhi(Phi2D* inputPhi)
{
	if ((inputPhi == NULL) || (inputPhi->PhiBottom == NULL))
	{
		return 0.0;
	}
	else
	{
		double dyBottom = inputPhi->y - inputPhi->PhiBottom->y;
		return (inputPhi->phi - inputPhi->PhiBottom->phi) / dyBottom + dyBottom/2.0*minmod(dxxPhi(inputPhi),dxxPhi(inputPhi->PhiBottom));
	}
}

double dyPlusPhi(Phi2D* inputPhi)
{
	if ((inputPhi == NULL) || (inputPhi->PhiTop == NULL))
	{
		return 0.0;
	}
	else
	{
		double dyTop = inputPhi->PhiTop->y - inputPhi->y;
		return (inputPhi->PhiTop->phi - inputPhi->phi) / dyTop - dyTop/2.0*minmod(dxxPhi(inputPhi),dxxPhi(inputPhi->PhiTop));
	}
}

double dyPhi(Phi2D* inputPhi)
{
	if (inputPhi == NULL)
	{
		return 0.0;
	}
	else if ((inputPhi->PhiTop == NULL) && (inputPhi->PhiBottom == NULL))
	{
		return 0.0;
	}
	else if (inputPhi->PhiTop == NULL)
	{
		return dyMinusPhi(inputPhi);
	}
	else if (inputPhi->PhiBottom == NULL)
	{
		return dyPlusPhi(inputPhi);
	}
	else
	{
		double dyBottom  = inputPhi->y - inputPhi->PhiBottom->y;
		double dyTop     = inputPhi->PhiTop->y - inputPhi->y;
		return (inputPhi->phi - inputPhi->PhiBottom->phi) / dyBottom *dyTop/(dyBottom + dyTop) + (inputPhi->PhiTop->phi - inputPhi->phi) / dyTop *dyBottom/(dyTop + dyBottom);
	}
}

//double dxyPhi(Phi2D* inputPhi)
//{
//	double dx = inputPhi->PhiRight->x - inputPhi->PhiLeft->x;
//	double dy = inputPhi->PhiTop->y - inputPhi->PhiBottom->y;
//	return (inputPhi->PhiRight->PhiTop->phi - inputPhi->PhiRight->PhiBottom->phi - inputPhi->PhiLeft->PhiTop->phi - inputPhi->PhiLeft->PhiBottom->phi)/(4.0*dx*dy);
//}



double wenoApproximationXMinus(Phi2D* inputPhi)
{
	double v1; double v2; double v3; double v4; double v5;

	if (inputPhi == NULL || inputPhi->PhiRight == NULL || inputPhi->PhiRight->PhiRight == NULL || inputPhi->PhiLeft == NULL || inputPhi->PhiLeft->PhiLeft == NULL || inputPhi->PhiLeft->PhiLeft->PhiLeft == NULL)
	{
		return dxMinusPhi(inputPhi);
	}

	v1 = (inputPhi->PhiLeft->PhiLeft->phi - inputPhi->PhiLeft->PhiLeft->PhiLeft->phi) / (inputPhi->PhiLeft->PhiLeft->x - inputPhi->PhiLeft->PhiLeft->PhiLeft->x);
	v2 = (inputPhi->PhiLeft->phi - inputPhi->PhiLeft->PhiLeft->phi) / (inputPhi->PhiLeft->x - inputPhi->PhiLeft->PhiLeft->x);
	v3 = (inputPhi->phi - inputPhi->PhiLeft->phi) / (inputPhi->x - inputPhi->PhiLeft->x);
	v4 = (inputPhi->PhiRight->phi - inputPhi->phi) / (inputPhi->PhiRight->x - inputPhi->x);
	v5 = (inputPhi->PhiRight->PhiRight->phi - inputPhi->PhiRight->phi) / (inputPhi->PhiRight->PhiRight->x - inputPhi->PhiRight->x);

	double s1 = 13.0/12.0*(v1-2.0*v2+v3)*(v1-2.0*v2+v3) + 1.0/4.0*(v1-4.0*v2+3.0*v3)*(v1-4.0*v2+3.0*v3);
	double s2 = 13.0/12.0*(v2-2.0*v3+v4)*(v2-2.0*v3+v4) + 1.0/4.0*(v2-v4)*(v2-v4);
	double s3 = 13.0/12.0*(v3-2.0*v4+v5)*(v3-2.0*v4+v5) + 1.0/4.0*(3.0*v3-4.0*v4+v5)*(3.0*v3-4.0*v4+v5);

	double eps = 1e-6;

	double a1 = 1.0/10.0 * 1.0/((DBL_EPSILON + s1)*(DBL_EPSILON + s1));
	double a2 = 6.0/10.0 * 1.0/((DBL_EPSILON + s2)*(DBL_EPSILON + s2));
	double a3 = 3.0/10.0 * 1.0/((DBL_EPSILON + s3)*(DBL_EPSILON + s3));
	double w1 = a1/(a1+a2+a3);
	double w2 = a2/(a1+a2+a3);
	double w3 = a3/(a1+a2+a3);

	return w1*(1.0/3.0*v1-7.0/6.0*v2+11.0/6.0*v3) + w2*(-1.0/6.0*v2+5.0/6.0*v3+1.0/3.0*v4) + w3*(1.0/3.0*v3+5.0/6.0*v4-1.0/6.0*v5);

}

double wenoApproximationXPlus(Phi2D* inputPhi)
{
	double v1; double v2; double v3; double v4; double v5;

	if (inputPhi == NULL || inputPhi->PhiRight == NULL || inputPhi->PhiRight->PhiRight == NULL || inputPhi->PhiRight->PhiRight->PhiRight == NULL || inputPhi->PhiLeft == NULL || inputPhi->PhiLeft->PhiLeft == NULL)
	{
		return dxPlusPhi(inputPhi);
	}

	v1 = (inputPhi->PhiRight->PhiRight->PhiRight->phi - inputPhi->PhiRight->PhiRight->phi) / (inputPhi->PhiRight->PhiRight->PhiRight->x - inputPhi->PhiRight->PhiRight->x);
	v2 = (inputPhi->PhiRight->PhiRight->phi - inputPhi->PhiRight->phi) / (inputPhi->PhiRight->PhiRight->x - inputPhi->PhiRight->x);
	v3 = (inputPhi->PhiRight->phi - inputPhi->phi) / (inputPhi->PhiRight->x - inputPhi->x);
	v4 = (inputPhi->phi - inputPhi->PhiLeft->phi) / (inputPhi->x - inputPhi->PhiLeft->x);
	v5 = (inputPhi->PhiLeft->phi - inputPhi->PhiLeft->PhiLeft->phi) / (inputPhi->PhiLeft->x - inputPhi->PhiLeft->PhiLeft->x);

	double s1 = 13.0/12.0*(v1-2.0*v2+v3)*(v1-2.0*v2+v3) + 1.0/4.0*(v1-4.0*v2+3.0*v3)*(v1-4.0*v2+3.0*v3);
	double s2 = 13.0/12.0*(v2-2.0*v3+v4)*(v2-2.0*v3+v4) + 1.0/4.0*(v2-v4)*(v2-v4);
	double s3 = 13.0/12.0*(v3-2.0*v4+v5)*(v3-2.0*v4+v5) + 1.0/4.0*(3.0*v3-4.0*v4+v5)*(3.0*v3-4.0*v4+v5);

	double eps = 1e-6;

	double a1 = 1.0/10.0 * 1.0/((DBL_EPSILON + s1)*(DBL_EPSILON + s1));
	double a2 = 6.0/10.0 * 1.0/((DBL_EPSILON + s2)*(DBL_EPSILON + s2));
	double a3 = 3.0/10.0 * 1.0/((DBL_EPSILON + s3)*(DBL_EPSILON + s3));
	double w1 = a1/(a1+a2+a3);
	double w2 = a2/(a1+a2+a3);
	double w3 = a3/(a1+a2+a3);

	return w1*(1.0/3.0*v1-7.0/6.0*v2+11.0/6.0*v3) + w2*(-1.0/6.0*v2+5.0/6.0*v3+1.0/3.0*v4) + w3*(1.0/3.0*v3+5.0/6.0*v4-1.0/6.0*v5);

}

double wenoApproximationYMinus(Phi2D* inputPhi)
{
	double v1; double v2; double v3; double v4; double v5;

	if (inputPhi == NULL || inputPhi->PhiTop == NULL || inputPhi->PhiTop->PhiTop == NULL || inputPhi->PhiBottom == NULL || inputPhi->PhiBottom->PhiBottom == NULL || inputPhi->PhiBottom->PhiBottom->PhiBottom == NULL)
	{
		return dyMinusPhi(inputPhi);
	}

	v1 = (inputPhi->PhiBottom->PhiBottom->phi - inputPhi->PhiBottom->PhiBottom->PhiBottom->phi) / (inputPhi->PhiBottom->PhiBottom->y - inputPhi->PhiBottom->PhiBottom->PhiBottom->y);
	v2 = (inputPhi->PhiBottom->phi - inputPhi->PhiBottom->PhiBottom->phi) / (inputPhi->PhiBottom->y - inputPhi->PhiBottom->PhiBottom->y);
	v3 = (inputPhi->phi - inputPhi->PhiBottom->phi) / (inputPhi->y - inputPhi->PhiBottom->y);
	v4 = (inputPhi->PhiTop->phi - inputPhi->phi) / (inputPhi->PhiTop->y - inputPhi->y);
	v5 = (inputPhi->PhiTop->PhiTop->phi - inputPhi->PhiTop->phi) / (inputPhi->PhiTop->PhiTop->y - inputPhi->PhiTop->y);

	double s1 = 13.0/12.0*(v1-2.0*v2+v3)*(v1-2.0*v2+v3) + 1.0/4.0*(v1-4.0*v2+3.0*v3)*(v1-4.0*v2+3.0*v3);
	double s2 = 13.0/12.0*(v2-2.0*v3+v4)*(v2-2.0*v3+v4) + 1.0/4.0*(v2-v4)*(v2-v4);
	double s3 = 13.0/12.0*(v3-2.0*v4+v5)*(v3-2.0*v4+v5) + 1.0/4.0*(3.0*v3-4.0*v4+v5)*(3.0*v3-4.0*v4+v5);

	double eps = 1e-6;

	double a1 = 1.0/10.0 * 1.0/((DBL_EPSILON + s1)*(DBL_EPSILON + s1));
	double a2 = 6.0/10.0 * 1.0/((DBL_EPSILON + s2)*(DBL_EPSILON + s2));
	double a3 = 3.0/10.0 * 1.0/((DBL_EPSILON + s3)*(DBL_EPSILON + s3));
	double w1 = a1/(a1+a2+a3);
	double w2 = a2/(a1+a2+a3);
	double w3 = a3/(a1+a2+a3);

	return w1*(1.0/3.0*v1-7.0/6.0*v2+11.0/6.0*v3) + w2*(-1.0/6.0*v2+5.0/6.0*v3+1.0/3.0*v4) + w3*(1.0/3.0*v3+5.0/6.0*v4-1.0/6.0*v5);

}

double wenoApproximationYPlus(Phi2D* inputPhi)
{
	double v1; double v2; double v3; double v4; double v5;

	if (inputPhi == NULL || inputPhi->PhiTop == NULL || inputPhi->PhiTop->PhiTop == NULL || inputPhi->PhiTop->PhiTop->PhiTop == NULL || inputPhi->PhiBottom == NULL || inputPhi->PhiBottom->PhiBottom == NULL)
	{
		return dyPlusPhi(inputPhi);
	}

	v1 = (inputPhi->PhiTop->PhiTop->PhiTop->phi - inputPhi->PhiTop->PhiTop->phi) / (inputPhi->PhiTop->PhiTop->PhiTop->y - inputPhi->PhiTop->PhiTop->y);
	v2 = (inputPhi->PhiTop->PhiTop->phi - inputPhi->PhiTop->phi) / (inputPhi->PhiTop->PhiTop->y - inputPhi->PhiTop->y);
	v3 = (inputPhi->PhiTop->phi - inputPhi->phi) / (inputPhi->PhiTop->y - inputPhi->y);
	v4 = (inputPhi->phi - inputPhi->PhiBottom->phi) / (inputPhi->y - inputPhi->PhiBottom->y);
	v5 = (inputPhi->PhiBottom->phi - inputPhi->PhiBottom->PhiBottom->phi) / (inputPhi->PhiBottom->y - inputPhi->PhiBottom->PhiBottom->y);

	double s1 = 13.0/12.0*(v1-2.0*v2+v3)*(v1-2.0*v2+v3) + 1.0/4.0*(v1-4.0*v2+3.0*v3)*(v1-4.0*v2+3.0*v3);
	double s2 = 13.0/12.0*(v2-2.0*v3+v4)*(v2-2.0*v3+v4) + 1.0/4.0*(v2-v4)*(v2-v4);
	double s3 = 13.0/12.0*(v3-2.0*v4+v5)*(v3-2.0*v4+v5) + 1.0/4.0*(3.0*v3-4.0*v4+v5)*(3.0*v3-4.0*v4+v5);

	double eps = 1e-6;

	double a1 = 1.0/10.0 * 1.0/((DBL_EPSILON + s1)*(DBL_EPSILON + s1));
	double a2 = 6.0/10.0 * 1.0/((DBL_EPSILON + s2)*(DBL_EPSILON + s2));
	double a3 = 3.0/10.0 * 1.0/((DBL_EPSILON + s3)*(DBL_EPSILON + s3));
	double w1 = a1/(a1+a2+a3);
	double w2 = a2/(a1+a2+a3);
	double w3 = a3/(a1+a2+a3);

	return w1*(1.0/3.0*v1-7.0/6.0*v2+11.0/6.0*v3) + w2*(-1.0/6.0*v2+5.0/6.0*v3+1.0/3.0*v4) + w3*(1.0/3.0*v3+5.0/6.0*v4-1.0/6.0*v5);
}


double reinitialGodonov(double a, double b, double	c, double d, double inputPhi)
{
	if (inputPhi<=0)
	{
		double aPlus  = max(a,0.0);
		double bMinus = min(b,0.0);
		double cPlus  = max(c,0.0);
		double dMinus = min(d,0.0);

		return sqrt(max(aPlus*aPlus, bMinus*bMinus) + max(cPlus*cPlus, dMinus*dMinus))-1.0;
	}
	else
	{
		double aMinus = min(a,0.0);
		double bPlus  = max(b,0.0);
		double cMinus = min(c,0.0);
		double dPlus  = max(d,0.0);

		return sqrt(max(aMinus*aMinus, bPlus*bPlus) + max(cMinus*cMinus, dPlus*dPlus))-1.0;
	}
}

double propagatingGodonov(double a, double b, double	c, double d, double inputPhi)
{
	if (inputPhi<=0)
	{
		double aPlus  = max(a,0.0);
		double bMinus = min(b,0.0);
		double cPlus  = max(c,0.0);
		double dMinus = min(d,0.0);

		return sqrt(max(aPlus*aPlus, bMinus*bMinus) + max(cPlus*cPlus, dMinus*dMinus));
	}
	else
	{
		double aMinus = min(a,0.0);
		double bPlus  = max(b,0.0);
		double cMinus = min(c,0.0);
		double dPlus  = max(d,0.0);

		return sqrt(max(aMinus*aMinus, bPlus*bPlus) + max(cMinus*cMinus, dPlus*dPlus));
	}
}