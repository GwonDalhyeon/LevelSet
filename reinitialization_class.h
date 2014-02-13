#include <iostream>
#include <stdio.h>
#include <cmath>
#include <fstream>
#include <cassert>
#include <algorithm>
#include <sstream>
#include <tuple>

#define PI 3.141592653589793238462643383279502884197169399375105820974944592307816406
using namespace std;



class LEVELSET
{
public:
	int dimension; // inicate 1D, 2D, 3D.
	int num_mesh, num_xmesh, num_ymesh, num_step;  // Decide the number of mesh , not size. And times.
	double x1, x2, dx, y1,y2, dy; // 1 is a starting point. 2 is a ending point. and times.
	double side_xlength, side_ylength; // length foreach axis. 
	// codinate vectors for each axis, and TImes.
	double* X;
	double* Y;

	// Level sets.
	double** Psi;
	//double** zero_level_set;

	LEVELSET(double xlength, double ylength, int num_x, int num_y);
	LEVELSET(double xlength, double ylength, double delta_x, double delta_y);
	LEVELSET(const LEVELSET& original_level_set);
	LEVELSET operator=(LEVELSET orignal_LEVELSET)
	{
		dimension = orignal_LEVELSET.dimension; // inicate 1D, 2D, 3D.
		num_mesh = orignal_LEVELSET.num_mesh;
		num_xmesh =orignal_LEVELSET.num_xmesh;
		num_ymesh = orignal_LEVELSET.num_ymesh;
		num_step = orignal_LEVELSET.num_step;  // Decide the number of mesh , not size. And times.

		x1 = orignal_LEVELSET.x1;
		x2 = orignal_LEVELSET.x2;
		dx = orignal_LEVELSET.dx;
		y1 = orignal_LEVELSET.y1;
		y2 = orignal_LEVELSET.y2;
		dy = orignal_LEVELSET.dy; // 1 is a starting point. 2 is a ending point. and times.

		side_xlength = orignal_LEVELSET.side_xlength;
		side_ylength = orignal_LEVELSET.side_ylength; // length foreach axis. 

		X = new double [num_xmesh];
		Y = new double [num_ymesh];

		Psi = new double*[num_xmesh];

		for (int i = 0; i < num_xmesh; i++)
		{
			X[i]=orignal_LEVELSET.X[i];
			Y[i]=orignal_LEVELSET.Y[i];		
			Psi[i] = new double [num_ymesh];
			for (int j = 0; j < num_ymesh; j++)
			{
				Psi[i][j] = orignal_LEVELSET.Psi[i][j];
			}
		}
		return *this;
	}
	~LEVELSET();
private:

};

LEVELSET::LEVELSET(double xlength, double ylength, int num_x, int num_y)
{
	side_xlength = xlength;
	side_ylength = ylength;
	num_xmesh = num_x;
	num_ymesh = num_y;
	x1 = -side_xlength/2; x2 = side_xlength/2;
	y1 = -side_ylength/2; y2 = side_ylength/2;
	dx = ((double)side_xlength)/num_xmesh;
	dy = ((double)side_ylength)/num_ymesh;
	X = new double [num_xmesh];
	Y = new double [num_ymesh];

	Psi = new double*[num_xmesh];

	for (int i = 0; i < num_xmesh; i++)
	{
		X[i]=x1+dx*((double)i);
		Y[i]=y1+dy*((double)i);		
		Psi[i] = new double [num_ymesh];
	}
}

LEVELSET::LEVELSET(double xlength, double ylength, double delta_x, double delta_y)
{
	side_xlength = xlength;
	side_ylength = ylength;
	dx = delta_x;//((double)side_xlength)/num_xmesh;
	dy = delta_y;//((double)side_ylength)/num_ymesh;

	num_xmesh = ceil( side_xlength/dx);
	num_ymesh = ceil( side_ylength/dy);

	x1 = -side_xlength/2; x2 = side_xlength/2;
	y1 = -side_ylength/2; y2 = side_ylength/2;

	X = new double [num_xmesh];
	Y = new double [num_ymesh];

	Psi = new double*[num_xmesh];

	for (int i = 0; i < num_xmesh; i++)
	{
		X[i]=x1+dx*((double)i);
		Y[i]=y1+dy*((double)i);		
		Psi[i] = new double [num_ymesh];
	}
}


LEVELSET :: LEVELSET(const LEVELSET& original_level_set)
{
	dimension = original_level_set.dimension; // inicate 1D, 2D, 3D.
	num_mesh = original_level_set.num_mesh;
	num_xmesh =original_level_set.num_xmesh;
	num_ymesh = original_level_set.num_ymesh;
	num_step = original_level_set.num_step;  // Decide the number of mesh , not size. And times.

	x1 = original_level_set.x1;
	x2 = original_level_set.x2;
	dx = original_level_set.dx;
	y1 = original_level_set.y1;
	y2 = original_level_set.y2;
	dy = original_level_set.dy; // 1 is a starting point. 2 is a ending point. and times.

	side_xlength = original_level_set.side_xlength;
	side_ylength = original_level_set.side_ylength; // length foreach axis. 

	X = new double [num_xmesh];
	Y = new double [num_ymesh];

	Psi = new double*[num_xmesh];

	for (int i = 0; i < num_xmesh; i++)
	{
		X[i]=original_level_set.X[i];
		Y[i]=original_level_set.Y[i];		
		Psi[i] = new double [num_ymesh];
		for (int j = 0; j < num_ymesh; j++)
		{
			Psi[i][j] = original_level_set.Psi[i][j];
		}
	}

}

LEVELSET::~LEVELSET()
{
	delete []X;
	delete []Y;

	for (int i = 0; i < num_xmesh; i++)
	{
		delete [] Psi[i];
		//delete [] zero_level_set[i];
	}

}



class point_set
{
public:
	double* pointX;
	double* pointY;
	//int num_point;
	point_set();
	point_set(LEVELSET input_level_set, int num_xmax, int num_ymax);
	~point_set();

private:

};




class space_domain_info
{

public:
	double side_xlength;
	double side_ylength;
	int num_xmesh;
	int num_ymesh;
	double dx;
	double dy;

	space_domain_info(double side_xlength, double side_ylength, int num_xmesh, int num_ymesh);
	space_domain_info(double xlength, double ylength, double delta_x, double delta_y);
	~space_domain_info();

private:

};

space_domain_info::space_domain_info(double xlength, double ylength, int num_x, int num_y)
{
	double side_xlength = xlength;
	double side_ylength = ylength;
	int num_xmesh = num_x;
	int num_ymesh = num_y;
}

space_domain_info::space_domain_info(double xlength, double ylength, double delta_x, double delta_y)
{
	double side_xlength = xlength;
	double side_ylength = ylength;
	double dx = delta_x;
	double dy = delta_y;
}

space_domain_info::~space_domain_info()
{
}
