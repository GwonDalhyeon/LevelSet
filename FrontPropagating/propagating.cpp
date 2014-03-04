#include <iostream>
#include <stdio.h>
#include <cmath>
#include <fstream>
#include <cassert>
#include <algorithm>
#include <sstream>
#include "propagating.h"
#define PI 3.141592653589793238462643383279502884197169399375105820974944592307816406
using namespace std;

int main()
{
	int dimension; // inicate 1D, 2D, 3D.
	int num_mesh, num_xmesh, num_ymesh, num_zmesh, num_step;  // Decide the number of mesh , not size. And times.
	double x1, x2, dx, y1,y2, dy, z1,z2, dz, dt, t_total; // 1 is a starting point. 2 is a ending point. and times.
	double side_length, side_xlength, side_ylength, side_zlength; // length foreach axis. 
	 // codinate vectors for each axis, and TImes.
	double* X;
	double* Y;
	double* Z;
	double* T;
	 // Level sets.
	double** initial_level_set;
	double** Psi1;
	double** Psi2;
	double** Psi_temp;
	 // give initial surface.
	double** initial_surface;
	
	dimension = 2;
	side_length = 1;
	num_mesh = 601; // use 600 mesh points per side.
	dt = 0.0005; // time step.
	t_total = 2.0;
	num_step = ceil(t_total/dt);
	T = new double [num_step];



	 // The computational domain is a square centered at the origin of side length 1.
	/*if (dimension == 2)
	{*/
		num_xmesh = num_mesh;
		num_ymesh = num_mesh;
		side_xlength = side_length;
		side_ylength = side_length;
		x1 = -side_xlength/2; x2 = side_xlength/2;
		y1 = -side_ylength/2; y2 = side_ylength/2;
		dx = ((double)side_length)/num_xmesh;
		dy = ((double)side_length)/num_ymesh;
		X = new double [num_xmesh];
		Y = new double [num_ymesh];

		initial_level_set = new double*[num_xmesh];
		Psi1 = new double* [num_xmesh];
		Psi2 = new double* [num_xmesh];
		Psi_temp = new double* [num_xmesh];
		for (int i = 0; i < num_mesh; i++)
		{
			X[i]=x1+dx*((double)i);
			Y[i]=y1+dy*((double)i);		
			initial_level_set[i] = new double [num_ymesh];
			Psi1[i] = new double [num_ymesh];				
			Psi2[i] = new double [num_ymesh];				
			Psi_temp[i]=new double [num_ymesh];
			
		}

		
	/*}
	else if (dimension == 3)
	{
		num_xmesh = num_mesh;
		num_ymesh = num_mesh;
		num_zmesh = num_mesh;
		side_xlength = side_length;
		side_ylength = side_length;
		side_zlength = side_length;
		x1 = -side_xlength/2; x2 = side_xlength/2;
		y1 = -side_ylength/2; y2 = side_ylength/2;
		z1 = -side_zlength/2; z2 = side_zlength/2;
		dx=((double)side_length)/num_xmesh;
		dy=((double)side_length)/num_ymesh;
		dz=((double)side_length)/num_zmesh;
		X = new double [num_xmesh];		
		Y = new double [num_ymesh];
		Z = new double [num_zmesh];

		for (int i = 0; i < num_mesh; i++)
		{
			X[i]=x1+dx*((double)i);
			Y[i]=y1+dy*((double)i);
			Z[i]=z1+dz*((double)i);
		}
	}*/
	



	// Define initial surface.

	double* S;
	int num_initial_point; // the number of point on initial surface.
	num_initial_point = 1000;

	S = new double [num_initial_point];
	for (int i = 0; i < num_initial_point; i++)
	{
		S[i] = 1.0/num_initial_point *i; 
	}
	
	ofstream initial_surface_output("initial surface.dat"); // Write data file.
	assert(initial_surface_output.is_open());
	initial_surface = new double* [num_initial_point]; // contain initial surface coordinate.
	for (int i = 0; i < num_initial_point; i++)
	{
		//seven point star : start
		initial_surface [i]= new double [dimension];
		initial_surface [i][0] =  0.065*sin(7*2*PI*S[i])*cos(2*PI*S[i]);
		initial_surface [i][1] =  0.065*sin(7*2*PI*S[i])*sin(2*PI*S[i]);
		//seven point star : end

		initial_surface_output << initial_surface [i][0] << " " << initial_surface [i][1] << "\n";
	}
	delete[] S;
	initial_surface_output.close();









	// Define initial level set using sign distance function

	ofstream initial_levelset_output("initial level set.dat"); // Write data file.
	assert(initial_levelset_output.is_open());
	double closed_coordinate[1][2], temp_coordinate[1][2];

	for (int i = 0; i < num_xmesh; i++) // Loop for X Axis
	{
		
		for (int j = 0; j < num_ymesh; j++) // Loop for Y axis
		{
			double shortest_distance = sqrt((X[i]-initial_surface[0][0])*(X[i]-initial_surface[0][0])+(Y[j]-initial_surface[0][1])*(Y[j]-initial_surface[0][1]));

			for (int k = 1; k < num_initial_point; k++)		// Loop for initla surface number
			{			
				double temp_distance = 0;
				
				temp_distance = sqrt((X[i]-initial_surface[k][0])*(X[i]-initial_surface[k][0])+(Y[j]-initial_surface[k][1])*(Y[j]-initial_surface[k][1]));
			
				if (temp_distance <= shortest_distance)			
				{
					shortest_distance = temp_distance;	// update shortest distance
					closed_coordinate[0][0]=i;
					closed_coordinate[0][1]=j;
				}
			}
			// Give positive sign on inside initial contour, and negative sign on outside
			if ((X[i])*(X[i])+(Y[j])*(Y[j])>0.065*0.065*sin(7*atan2(Y[j],X[i]))*sin(7*atan2(Y[j],X[i])))
			{
				initial_level_set[i][j] = -shortest_distance;			
				initial_levelset_output << -shortest_distance << " ";
			} 
			else if ((0<=atan2(Y[j],X[i]) && atan2(Y[j],X[i])<PI/7.0) || (2*PI/7<=atan2(Y[j],X[i])&& atan2(Y[j],X[i])<3*PI/7.0) || (4*PI/7<=atan2(Y[j],X[i])&& atan2(Y[j],X[i])<5*PI/7.0) || (6*PI/7<=atan2(Y[j],X[i])&& atan2(Y[j],X[i])<PI) || (-2*PI/7.0<=atan2(Y[j],X[i])&& atan2(Y[j],X[i])<-1*PI/7.0) || (-4*PI/7.0<=atan2(Y[j],X[i])&& atan2(Y[j],X[i])<-3*PI/7.0) || (-6*PI/7.0<=atan2(Y[j],X[i])&& atan2(Y[j],X[i])<-5*PI/7.0))
			{					
				initial_level_set[i][j] = shortest_distance;								
				initial_levelset_output << shortest_distance << " ";				
			}				
			else				
			{			
				initial_level_set[i][j] = -shortest_distance;						
				initial_levelset_output << -shortest_distance << " ";
			}
		}
		initial_levelset_output << " "<<"\n";
	}
	initial_levelset_output.close();


	Psi1 = initial_level_set;

	// Define Psi1.
	/*for (int i = 0; i < num_xmesh; i++)
	{
		for (int j = 0; j < num_ymesh; j++)
		{
			Psi1;
		}
	}*/



	// Propgating part

	for (int t = 0; t  < floor(t_total/dt); t ++)
	{
		cout<<"time step "<< t <<"\n";

		ostringstream file_name;
		file_name <<  "D:\level set "<<t<<".dat";
		
			ofstream levelset_output(file_name.str());
			assert(levelset_output.is_open());
		
		

		for (int i = 0; i < num_xmesh; i++)
		{
			for (int j = 0; j < num_ymesh; j++)
			{
				double dxminus=0,dxplus=0, dyminus=0, dyplus=0;
				double gHJ=0, temp_gHJ=0;

				if (i>0 && i<num_xmesh-1)
				{
					dxminus = Psi1[i][j]-Psi1[i-1][j];
					dxplus = Psi1[i+1][j]-Psi1[i][j];
				}
				else if (i==0)
				{
					dxminus = 0;
					dxplus = Psi1[i+1][j]-Psi1[i][j];
				}
				else
				{
					dxminus = Psi1[i][j]-Psi1[i-1][j];
					dxplus = 0;
				}
				
				if (j>0 && j<num_ymesh-1)
				{
					dyminus = Psi1[i][j]-Psi1[i][j-1];
					dyplus = Psi1[i][j+1]-Psi1[i][j];
				}
				else if (j==0)
				{
					dyminus = 0;
					dyplus = Psi1[i][j+1]-Psi1[i][j];
				}
				else
				{
					dyminus = Psi1[i][j]-Psi1[i][j-1];
					dyplus = 0;
				}
				

				temp_gHJ = min(dxminus,0.0)*min(dxminus,0.0) + max(dxplus,0.0)*max(dxplus,0.0) + min(dyminus,0.0)*min(dyminus,0.0) + max(dyplus,0.0)*max(dyplus,0.0);
				gHJ = - sqrt(temp_gHJ);
				Psi2[i][j] = Psi1[i][j] - dt*gHJ;
				
				levelset_output << Psi2[i][j] << " ";
				
			}
			levelset_output << " "<<"\n";
		}
		levelset_output.close();

		// Update Psi2 -> Psi1.
		Psi1 = Psi2;
	}

}