#include <iostream>
#include <stdio.h>
#include <cmath>
#include <fstream>
#include <cassert>
#include <algorithm>
#include <sstream>
#include <tuple>
#include "reinitialization.h"

#define PI 3.141592653589793238462643383279502884197169399375105820974944592307816406
using namespace std;


int main()
{
	int dimension; // inicate 1D, 2D, 3D.
	int num_mesh, num_xmesh, num_ymesh, num_step;  // Decide the number of mesh , not size. And times.
	double x1, x2, dx, y1,y2, dy, dt, t_total; // 1 is a starting point. 2 is a ending point. and times.
	double side_length, side_xlength, side_ylength; // length foreach axis. 
	// codinate vectors for each axis, and TImes.
	double* X;
	double* Y;
	double* L0;
	double* L1;


	// Level sets.
	double** initial_level_set;
	double** Psi1;
	double** Psi2;
	double** Psi_temp1;
	double** Psi_temp2;
	double** sign_distance_level_set;
	// give initial surface.

	int** closed_coordinate;
	double cfl = 0.45;
	double dt_ij;
	double l0_error;
	double l1_error;
	double error;

	dimension = 2;
	side_length = 4.0;
	num_mesh = 151;
	//num_mesh = 41; 

	dt = 0.0005; // time step.

	t_total = 0.2;
	num_step = ceil(t_total/dt);



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

	L1 = new double [num_step];
	L0 = new double [num_step];
	initial_level_set = new double*[num_xmesh];
	Psi1 = new double* [num_xmesh];
	Psi2 = new double* [num_xmesh];
	Psi_temp2 = new double* [num_xmesh];
	sign_distance_level_set = new double* [num_xmesh];
	Psi_temp1 = new double* [num_xmesh];
	closed_coordinate = new int* [num_xmesh];
	for (int i = 0; i < num_mesh; i++)
	{
		X[i]=x1+dx*((double)i);
		Y[i]=y1+dy*((double)i);		
		initial_level_set[i] = new double [num_ymesh];
		Psi1[i] = new double [num_ymesh];				
		Psi2[i] = new double [num_ymesh];	
		Psi_temp2[i]= new double [num_ymesh];	
		sign_distance_level_set[i]= new double [num_ymesh];	
		Psi_temp1[i]=new double [num_ymesh];
		//closed_coordinate = new double [num_ymesh];

	}



	// Define initial level set which is not a signed distance function

	ofstream initial_levelset_output("initial level set.dat"); // Write data file.
	assert(initial_levelset_output.is_open());
	ofstream distance_levelset_output("sign distance function.dat");
	assert(distance_levelset_output.is_open());
	double a=0.7;
	double r=1.0;
	for (int i = 0; i < num_xmesh; i++) // Loop for X Axis
	{
		for (int j = 0; j < num_ymesh; j++) // Loop for Y axis
		{
			//// A circle with center at the origen and radius 1.
			//initial_level_set[i][j] = (sqrt(X[i]*X[i]+Y[j]*Y[j])-1.0)*((X[i]-1.0)*(X[i]-1.0)+(Y[j]-1.0)*(Y[j]-1.0)+0.1);
			//sign_distance_level_set[i][j] = sqrt(X[i]*X[i]+Y[j]*Y[j])-1.0;


			//// Two circles of radius r are placed at (+-a,0) on the plane. Let 0<a<r, sh that the two circles intersect each other.
			double temp1 = (a-X[i])/sqrt((a-X[i])*(a-X[i])+Y[j]*Y[j]);
			double temp2 = (a+X[i])/sqrt((a+X[i])*(a+X[i])+Y[j]*Y[j]);
			if (temp1>= a/r && temp2 >=a/r)
			{
				double temp3 = min(  X[i]*X[i] + (Y[j]+sqrt(r*r-a*a))*(Y[j]+sqrt(r*r-a*a)),   X[i]*X[i] + (Y[j]-sqrt(r*r-a*a))*(Y[j]-sqrt(r*r-a*a)));
				initial_level_set[i][j] =sqrt(temp3) * ((X[i]-1.0)*(X[i]-1.0)+(Y[j]-1.0)*(Y[j]-1.0)+0.1);
				sign_distance_level_set[i][j] =sqrt(temp3);
			}
			else
			{
				double temp3 = min(  sqrt((X[i]+a)*(X[i]+a) + Y[j]*Y[j]) - r ,   sqrt((X[i]-a)*(X[i]-a) + Y[j]*Y[j]) - r);
				initial_level_set[i][j] =temp3 * ((X[i]-1.0)*(X[i]-1.0)+(Y[j]-1.0)*(Y[j]-1.0)+0.1); 
				sign_distance_level_set[i][j] =temp3 ;
			}


			initial_levelset_output << initial_level_set[i][j] << " ";
			distance_levelset_output << sign_distance_level_set[i][j] << " ";
		}
		initial_levelset_output << " "<<"\n";
		distance_levelset_output << " "<<"\n";
	}
	initial_levelset_output.close();
	distance_levelset_output.close();

	Psi1 = initial_level_set;





	//////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	//                       reinitialization
	//
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	double a1;
	double b;
	double c;
	char string_x='x';
	char string_y='y';


	////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	////
	////                       TVD Runge_Kutta
	////
	////                             START
	////
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////

	for (int t = 0; t < floor(t_total/dt); t++)
	{
		cout<<"TVD Runge-Kutta time step " <<t<<"\n";

		l1_error = 0;
		l0_error = 0;
		error = 0;
		//ostringstream file_name1;
		//file_name1 <<  "D:/TVD RK/reinitialization set temp "<<t<<".dat";
		//ofstream levelset_output1(file_name1.str());
		//assert(levelset_output1.is_open());

		// Update Psi_temp1.
		for (int i = 0; i < num_xmesh; i++)
		{
			for (int j = 0; j < num_ymesh; j++)
			{

				tuple<double,double> dxplus = derivation1x_plus(i,j,Psi1,dx,num_xmesh);
				tuple<double,double> dxminus = derivation1x_minus(i,j,Psi1,dx,num_xmesh);
				tuple<double,double> dyplus = derivation1y_plus(i,j,Psi1,dy,num_ymesh);
				tuple<double,double> dyminus = derivation1y_minus(i,j,Psi1,dy,num_ymesh);

				double derv_xplus, derv_xminus, derv_yplus, derv_yminus;
				double delta1, delta2, delta3, delta4;

				// D+-xy .   i.e. difference value
				derv_xplus = get<0>(dxplus) ;
				derv_xminus = get<0>(dxminus) ;
				derv_yplus = get<0>(dyplus) ;
				derv_yminus = get<0>(dyminus) ;
				// delta x+,delta x-,delta y+,delta y-.
				delta1 = get<1>(dxplus) ;
				delta2 = get<1>(dxminus) ;
				delta3 = get<1>(dyplus) ;
				delta4 = get<1>(dyminus) ;

				dt_ij = cfl*min(min(delta1, delta2),min(delta3,delta4));

				//Psi_temp1[i][j] = Psi1[i][j] - dt_ij*sign_func(Psi1[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,Psi1[i][j])-1.0);
				Psi_temp1[i][j] = Psi1[i][j] - dt_ij*sign_func(Psi1[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,Psi1[i][j])-1.0);
				a = Psi_temp1[i][j];
				a1 = Psi1[i][j];
				b = 0.0;
				c = HG(derv_xplus,derv_xminus, derv_yplus, derv_yminus, Psi1[i][j])-1.0;
				//levelset_output1 << Psi_temp1[i][j] << " ";
			}

			//levelset_output1 << " "<<"\n";
		}
		//levelset_output1.close();
		// Update Psi_temp1 -> Psi1.
		//Psi1=Psi_temp1;



		//////////////////////////////////////////////////////////////////////////////////////////////////
		//
		//                     Psi2 is updated from Psi_temp1
		//
		///////////////////////////////////////////////////////////////////////////////////////////////////
		ostringstream file_name2;
		file_name2 <<  "D:/TVD RK/reinitialization set "<<t<<".dat";
		ofstream levelset_output2(file_name2.str());
		assert(levelset_output2.is_open());

		ostringstream file_name3;
		file_name3 <<  "D:/TVD RK/Error "<<t<<".dat";
		ofstream levelset_output3(file_name3.str());
		assert(levelset_output3.is_open());

		// Update Psi_temp2
		for (int i = 0; i < num_xmesh; i++)
		{
			for (int j = 0; j < num_ymesh; j++)
			{

				tuple<double,double> dxplus = derivation1x_plus(i,j,Psi_temp1,dx,num_xmesh);
				tuple<double,double> dxminus = derivation1x_minus(i,j,Psi_temp1,dx,num_xmesh);
				tuple<double,double> dyplus = derivation1y_plus(i,j,Psi_temp1,dy,num_ymesh);
				tuple<double,double> dyminus = derivation1y_minus(i,j,Psi_temp1,dy,num_ymesh);
				double derv_xplus,derv_xminus,derv_yplus,derv_yminus;
				double delta1,delta2,delta3,delta4;

				// D+-xy .   i.e. difference value
				derv_xplus = get<0>(dxplus) ;
				derv_xminus = get<0>(dxminus) ;
				derv_yplus = get<0>(dyplus) ;
				derv_yminus = get<0>(dyminus) ;
				// delta x+,delta x-,delta y+,delta y-.
				delta1 = get<1>(dxplus) ;
				delta2 = get<1>(dxminus) ;
				delta3 = get<1>(dyplus) ;
				delta4 = get<1>(dyminus) ;

				dt_ij = cfl * min(min(delta1, delta2),min(delta3,delta4));


				Psi_temp2[i][j] = Psi_temp1[i][j] - dt_ij*sign_func(Psi_temp1[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,Psi_temp1[i][j])-1.0);
				double a=Psi_temp2[i][j];
				double b=0;

				Psi2[i][j]= (Psi1[i][j]+Psi_temp2[i][j])/2.0;
				a=Psi2[i][j];
				b=0;

				levelset_output2 << Psi2[i][j] << " ";
				error =sign_distance_level_set[i][j] - Psi2[i][j];
				levelset_output3 << error << " ";

				// Compute L1 and L2 error.
				if (l0_error<abs(error))
					{
						l0_error = abs(error);
					}
				l1_error += abs(error);

			} //end 'for' with j

			levelset_output2 << " "<<"\n";
			levelset_output3 << " "<<"\n";
		} // emd 'for' with i
		levelset_output2.close();
		levelset_output3.close();


		L1[t] = l1_error*(dx*dy);
		L0[t] = l0_error;

		ostringstream file_name4;
		file_name4 <<  "D:/TVD RK/L0 "<<t<<".dat";
		ofstream L0_output4(file_name4.str());
		assert(L0_output4.is_open());
		L0_output4 << L0[t] << " ";
		L0_output4.close();

		ostringstream file_name5;
		file_name5 <<  "D:/TVD RK/L1 "<<t<<".dat";
		ofstream L1_output(file_name5.str());
		assert(L1_output.is_open());
		L1_output << L1[t] << " ";
		L1_output.close();

		// Update Psi2 -> Psi1.
		Psi1 = Psi2;
	} // end 'for' with t
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	////
	////                       TVD Runge_Kutta
	////
	////                             END
	////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////



	////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	////
	////                       Forward Euler method
	////
	////                             START
	////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	Psi1 = initial_level_set;
	for (int t = 0; t < floor(t_total/dt); t++)
	{
		cout<<"Forward Euler time step " <<t<<"\n";

		l1_error = 0;
		l0_error = 0;
		error = 0;
		//////////////////////////////////////////////////////////////////////////////////////////////////
		//
		//                     psi2 is updated from psi1
		//
		///////////////////////////////////////////////////////////////////////////////////////////////////
		ostringstream file_name1;
		file_name1 <<  "D:/FE/reinitialization set "<<t<<".dat";
		ofstream levelset_output1(file_name1.str());
		assert(levelset_output1.is_open());

		ostringstream file_name3;
		file_name3 <<  "D:/FE/Error "<<t<<".dat";
		ofstream levelset_output3(file_name3.str());
		assert(levelset_output3.is_open());

		// update psi2.
		for (int i = 0; i < num_xmesh; i++)
		{
			for (int j = 0; j < num_ymesh; j++)
			{

				tuple<double,double> dxplus = derivation1x_plus(i,j,Psi1,dx,num_xmesh);
				tuple<double,double> dxminus = derivation1x_minus(i,j,Psi1,dx,num_xmesh);
				tuple<double,double> dyplus = derivation1y_plus(i,j,Psi1,dy,num_ymesh);
				tuple<double,double> dyminus = derivation1y_minus(i,j,Psi1,dy,num_ymesh);

				double derv_xplus, derv_xminus, derv_yplus, derv_yminus;
				double delta1, delta2, delta3, delta4;

				// d+-xy .   i.e. difference value
				derv_xplus = get<0>(dxplus) ;
				derv_xminus = get<0>(dxminus) ;
				derv_yplus = get<0>(dyplus) ;
				derv_yminus = get<0>(dyminus) ;
				// delta x+,delta x-,delta y+,delta y-.
				delta1 = get<1>(dxplus) ;
				delta2 = get<1>(dxminus) ;
				delta3 = get<1>(dyplus) ;
				delta4 = get<1>(dyminus) ;

				dt_ij = cfl*min(min(delta1, delta2),min(delta3,delta4));

				//psi_temp1[i][j] = psi1[i][j] - dt_ij*sign_func(psi1[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,psi1[i][j])-1.0);
				Psi2[i][j] = Psi1[i][j] - dt_ij*sign_func(Psi1[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,Psi1[i][j])-1.0);
				a = Psi2[i][j];
				a1 = Psi1[i][j];
				b = 0.0;
				c = HG(derv_xplus,derv_xminus, derv_yplus, derv_yminus, Psi1[i][j])-1.0;
				levelset_output1 << Psi2[i][j] << " ";
				error = sign_distance_level_set[i][j] - Psi2[i][j];
				levelset_output3 << error << " ";

				// Compute L1 and L2 error.
				if (l0_error<abs(error))
					{
						l0_error = abs(error);
					}
				l1_error += abs(error);
			}

			levelset_output1 << " "<<"\n";
			levelset_output3 << " "<<"\n";
		}
		levelset_output1.close();
		levelset_output3.close();


		L1[t] = l1_error*(dx*dy);
		L0[t] = l0_error;

		ostringstream file_name4;
		file_name4 <<  "D:/FE/L0 "<<t<<".dat";
		ofstream L0_output4(file_name4.str());
		assert(L0_output4.is_open());
		L0_output4 << L0[t] << " ";
		L0_output4.close();

		ostringstream file_name5;
		file_name5 <<  "D:/FE/L1 "<<t<<".dat";
		ofstream L1_output(file_name5.str());
		assert(L1_output.is_open());
		L1_output << L1[t] << " ";
		L1_output.close();


		// update psi_temp1 -> psi1.
		Psi1 = Psi2;
	}
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	////
	////                       Forward Euler method
	////
	////                             END
	////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////



	////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	////
	////                       Gauss-Seidel iteration
	////
	////                             START
	////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	Psi1 = initial_level_set;
	for (int t = 0; t < floor(t_total/dt); t++)
	{
		cout<<"Gauss-Seidel iteration time step " <<t<<"\n";

		l1_error = 0;
		l0_error = 0;
		error = 0;
		int raster_visiting;
		raster_visiting = t%4;


		//////////////////////////////////////////////////////////////////////////////////////////////////
		//
		//                     Psi1 is updated from Psi1
		//
		///////////////////////////////////////////////////////////////////////////////////////////////////
		ostringstream file_name1;
		file_name1 <<  "D:/Gauss Seidel/reinitialization set "<<t<<".dat";
		ofstream levelset_output1(file_name1.str());
		assert(levelset_output1.is_open());

		ostringstream file_name3;
		file_name3 <<  "D:/Gauss Seidel/Error "<<t<<".dat";
		ofstream levelset_output3(file_name3.str());
		assert(levelset_output3.is_open());

		// Update Psi1, itself.
		switch (raster_visiting)
		{

		case 0 : 
			for (int i = 0; i < num_xmesh; i++)
			{
				for (int j = 0; j < num_ymesh; j++)
				{
					tuple<double,double> dxplus = derivation1x_plus(i,j,Psi1,dx,num_xmesh);
					tuple<double,double> dxminus = derivation1x_minus(i,j,Psi1,dx,num_xmesh);
					tuple<double,double> dyplus = derivation1y_plus(i,j,Psi1,dy,num_ymesh);
					tuple<double,double> dyminus = derivation1y_minus(i,j,Psi1,dy,num_ymesh);

					double derv_xplus, derv_xminus, derv_yplus, derv_yminus;
					double delta1, delta2, delta3, delta4;

					// D+-xy .   i.e. difference value
					derv_xplus = get<0>(dxplus) ;
					derv_xminus = get<0>(dxminus) ;
					derv_yplus = get<0>(dyplus) ;
					derv_yminus = get<0>(dyminus) ;
					// delta x+,delta x-,delta y+,delta y-.
					delta1 = get<1>(dxplus) ;
					delta2 = get<1>(dxminus) ;
					delta3 = get<1>(dyplus) ;
					delta4 = get<1>(dyminus) ;

					dt_ij = cfl*min(min(delta1, delta2),min(delta3,delta4));

					//Psi_temp1[i][j] = Psi1[i][j] - dt_ij*sign_func(Psi1[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,Psi1[i][j])-1.0);
					Psi1[i][j] = Psi1[i][j] - dt_ij*sign_func(Psi1[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,Psi1[i][j])-1.0);
					a1 = Psi1[i][j];
					b = 0.0;
					c = HG(derv_xplus,derv_xminus, derv_yplus, derv_yminus, Psi1[i][j])-1.0;
					levelset_output1 << Psi1[i][j] << " ";
					error = sign_distance_level_set[i][j] - Psi1[i][j];
					levelset_output3 << error << " ";

					// Compute L1 and L2 error.
					if (l0_error<abs(error))
					{
						l0_error = abs(error);
					}
					l1_error += abs(error);

				} //end 'for' with j

				levelset_output1 << " "<<"\n";
				levelset_output3 << " "<<"\n";
			} //end 'for' with i
			break;




		case 1 : 
			for (int i = 0; i < num_xmesh; i++)
			{
				for (int j = num_ymesh-1; j > -1; j--)
				{
					tuple<double,double> dxplus = derivation1x_plus(i,j,Psi1,dx,num_xmesh);
					tuple<double,double> dxminus = derivation1x_minus(i,j,Psi1,dx,num_xmesh);
					tuple<double,double> dyplus = derivation1y_plus(i,j,Psi1,dy,num_ymesh);
					tuple<double,double> dyminus = derivation1y_minus(i,j,Psi1,dy,num_ymesh);

					double derv_xplus, derv_xminus, derv_yplus, derv_yminus;
					double delta1, delta2, delta3, delta4;

					// D+-xy .   i.e. difference value
					derv_xplus = get<0>(dxplus) ;
					derv_xminus = get<0>(dxminus) ;
					derv_yplus = get<0>(dyplus) ;
					derv_yminus = get<0>(dyminus) ;
					// delta x+,delta x-,delta y+,delta y-.
					delta1 = get<1>(dxplus) ;
					delta2 = get<1>(dxminus) ;
					delta3 = get<1>(dyplus) ;
					delta4 = get<1>(dyminus) ;

					dt_ij = cfl*min(min(delta1, delta2),min(delta3,delta4));

					//Psi_temp1[i][j] = Psi1[i][j] - dt_ij*sign_func(Psi1[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,Psi1[i][j])-1.0);
					Psi1[i][j] = Psi1[i][j] - dt_ij*sign_func(Psi1[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,Psi1[i][j])-1.0);
					a1 = Psi1[i][j];
					b = 0.0;
					c = HG(derv_xplus,derv_xminus, derv_yplus, derv_yminus, Psi1[i][j])-1.0;
					/*levelset_output1 << Psi1[i][j] << " ";
					levelset_output3 << sign_distance_level_set[i][j] - Psi1[i][j] << " ";*/
				} //end 'for' with j

				for (int j = 0; j  < num_ymesh; j ++)
				{
					levelset_output1 << Psi1[i][j] << " ";
					error = sign_distance_level_set[i][j] - Psi1[i][j];
					levelset_output3 << error << " ";

					// Compute L1 and L2 error.
					if (l0_error<abs(error))
					{
						l0_error = abs(error);
					}
					l1_error += abs(error);
				}

				levelset_output1 << " "<<"\n";
				levelset_output3 << " "<<"\n";
			} //end 'for' with i
			break;




		case 2 : 
			for (int i = num_xmesh-1; i>-1 ; i--)
			{
				for (int j = 0; j < num_ymesh; j++)
				{
					tuple<double,double> dxplus = derivation1x_plus(i,j,Psi1,dx,num_xmesh);
					tuple<double,double> dxminus = derivation1x_minus(i,j,Psi1,dx,num_xmesh);
					tuple<double,double> dyplus = derivation1y_plus(i,j,Psi1,dy,num_ymesh);
					tuple<double,double> dyminus = derivation1y_minus(i,j,Psi1,dy,num_ymesh);

					double derv_xplus, derv_xminus, derv_yplus, derv_yminus;
					double delta1, delta2, delta3, delta4;

					// D+-xy .   i.e. difference value
					derv_xplus = get<0>(dxplus) ;
					derv_xminus = get<0>(dxminus) ;
					derv_yplus = get<0>(dyplus) ;
					derv_yminus = get<0>(dyminus) ;
					// delta x+,delta x-,delta y+,delta y-.
					delta1 = get<1>(dxplus) ;
					delta2 = get<1>(dxminus) ;
					delta3 = get<1>(dyplus) ;
					delta4 = get<1>(dyminus) ;

					dt_ij = cfl*min(min(delta1, delta2),min(delta3,delta4));

					//Psi_temp1[i][j] = Psi1[i][j] - dt_ij*sign_func(Psi1[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,Psi1[i][j])-1.0);
					Psi1[i][j] = Psi1[i][j] - dt_ij*sign_func(Psi1[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,Psi1[i][j])-1.0);
					a1 = Psi1[i][j];
					b = 0.0;
					c = HG(derv_xplus,derv_xminus, derv_yplus, derv_yminus, Psi1[i][j])-1.0;
					//levelset_output1 << Psi1[i][j] << " ";
					//levelset_output3 << sign_distance_level_set[i][j] - Psi1[i][j] << " ";
				} //end 'for' with j
				//levelset_output1 << " "<<"\n";
				//levelset_output3 << " "<<"\n";
			} //end 'for' with i

			for (int i = 0; i < num_xmesh; i++)
			{
				for (int j = 0; j < num_ymesh; j++)
				{
					levelset_output1 << Psi1[i][j] << " ";
					error = sign_distance_level_set[i][j] - Psi1[i][j];
					levelset_output3 << error << " ";

					// Compute L1 and L2 error.
					if (l0_error<abs(error))
					{
						l0_error = abs(error);
					}
					l1_error += abs(error);
				}
				levelset_output1 << " "<<"\n";
				levelset_output3 << " "<<"\n";
			}
			break;



		case 3 : 
			for (int i = num_xmesh-1; i>-1 ; i--)
			{
				for (int j = num_ymesh-1; j > -1; j--)
				{

					tuple<double,double> dxplus = derivation1x_plus(i,j,Psi1,dx,num_xmesh);
					tuple<double,double> dxminus = derivation1x_minus(i,j,Psi1,dx,num_xmesh);
					tuple<double,double> dyplus = derivation1y_plus(i,j,Psi1,dy,num_ymesh);
					tuple<double,double> dyminus = derivation1y_minus(i,j,Psi1,dy,num_ymesh);

					double derv_xplus, derv_xminus, derv_yplus, derv_yminus;
					double delta1, delta2, delta3, delta4;

					// D+-xy .   i.e. difference value
					derv_xplus = get<0>(dxplus) ;
					derv_xminus = get<0>(dxminus) ;
					derv_yplus = get<0>(dyplus) ;
					derv_yminus = get<0>(dyminus) ;
					// delta x+,delta x-,delta y+,delta y-.
					delta1 = get<1>(dxplus) ;
					delta2 = get<1>(dxminus) ;
					delta3 = get<1>(dyplus) ;
					delta4 = get<1>(dyminus) ;

					dt_ij = cfl*min(min(delta1, delta2),min(delta3,delta4));

					//Psi_temp1[i][j] = Psi1[i][j] - dt_ij*sign_func(Psi1[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,Psi1[i][j])-1.0);
					Psi1[i][j] = Psi1[i][j] - dt_ij*sign_func(Psi1[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,Psi1[i][j])-1.0);
					a1 = Psi1[i][j];
					b = 0.0;
					c = HG(derv_xplus,derv_xminus, derv_yplus, derv_yminus, Psi1[i][j])-1.0;
					//levelset_output1 << Psi1[i][j] << " ";
					//levelset_output3 << sign_distance_level_set[i][j] - Psi1[i][j] << " ";
				} //end 'for' with j
				//levelset_output1 << " "<<"\n";
				//levelset_output3 << " "<<"\n";
			} //end 'for' with i

			for (int i = 0; i < num_xmesh; i++)
			{
				for (int j = 0; j < num_ymesh; j++)
				{
					levelset_output1 << Psi1[i][j] << " ";
					error = sign_distance_level_set[i][j] - Psi1[i][j];
					levelset_output3 << error << " ";

					// Compute L1 and L2 error.
					if (l0_error<abs(error))
					{
						l0_error = abs(error);
					}
					l1_error += abs(error);
				}
				levelset_output1 << " "<<"\n";
				levelset_output3 << " "<<"\n";
			}
			break;

		} // end 'switch'

		levelset_output1.close();


		L1[t] = l1_error*(dx*dy);
		L0[t] = l0_error;

		ostringstream file_name4;
		file_name4 <<  "D:/Gauss Seidel/L0 "<<t<<".dat";
		ofstream L0_output4(file_name4.str());
		assert(L0_output4.is_open());
		L0_output4 << L0[t] << " ";
		L0_output4.close();

		ostringstream file_name5;
		file_name5 <<  "D:/Gauss Seidel/L1 "<<t<<".dat";
		ofstream L1_output(file_name5.str());
		assert(L1_output.is_open());
		L1_output << L1[t] << " ";
		L1_output.close();

	}//end 'for' with t

	////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	////
	////                       Gauss-Seidel iteration
	////
	////                             END
	////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////

}