#include <iostream>
#include <stdio.h>
#include <cmath>
#include <fstream>
#include <cassert>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include "reinitialization_class.h"



///////////////////////////////////////////////////////////////////////////////////
//
//			sign function
//
//////////////////////////////////////////////////////////////////////////////////
double sign_func(double a)
{
	if (abs(a)<DBL_EPSILON)
	{
		return a=0;
	}
	else
	{
		return a/abs(a);
	}
}


///////////////////////////////////////////////////////////////////////////////////
//
//			MINMOD function
//
//////////////////////////////////////////////////////////////////////////////////
double minmod(double a, double b)
{
	if (a*b<0)
	{
		return 0;
	}
	else if (abs(a)>=abs(b))
	{
		return b;
	}
	else
	{
		return a;
	}
}

///////////////////////////////////////////////////////////////////////////////////
//
//			(a)+, (a)- function
//
//////////////////////////////////////////////////////////////////////////////////
double plus_function(double a)
{
	return max(a,0.0);
}


double minus_function(double a)
{
	return min(a,0.0);
}



///////////////////////////////////////////////////////////////////////////////////
//
//			Godunov Hamiltonian
//
//////////////////////////////////////////////////////////////////////////////////
double HG(double a, double b, double c, double d, double signal)
{
	double temp;
	double result;

	if (signal>=0)
	{
		temp = max(minus_function(a)*minus_function(a),plus_function(b)*plus_function(b))+max(minus_function(c)*minus_function(c),plus_function(d)*plus_function(d));
		result = sqrt(temp);
	}
	else
	{
		temp = max(plus_function(a)*plus_function(a),minus_function(b)*minus_function(b))+max(plus_function(c)*plus_function(c),minus_function(d)*minus_function(d));
		result = sqrt(temp);
	}

	return result;
}



///////////////////////////////////////////////////////////////////////////////////
//
//			second derivation function
//
//////////////////////////////////////////////////////////////////////////////////
double derivation2x(int i,int j, double** A, double delta , int max_num)
{
	double result;

	if (i>0 && i<(max_num-1))
	{
		result = A[i-1][j]-2.0*A[i][j]+A[i+1][j];
	}
	else if (i==0)
	{
		result = A[i][j]-2.0*A[i+1][j]+A[i+2][j];
	}
	else if (i==(max_num-1))
	{
		result = A[i][j]-2.0*A[i-1][j]+A[i-2][j];
	}
	else
	{
		result = 0;
	}

	result = result/(delta*delta);
	return result;

}


double derivation2y(int i,int j, double** A, double delta , int max_num)
{
	double result;

	if (j>0 && j<(max_num-1))
	{
		result = A[i][j-1]-2.0*A[i][j]+A[i][j+1];
	}
	else if (j==0)
	{
		result = A[i][j]-2.0*A[i][j+1]+A[i][j+2];
	}
	else if (j==(max_num-1))
	{
		result = A[i][j]-2.0*A[i][j-1]+A[i][j-2];
	}
	else
	{
		result = 0;
	}

	result = result/(delta*delta);
	return result;



}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       computing Dx+,Dy+.
///
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// first derivation plus function 
// (difference value, delta plus)
tuple<double,double> derivation1x_plus(int i,int j, double** A, double delta, int max_num)
{
	double result;
	double Aij;

	// calculate D+xPSI_ij
	Aij=A[i][j];
	if (i<max_num-1)
	{
		result = (A[i+1][j] -A[i][j])/delta - delta/2.0*minmod(derivation2x(i,j,A,delta, max_num),derivation2x(i+1,j,A,delta,max_num));
	}
	else
	{
		result = (A[i][j] -A[i-1][j])/delta - delta/2.0*minmod(derivation2x(i,j,A,delta,max_num),derivation2x(i+1,j,A,delta,max_num));

	}
	return make_tuple(result, delta);
}


tuple<double,double> derivation1y_plus(int i,int j, double** A, double delta,  int max_num)
{
	double result;
	double Aij;

	// calculate D+xPSI_ij
	Aij=A[i][j];
	if (j<max_num-1)
	{
		result = (A[i][j+1]-A[i][j])/delta - delta/2.0*minmod(derivation2y(i,j,A,delta,max_num),derivation2y(i,j+1,A,delta,max_num));
	}
	else
	{
		result = (A[i][j]-A[i][j-1])/delta - delta/2.0*minmod(derivation2y(i,j,A,delta,max_num),derivation2y(i,j+1,A,delta,max_num));

	}
	return make_tuple(result, delta);

}




////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///      computing Dx-,Dy-.
///
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// first derivation minus function
// (difference value, delta minus)
tuple<double,double> derivation1x_minus(int i,int j, double** A, double delta, int max_num)
{
	double result;
	double Aij;

	// calculate D+xPSI_ij
	Aij=A[i][j];
	if (i>0)
	{
		result = (A[i][j]-A[i-1][j])/delta + delta/2.0*minmod(derivation2x(i,j,A,delta,max_num),derivation2x(i-1,j,A,delta,max_num));
	}
	else
	{
		result = (A[i+1][j]-A[i][j])/delta +delta/2.0*minmod(derivation2x(i,j,A,delta,max_num),derivation2x(i-1,j,A,delta,max_num));
	}
	return make_tuple(result, delta);
}

tuple<double,double> derivation1y_minus(int i,int j, double** A, double delta,  int max_num)
{
	double result;
	double Aij;

	// calculate D+xPSI_ij
	Aij=A[i][j];
	if (j>0)
	{
		result = (A[i][j]-A[i][j-1])/delta + delta/2.0*minmod(derivation2y(i,j,A,delta,max_num),derivation2y(i,j-1,A,delta,max_num));
	}
	else
	{
		result = (A[i][j+1]-A[i][j])/delta + delta/2.0*minmod(derivation2y(i,j,A,delta,max_num),derivation2y(i,j-1,A,delta,max_num));
	}

	return make_tuple(result, delta);
}






////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       Near interface, computing Dx+,Dy+.
///
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// first derivation plus function near interface
// (difference value, delta plus)
tuple<double,double> derivation1_plus_near(int i,int j, double** A, double delta, char axis, int max_num)
{
	double delta_plus;
	double temp1;
	double temp2;
	double psi0_xxyy;
	double discriminant;
	double result;

	double Aij;


	switch (axis)
	{
	case 'x' : 
		if (i>0 && i<(max_num-2))
		{
			temp1 = (A[i-1][j]-2.0*A[i][j]+A[i+1][j])/(delta*delta);
			temp2 = A[i][j]-2.0*A[i+1][j]+A[i+2][j];
			psi0_xxyy = minmod(temp1,temp2);
		}
		else
		{
			psi0_xxyy = 0;
		}

		// Calculate deltaX+
		if (psi0_xxyy > DBL_EPSILON) // almost linear case
		{

			if (i<max_num-1)
			{
				discriminant = (psi0_xxyy/2.0-A[i][j]-A[i+1][j])*(psi0_xxyy/2.0-A[i][j]-A[i+1][j]) - 4.0*A[i][j]*A[i+1][j] ;
				discriminant = abs(discriminant);
				delta_plus = delta * (1.0/2.0 + A[i][j]-A[i+1][j]-sign_func(A[i][j]-A[i+1][j])*sqrt(discriminant));
			}
			else
			{
				discriminant = (psi0_xxyy/2.0-A[i][j]-A[i][j])*(psi0_xxyy/2.0-A[i][j]-A[i][j]) - 4.0*A[i][j]*A[i][j] ;
				discriminant = abs(discriminant);
				delta_plus = delta * (1.0/2.0 + A[i][j]-A[i][j]-sign_func(A[i][j]-A[i][j])*sqrt(discriminant));
			}
		}
		else // non linear case
		{
			if (i<max_num-1)
			{
				delta_plus = delta * A[i][j]/(A[i][j]-A[i+1][j]);

			}
			else
			{
				delta_plus = delta * A[i][j]/(2.0*A[i][j]);
			}
		}

		// calculate D+xPSI_ij
		Aij=A[i][j];
		result = -A[i][j]/delta_plus - delta_plus/2.0*minmod(derivation2x(i,j,A,delta,max_num),derivation2x(i+1,j,A,delta,max_num));
		return make_tuple(result, delta_plus);
		break;



	case 'y' : 
		if (j>0 && j<max_num-2)
		{
			temp1 = (A[i][j-1]-2.0*A[i][j]+A[i][j+1])/(delta*delta);
			temp2 = A[i][j]-2.0*A[i][j+1]+A[i][j+2];
			psi0_xxyy = minmod(temp1, temp2);
		}
		else
		{
			psi0_xxyy = 0;
		}

		// Calculate deltaX+
		if (psi0_xxyy > DBL_EPSILON) // almost linear case
		{

			if (j<max_num-1)
			{
				discriminant = (psi0_xxyy/2.0-A[i][j]-A[i][j+1])*(psi0_xxyy/2.0-A[i][j]-A[i][j+1]) - 4.0*A[i][j]*A[i][j+1] ;
				discriminant = abs(discriminant);
				delta_plus = delta * (1.0/2.0 + A[i][j]-A[i][j+1]-sign_func(A[i][j]-A[i][j+1])*sqrt(discriminant));
			}
			else
			{
				discriminant = (psi0_xxyy/2.0-A[i][j]-A[i][j])*(psi0_xxyy/2.0-A[i][j]-A[i][j]) - 4.0*A[i][j]*A[i][j] ;
				discriminant = abs(discriminant);
				delta_plus = delta * (1.0/2.0 + A[i][j]-A[i][j]-sign_func(A[i][j]-A[i][j])*sqrt(discriminant));
			}
		}
		else // non linear case
		{
			if (j<max_num-1)
			{
				delta_plus = delta * A[i][j]/(A[i][j]-A[i][j+1]);
			}
			else
			{
				delta_plus = delta * A[i][j]/(2.0*A[i][j]);
			}
		}

		// calculate D+xPSI_ij
		Aij=A[i][j];
		result = -A[i][j]/delta_plus - delta_plus/2.0*minmod(derivation2y(i,j,A,delta,max_num),derivation2y(i,j+1,A,delta,max_num));
		return make_tuple(result, delta_plus);
		break;


	default : cout<<"axis should be 'x' or 'y'."<<"\n";
		break;
	}

	//return result;
}







////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       Near interface, computing Dx-,Dy-.
///
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// first derivation minus function
// (difference value, delta minus)
tuple<double,double> derivation1_minus_near(int i,int j, double** A, double delta, char axis, int max_num)
{
	double delta_minus;
	double temp1;
	double temp2;
	double psi0_xxyy;
	double discriminant;
	double result;
	double Aij;


	switch (axis)
	{
	case 'x' : 
		if (i>1 && i<max_num-1)
		{
			temp1 = (A[i-1][j]-2.0*A[i][j]+A[i+1][j])/(delta*delta);
			temp2 = A[i][j]-2.0*A[i-1][j]+A[i-2][j];
			psi0_xxyy = minmod(temp1,temp2);
		}
		else
		{
			psi0_xxyy = 0;
		}

		// Calculate deltaX+
		if (psi0_xxyy > DBL_EPSILON) // almost linear case
		{

			if (i>0)
			{
				discriminant = (psi0_xxyy/2.0-A[i][j]-A[i-1][j])*(psi0_xxyy/2.0-A[i][j]-A[i-1][j]) - 4.0*A[i][j]*A[i-1][j] ;
				discriminant = abs(discriminant);
				delta_minus = delta * (1.0/2.0 + A[i][j]-A[i-1][j]-sign_func(A[i][j]-A[i-1][j])*sqrt(discriminant));
			}
			else
			{
				discriminant = (psi0_xxyy/2.0-A[i][j]-A[i][j])*(psi0_xxyy/2.0-A[i][j]-A[i][j]) - 4.0*A[i][j]*A[i][j] ;
				discriminant = abs(discriminant);
				delta_minus = delta * (1.0/2.0 + A[i][j]-A[i][j]-sign_func(A[i][j]-A[i][j])*sqrt(discriminant));
			}
		}
		else // non linear case
		{
			if (i>0)
			{
				delta_minus = delta * A[i][j]/(A[i][j]-A[i-1][j]);
			}
			else
			{
				delta_minus = delta * A[i][j]/(2.0*A[i][j]);
			}
		}

		// calculate D+xPSI_ij
		Aij=A[i][j];
		result = A[i][j]/delta_minus + delta_minus/2.0*minmod(derivation2x(i,j,A,delta,max_num),derivation2x(i-1,j,A,delta,max_num));
		return make_tuple(result, delta_minus);
		break;



	case 'y' : 
		if (j>1 && j<max_num-1)
		{
			temp1 = (A[i][j-1]-2.0*A[i][j]+A[i][j+1])/(delta*delta);
			temp2 = A[i][j]-2.0*A[i][j-1]+A[i][j-2];
			psi0_xxyy = minmod(temp1, temp2);
		}
		else
		{
			psi0_xxyy = 0;
		}

		// Calculate deltaX+
		if (psi0_xxyy > DBL_EPSILON) // almost linear case
		{

			if (j>0)
			{
				discriminant = (psi0_xxyy/2.0-A[i][j]-A[i][j-1])*(psi0_xxyy/2.0-A[i][j]-A[i][j-1]) - 4.0*A[i][j]*A[i][j-1] ;
				discriminant = abs(discriminant);
				delta_minus = delta * (1.0/2.0 + A[i][j]-A[i][j-1]-sign_func(A[i][j]-A[i][j-1])*sqrt(discriminant));
			}
			else
			{
				discriminant = (psi0_xxyy/2.0-A[i][j]-A[i][j])*(psi0_xxyy/2.0-A[i][j]-A[i][j]) - 4.0*A[i][j]*A[i][j] ;
				discriminant = abs(discriminant);
				delta_minus = delta * (1.0/2.0 + A[i][j]-A[i][j]-sign_func(A[i][j]-A[i][j])*sqrt(discriminant));
			}
		}
		else // non linear case
		{
			if (j>0)
			{
				delta_minus = delta * A[i][j]/(A[i][j]-A[i][j-1]);
			}
			else
			{
				delta_minus = delta * A[i][j]/(2.0*A[i][j]);
			}
		}

		// calculate D+xPSI_ij
		Aij=A[i][j];
		result = -A[i][j]/delta_minus - delta_minus/2.0*minmod(derivation2y(i,j,A,delta,max_num),derivation2y(i,j-1,A,delta,max_num));

		return make_tuple(result, delta_minus);
		break;


	default : cout<<"axis should be 'x' or 'y'."<<"\n";
		break;
	}


}




////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////
////                       TVD Runge-Kutta
////
////                             START
////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
void TVD_RK(level_set InputLevelSet, double total_time ,double delta_t)
{
	double cfl = 0.45;
	double dt_ij;

	level_set LevelSet1 = level_set(InputLevelSet);
	level_set LevelSet2 = level_set(InputLevelSet);
	level_set LevelSetTemp1 = level_set(InputLevelSet);
	level_set LevelSetTemp2 = level_set(InputLevelSet);

	double dt = delta_t; // time step.
	double t_total = total_time;
	int num_step = ceil(t_total/dt);

	for (int t = 0; t < num_step; t++)
	{
		cout<<"TVD Runge-Kutta time step " <<t<<"\n";

		//ostringstream file_name1;
		//file_name1 <<  "D:/TVD RK/reinitialization set temp "<<t<<".dat";
		//ofstream levelset_output1(file_name1.str());
		//assert(levelset_output1.is_open());

		// Update Psi_temp1.
		for (int i = 0; i <LevelSet1.num_xmesh; i++)
		{
			for (int j = 0; j < LevelSet1.num_ymesh; j++)
			{

				tuple<double,double> dxplus = derivation1x_plus(i,j,LevelSet1.Psi,LevelSet1.dx,LevelSet1.num_xmesh);
				tuple<double,double> dxminus = derivation1x_minus(i,j,LevelSet1.Psi,LevelSet1.dx,LevelSet1.num_xmesh);
				tuple<double,double> dyplus = derivation1y_plus(i,j,LevelSet1.Psi,LevelSet1.dy,LevelSet1.num_ymesh);
				tuple<double,double> dyminus = derivation1y_minus(i,j,LevelSet1.Psi,LevelSet1.dy,LevelSet1.num_ymesh);

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
				LevelSetTemp1.Psi[i][j] = LevelSet1.Psi[i][j] - dt_ij*sign_func(LevelSet1.Psi[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,LevelSet1.Psi[i][j])-1.0);
				//levelset_output1 << Psi_temp1[i][j] << " ";
			}

			//levelset_output1 << " "<<"\n";
		}
		//levelset_output1.close();
		// Update Psi_temp1 . Psi1.
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

		// Update Psi_temp2
		for (int i = 0; i <LevelSet1.num_xmesh; i++)
		{
			for (int j = 0; j < LevelSet1.num_ymesh; j++)
			{

				tuple<double,double> dxplus = derivation1x_plus(i,j,LevelSetTemp1.Psi,LevelSetTemp1.dx,LevelSetTemp1.num_xmesh);
				tuple<double,double> dxminus = derivation1x_minus(i,j,LevelSetTemp1.Psi,LevelSetTemp1.dx,LevelSetTemp1.num_xmesh);
				tuple<double,double> dyplus = derivation1y_plus(i,j,LevelSetTemp1.Psi,LevelSetTemp1.dy,LevelSetTemp1.num_ymesh);
				tuple<double,double> dyminus = derivation1y_minus(i,j,LevelSetTemp1.Psi,LevelSetTemp1.dy,LevelSetTemp1.num_ymesh);
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


				LevelSetTemp2.Psi[i][j] = LevelSetTemp1.Psi[i][j] - dt_ij*sign_func(LevelSetTemp1.Psi[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,LevelSetTemp1.Psi[i][j])-1.0);

				LevelSet2.Psi[i][j]= (LevelSet1.Psi[i][j]+LevelSetTemp2.Psi[i][j])/2.0;

				levelset_output2 << LevelSet2.Psi[i][j] << " ";

			} //end 'for' with j

			levelset_output2 << " "<<"\n";
		} // emd 'for' with i
		levelset_output2.close();

		// Update Psi2 . Psi1.
		LevelSet1 = LevelSet2;

	} // end 'for' with t
}



//void TVD_RK_with_error(level_set InputLevelSet, double total_time ,double delta_t, class level_set sign_distance_level_set)
//{
//	double l0_error;
//	double l1_error;
//	double error;
//	double cfl = 0.45;
//	double dt_ij;
//	double* L1;
//	double* L0;
//
//	level_set LevelSet1= LevelSet;
//	level_set LevelSet2=LevelSet;
//	level_set LevelSetTemp1=LevelSet;
//	level_set LevelSetTemp2=LevelSet;
//
//	double dt = delta_t; // time step.
//	double t_total = total_time;
//	int num_step = ceil(t_total/dt);
//
//	L1 = new double [num_step];
//	L0 = new double [num_step];
//
//	for (int t = 0; t < num_step; t++)
//	{
//		cout<<"TVD Runge-Kutta time step " <<t<<"\n";
//
//		l1_error = 0;
//		l0_error = 0;
//		error = 0;
//		//ostringstream file_name1;
//		//file_name1 <<  "D:/TVD RK/reinitialization set temp "<<t<<".dat";
//		//ofstream levelset_output1(file_name1.str());
//		//assert(levelset_output1.is_open());
//
//		// Update Psi_temp1.
//		for (int i = 0; i <LevelSet1.num_xmesh; i++)
//		{
//			for (int j = 0; j < LevelSet1.num_ymesh; j++)
//			{
//
//				tuple<double,double> dxplus = derivation1x_plus(i,j,LevelSet1.Psi,LevelSet1.dx,LevelSet1.num_xmesh);
//				tuple<double,double> dxminus = derivation1x_minus(i,j,LevelSet1.Psi,LevelSet1.dx,LevelSet1.num_xmesh);
//				tuple<double,double> dyplus = derivation1y_plus(i,j,LevelSet1.Psi,LevelSet1.dy,LevelSet1.num_ymesh);
//				tuple<double,double> dyminus = derivation1y_minus(i,j,LevelSet1.Psi,LevelSet1.dy,LevelSet1.num_ymesh);
//
//				double derv_xplus, derv_xminus, derv_yplus, derv_yminus;
//				double delta1, delta2, delta3, delta4;
//
//				// D+-xy .   i.e. difference value
//				derv_xplus = get<0>(dxplus) ;
//				derv_xminus = get<0>(dxminus) ;
//				derv_yplus = get<0>(dyplus) ;
//				derv_yminus = get<0>(dyminus) ;
//				// delta x+,delta x-,delta y+,delta y-.
//				delta1 = get<1>(dxplus) ;
//				delta2 = get<1>(dxminus) ;
//				delta3 = get<1>(dyplus) ;
//				delta4 = get<1>(dyminus) ;
//
//				dt_ij = cfl*min(min(delta1, delta2),min(delta3,delta4));
//
//				//Psi_temp1[i][j] = Psi1[i][j] - dt_ij*sign_func(Psi1[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,Psi1[i][j])-1.0);
//				LevelSetTemp1.Psi[i][j] = LevelSet1.Psi[i][j] - dt_ij*sign_func(LevelSet1.Psi[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,LevelSet1.Psi[i][j])-1.0);
//				//levelset_output1 << Psi_temp1[i][j] << " ";
//			}
//
//			//levelset_output1 << " "<<"\n";
//		}
//		//levelset_output1.close();
//		// Update Psi_temp1 . Psi1.
//		//Psi1=Psi_temp1;
//
//
//
//		//////////////////////////////////////////////////////////////////////////////////////////////////
//		//
//		//                     Psi2 is updated from Psi_temp1
//		//
//		///////////////////////////////////////////////////////////////////////////////////////////////////
//		ostringstream file_name2;
//		file_name2 <<  "D:/TVD RK/reinitialization set "<<t<<".dat";
//		ofstream levelset_output2(file_name2.str());
//		assert(levelset_output2.is_open());
//
//		ostringstream file_name3;
//		file_name3 <<  "D:/TVD RK/Error "<<t<<".dat";
//		ofstream levelset_output3(file_name3.str());
//		assert(levelset_output3.is_open());
//
//		// Update Psi_temp2
//		for (int i = 0; i <LevelSet1.num_xmesh; i++)
//		{
//			for (int j = 0; j < LevelSet1.num_ymesh; j++)
//			{
//
//				tuple<double,double> dxplus = derivation1x_plus(i,j,LevelSetTemp1.Psi,LevelSetTemp1.dx,LevelSetTemp1.num_xmesh);
//				tuple<double,double> dxminus = derivation1x_minus(i,j,LevelSetTemp1.Psi,LevelSetTemp1.dx,LevelSetTemp1.num_xmesh);
//				tuple<double,double> dyplus = derivation1y_plus(i,j,LevelSetTemp1.Psi,LevelSetTemp1.dy,LevelSetTemp1.num_ymesh);
//				tuple<double,double> dyminus = derivation1y_minus(i,j,LevelSetTemp1.Psi,LevelSetTemp1.dy,LevelSetTemp1.num_ymesh);
//				double derv_xplus,derv_xminus,derv_yplus,derv_yminus;
//				double delta1,delta2,delta3,delta4;
//
//				// D+-xy .   i.e. difference value
//				derv_xplus = get<0>(dxplus) ;
//				derv_xminus = get<0>(dxminus) ;
//				derv_yplus = get<0>(dyplus) ;
//				derv_yminus = get<0>(dyminus) ;
//				// delta x+,delta x-,delta y+,delta y-.
//				delta1 = get<1>(dxplus) ;
//				delta2 = get<1>(dxminus) ;
//				delta3 = get<1>(dyplus) ;
//				delta4 = get<1>(dyminus) ;
//
//				dt_ij = cfl * min(min(delta1, delta2),min(delta3,delta4));
//
//
//				LevelSetTemp2.Psi[i][j] = LevelSetTemp1.Psi[i][j] - dt_ij*sign_func(LevelSetTemp1.Psi[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,LevelSetTemp1.Psi[i][j])-1.0);
//
//				LevelSet2.Psi[i][j]= (LevelSet1.Psi[i][j]+LevelSetTemp2.Psi[i][j])/2.0;
//
//				levelset_output2 << LevelSet2.Psi[i][j] << " ";
//				error =sign_distance_level_set.Psi[i][j] - LevelSet2.Psi[i][j];
//				levelset_output3 << error << " ";
//
//				// Compute L1 and L2 error.
//				if (l0_error<abs(error))
//				{
//					l0_error = abs(error);
//				}
//				l1_error += abs(error);
//
//			} //end 'for' with j
//
//			levelset_output2 << " "<<"\n";
//			levelset_output3 << " "<<"\n";
//		} // emd 'for' with i
//		levelset_output2.close();
//		levelset_output3.close();
//
//
//		L1[t] = l1_error*(LevelSet1.dx*LevelSet1.dy);
//		L0[t] = l0_error;
//
//		ostringstream file_name4;
//		file_name4 <<  "D:/TVD RK/L0 "<<t<<".dat";
//		ofstream L0_output4(file_name4.str());
//		assert(L0_output4.is_open());
//		L0_output4 << L0[t] << " ";
//		L0_output4.close();
//
//		ostringstream file_name5;
//		file_name5 <<  "D:/TVD RK/L1 "<<t<<".dat";
//		ofstream L1_output(file_name5.str());
//		assert(L1_output.is_open());
//		L1_output << L1[t] << " ";
//		L1_output.close();
//
//		// Update Psi2 . Psi1.
//		LevelSet1 = LevelSet2;
//	} // end 'for' with t
//}
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////
////                       TVD Runge-Kutta
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
void FE(level_set LevelSet, double total_time ,double delta_t)
{
	double cfl = 0.45;
	double dt_ij;
	double* L1;

	level_set LevelSet1=level_set(LevelSet);
	level_set LevelSet2=level_set(LevelSet);
	level_set LevelSetTemp1=level_set(LevelSet);
	level_set LevelSetTemp2=level_set(LevelSet);

	double dt = delta_t; // time step.
	double t_total = total_time;
	int num_step = ceil(t_total/dt);

	for (int t = 0; t < num_step; t++)
	{
		cout<<"Forward Euler time step " <<t<<"\n";

		//////////////////////////////////////////////////////////////////////////////////////////////////
		//
		//                     psi2 is updated from psi1
		//
		///////////////////////////////////////////////////////////////////////////////////////////////////
		ostringstream file_name1;
		file_name1 <<  "D:/FE/reinitialization set "<<t<<".dat";
		ofstream levelset_output1(file_name1.str());
		assert(levelset_output1.is_open());

		// update psi2.
		for (int i = 0; i < LevelSet1.num_xmesh; i++)
		{
			for (int j = 0; j < LevelSet1.num_ymesh; j++)
			{

				tuple<double,double> dxplus = derivation1x_plus(i,j,LevelSet1.Psi,LevelSet1.dx,LevelSet1.num_xmesh);
				tuple<double,double> dxminus = derivation1x_minus(i,j,LevelSet1.Psi,LevelSet1.dx,LevelSet1.num_xmesh);
				tuple<double,double> dyplus = derivation1y_plus(i,j,LevelSet1.Psi,LevelSet1.dy,LevelSet1.num_ymesh);
				tuple<double,double> dyminus = derivation1y_minus(i,j,LevelSet1.Psi,LevelSet1.dy,LevelSet1.num_ymesh);

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
				LevelSet2.Psi[i][j] = LevelSet1.Psi[i][j] - dt_ij*sign_func(LevelSet1.Psi[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,LevelSet1.Psi[i][j])-1.0);
				levelset_output1 << LevelSet2.Psi[i][j] << " ";

				// Compute L1 and L2 error.

			}

			levelset_output1 << " "<<"\n";
		}
		levelset_output1.close();

		// update psi_temp1 . psi1.
		LevelSet1 = level_set( LevelSet2);
	}
}




void FE_with_error(level_set LevelSet, double total_time ,double delta_t, class level_set sign_distance_level_set)
{
	double l0_error;
	double l1_error;
	double error;
	double cfl = 0.45;
	double dt_ij;
	double* L1;
	double* L0;

	level_set LevelSet1=LevelSet;
	level_set LevelSet2=LevelSet;
	level_set LevelSetTemp1=LevelSet;
	level_set LevelSetTemp2=LevelSet;

	double dt = delta_t; // time step.
	double t_total = total_time;
	int num_step = ceil(t_total/dt);

	L1 = new double [num_step];
	L0 = new double [num_step];

	for (int t = 0; t < num_step; t++)
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
		for (int i = 0; i < LevelSet1.num_xmesh; i++)
		{
			for (int j = 0; j < LevelSet1.num_ymesh; j++)
			{

				tuple<double,double> dxplus = derivation1x_plus(i,j,LevelSet1.Psi,LevelSet1.dx,LevelSet1.num_xmesh);
				tuple<double,double> dxminus = derivation1x_minus(i,j,LevelSet1.Psi,LevelSet1.dx,LevelSet1.num_xmesh);
				tuple<double,double> dyplus = derivation1y_plus(i,j,LevelSet1.Psi,LevelSet1.dy,LevelSet1.num_ymesh);
				tuple<double,double> dyminus = derivation1y_minus(i,j,LevelSet1.Psi,LevelSet1.dy,LevelSet1.num_ymesh);

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
				LevelSet2.Psi[i][j] = LevelSet1.Psi[i][j] - dt_ij*sign_func(LevelSet1.Psi[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,LevelSet1.Psi[i][j])-1.0);
				levelset_output1 << LevelSet2.Psi[i][j] << " ";
				error = sign_distance_level_set.Psi[i][j] - LevelSet2.Psi[i][j];
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


		L1[t] = l1_error*(LevelSet1.dx*LevelSet1.dy);
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


		// update psi_temp1 . psi1.
		LevelSet1.Psi = LevelSet2.Psi;
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////
////                      Forward Euler method
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
//void GS(level_set LevelSet, double total_time ,double delta_t)
//{
//	double cfl = 0.45;
//	double dt_ij;
//
//	level_set LevelSet1=LevelSet;
//	level_set LevelSet2=LevelSet;
//	level_set LevelSetTemp1=LevelSet;
//	level_set LevelSetTemp2=LevelSet;
//
//	double dt = delta_t; // time step.
//	double t_total = total_time;
//	int num_step = ceil(t_total/dt);
//
//	for (int t = 0; t < num_step; t++)
//	{
//		cout<<"Gauss-Seidel iteration time step " <<t<<"\n";
//
//		int raster_visiting;
//		raster_visiting = t%4;
//
//
//		//////////////////////////////////////////////////////////////////////////////////////////////////
//		//
//		//                     Psi1 is updated from Psi1
//		//
//		///////////////////////////////////////////////////////////////////////////////////////////////////
//		ostringstream file_name1;
//		file_name1 <<  "D:/Gauss Seidel/reinitialization set "<<t<<".dat";
//		ofstream levelset_output1(file_name1.str());
//		assert(levelset_output1.is_open());
//
//		// Update Psi1, itself, using four raster-scan visiting.
//		switch (raster_visiting)
//		{
//
//		case 0 : 
//			for (int i = 0; i < LevelSet1.num_xmesh; i++)
//			{
//				for (int j = 0; j < LevelSet1.num_ymesh; j++)
//				{
//					tuple<double,double> dxplus = derivation1x_plus(i,j,LevelSet1.Psi,LevelSet1.dx,LevelSet1.num_xmesh);
//					tuple<double,double> dxminus = derivation1x_minus(i,j,LevelSet1.Psi,LevelSet1.dx,LevelSet1.num_xmesh);
//					tuple<double,double> dyplus = derivation1y_plus(i,j,LevelSet1.Psi,LevelSet1.dy,LevelSet1.num_ymesh);
//					tuple<double,double> dyminus = derivation1y_minus(i,j,LevelSet1.Psi,LevelSet1.dy,LevelSet1.num_ymesh);
//
//					double derv_xplus, derv_xminus, derv_yplus, derv_yminus;
//					double delta1, delta2, delta3, delta4;
//
//					// D+-xy .   i.e. difference value
//					derv_xplus = get<0>(dxplus) ;
//					derv_xminus = get<0>(dxminus) ;
//					derv_yplus = get<0>(dyplus) ;
//					derv_yminus = get<0>(dyminus) ;
//					// delta x+,delta x-,delta y+,delta y-.
//					delta1 = get<1>(dxplus) ;
//					delta2 = get<1>(dxminus) ;
//					delta3 = get<1>(dyplus) ;
//					delta4 = get<1>(dyminus) ;
//
//					dt_ij = cfl*min(min(delta1, delta2),min(delta3,delta4));
//
//					//Psi_temp1[i][j] = Psi1[i][j] - dt_ij*sign_func(Psi1[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,Psi1[i][j])-1.0);
//					LevelSet1.Psi[i][j] = LevelSet1.Psi[i][j] - dt_ij*sign_func(LevelSet1.Psi[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,LevelSet1.Psi[i][j])-1.0);
//
//					levelset_output1 << LevelSet1.Psi[i][j] << " ";
//				} //end 'for' with j
//
//				levelset_output1 << " "<<"\n";
//			} //end 'for' with i
//			break;
//
//
//
//
//		case 1 : 
//			for (int i = 0; i < LevelSet1.num_xmesh; i++)
//			{
//				for (int j = LevelSet1.num_ymesh-1; j > -1; j--)
//				{
//					tuple<double,double> dxplus = derivation1x_plus(i,j,LevelSet1.Psi,LevelSet1.dx,LevelSet1.num_xmesh);
//					tuple<double,double> dxminus = derivation1x_minus(i,j,LevelSet1.Psi,LevelSet1.dx,LevelSet1.num_xmesh);
//					tuple<double,double> dyplus = derivation1y_plus(i,j,LevelSet1.Psi,LevelSet1.dy,LevelSet1.num_ymesh);
//					tuple<double,double> dyminus = derivation1y_minus(i,j,LevelSet1.Psi,LevelSet1.dy,LevelSet1.num_ymesh);
//
//					double derv_xplus, derv_xminus, derv_yplus, derv_yminus;
//					double delta1, delta2, delta3, delta4;
//
//					// D+-xy .   i.e. difference value
//					derv_xplus = get<0>(dxplus) ;
//					derv_xminus = get<0>(dxminus) ;
//					derv_yplus = get<0>(dyplus) ;
//					derv_yminus = get<0>(dyminus) ;
//					// delta x+,delta x-,delta y+,delta y-.
//					delta1 = get<1>(dxplus) ;
//					delta2 = get<1>(dxminus) ;
//					delta3 = get<1>(dyplus) ;
//					delta4 = get<1>(dyminus) ;
//
//					dt_ij = cfl*min(min(delta1, delta2),min(delta3,delta4));
//
//					//Psi_temp1[i][j] = LevelSet1.Psi[i][j] - dt_ij*sign_func(LevelSet1.Psi[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,LevelSet1.Psi[i][j])-1.0);
//					LevelSet1.Psi[i][j] = LevelSet1.Psi[i][j] - dt_ij*sign_func(LevelSet1.Psi[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,LevelSet1.Psi[i][j])-1.0);
//					/*levelset_output1 << LevelSet1.Psi[i][j] << " ";
//					levelset_output3 << sign_distance_level_set[i][j] - LevelSet1.Psi[i][j] << " ";*/
//				} //end 'for' with j
//
//				for (int j = 0; j  < LevelSet1.num_ymesh; j ++)
//				{
//					levelset_output1 << LevelSet1.Psi[i][j] << " ";
//				}
//
//				levelset_output1 << " "<<"\n";
//			} //end 'for' with i
//			break;
//
//
//
//
//		case 2 : 
//			for (int i = LevelSet1.num_xmesh-1; i>-1 ; i--)
//			{
//				for (int j = 0; j < LevelSet1.num_ymesh; j++)
//				{
//					tuple<double,double> dxplus = derivation1x_plus(i,j,LevelSet1.Psi,LevelSet1.dx,LevelSet1.num_xmesh);
//					tuple<double,double> dxminus = derivation1x_minus(i,j,LevelSet1.Psi,LevelSet1.dx,LevelSet1.num_xmesh);
//					tuple<double,double> dyplus = derivation1y_plus(i,j,LevelSet1.Psi,LevelSet1.dy,LevelSet1.num_ymesh);
//					tuple<double,double> dyminus = derivation1y_minus(i,j,LevelSet1.Psi,LevelSet1.dy,LevelSet1.num_ymesh);
//
//					double derv_xplus, derv_xminus, derv_yplus, derv_yminus;
//					double delta1, delta2, delta3, delta4;
//
//					// D+-xy .   i.e. difference value
//					derv_xplus = get<0>(dxplus) ;
//					derv_xminus = get<0>(dxminus) ;
//					derv_yplus = get<0>(dyplus) ;
//					derv_yminus = get<0>(dyminus) ;
//					// delta x+,delta x-,delta y+,delta y-.
//					delta1 = get<1>(dxplus) ;
//					delta2 = get<1>(dxminus) ;
//					delta3 = get<1>(dyplus) ;
//					delta4 = get<1>(dyminus) ;
//
//					dt_ij = cfl*min(min(delta1, delta2),min(delta3,delta4));
//
//					//Psi_temp1[i][j] = LevelSet1.Psi[i][j] - dt_ij*sign_func(LevelSet1.Psi[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,LevelSet1.Psi[i][j])-1.0);
//					LevelSet1.Psi[i][j] = LevelSet1.Psi[i][j] - dt_ij*sign_func(LevelSet1.Psi[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,LevelSet1.Psi[i][j])-1.0);
//					//levelset_output1 << LevelSet1.Psi[i][j] << " ";
//					//levelset_output3 << sign_distance_level_set[i][j] - LevelSet1.Psi[i][j] << " ";
//				} //end 'for' with j
//				//levelset_output1 << " "<<"\n";
//				//levelset_output3 << " "<<"\n";
//			} //end 'for' with i
//
//			for (int i = 0; i < LevelSet1.num_xmesh; i++)
//			{
//				for (int j = 0; j < LevelSet1.num_ymesh; j++)
//				{
//					levelset_output1 << LevelSet1.Psi[i][j] << " ";
//				}
//				levelset_output1 << " "<<"\n";
//			}
//			break;
//
//
//
//		case 3 : 
//			for (int i = LevelSet1.num_xmesh-1; i>-1 ; i--)
//			{
//				for (int j = LevelSet1.num_ymesh-1; j > -1; j--)
//				{
//
//					tuple<double,double> dxplus = derivation1x_plus(i,j,LevelSet1.Psi,LevelSet1.dx,LevelSet1.num_xmesh);
//					tuple<double,double> dxminus = derivation1x_minus(i,j,LevelSet1.Psi,LevelSet1.dx,LevelSet1.num_xmesh);
//					tuple<double,double> dyplus = derivation1y_plus(i,j,LevelSet1.Psi,LevelSet1.dy,LevelSet1.num_ymesh);
//					tuple<double,double> dyminus = derivation1y_minus(i,j,LevelSet1.Psi,LevelSet1.dy,LevelSet1.num_ymesh);
//
//					double derv_xplus, derv_xminus, derv_yplus, derv_yminus;
//					double delta1, delta2, delta3, delta4;
//
//					// D+-xy .   i.e. difference value
//					derv_xplus = get<0>(dxplus) ;
//					derv_xminus = get<0>(dxminus) ;
//					derv_yplus = get<0>(dyplus) ;
//					derv_yminus = get<0>(dyminus) ;
//					// delta x+,delta x-,delta y+,delta y-.
//					delta1 = get<1>(dxplus) ;
//					delta2 = get<1>(dxminus) ;
//					delta3 = get<1>(dyplus) ;
//					delta4 = get<1>(dyminus) ;
//
//					dt_ij = cfl*min(min(delta1, delta2),min(delta3,delta4));
//
//					//Psi_temp1[i][j] = LevelSet1.Psi[i][j] - dt_ij*sign_func(LevelSet1.Psi[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,LevelSet1.Psi[i][j])-1.0);
//					LevelSet1.Psi[i][j] = LevelSet1.Psi[i][j] - dt_ij*sign_func(LevelSet1.Psi[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,LevelSet1.Psi[i][j])-1.0);
//
//					//levelset_output1 << LevelSet1.Psi[i][j] << " ";
//					//levelset_output3 << sign_distance_level_set[i][j] - LevelSet1.Psi[i][j] << " ";
//				} //end 'for' with j
//				//levelset_output1 << " "<<"\n";
//				//levelset_output3 << " "<<"\n";
//			} //end 'for' with i
//
//			for (int i = 0; i < LevelSet1.num_xmesh; i++)
//			{
//				for (int j = 0; j < LevelSet1.num_ymesh; j++)
//				{
//					levelset_output1 << LevelSet1.Psi[i][j] << " ";
//				}
//				levelset_output1 << " "<<"\n";
//			}
//			break;
//
//		} // end 'switch'
//
//		levelset_output1.close();
//
//	}//end 'for' with t
//}
//
//
//
//
//void GS_with_error(level_set LevelSet, double total_time ,double delta_t, class level_set sign_distance_level_set)
//{
//	double l0_error;
//	double l1_error;
//	double error;
//	double cfl = 0.45;
//	double dt_ij;
//	double* L1;
//	double* L0;
//
//	level_set LevelSet1=LevelSet;
//	level_set LevelSet2=LevelSet;
//	level_set LevelSetTemp1=LevelSet;
//	level_set LevelSetTemp2=LevelSet;
//
//	double dt = delta_t; // time step.
//	double t_total = total_time;
//	int num_step = ceil(t_total/dt);
//
//	L1 = new double [num_step];
//	L0 = new double [num_step];
//
//	for (int t = 0; t < num_step; t++)
//	{
//		cout<<"Gauss-Seidel iteration time step " <<t<<"\n";
//
//		l1_error = 0;
//		l0_error = 0;
//		error = 0;
//		int raster_visiting;
//		raster_visiting = t%4;
//
//
//		//////////////////////////////////////////////////////////////////////////////////////////////////
//		//
//		//                     Psi1 is updated from Psi1
//		//
//		///////////////////////////////////////////////////////////////////////////////////////////////////
//		ostringstream file_name1;
//		file_name1 <<  "D:/Gauss Seidel/reinitialization set "<<t<<".dat";
//		ofstream levelset_output1(file_name1.str());
//		assert(levelset_output1.is_open());
//
//		ostringstream file_name3;
//		file_name3 <<  "D:/Gauss Seidel/Error "<<t<<".dat";
//		ofstream levelset_output3(file_name3.str());
//		assert(levelset_output3.is_open());
//
//		// Update Psi1, itself.
//		switch (raster_visiting)
//		{
//
//		case 0 : 
//			for (int i = 0; i < LevelSet1.num_xmesh; i++)
//			{
//				for (int j = 0; j < LevelSet1.num_ymesh; j++)
//				{
//					tuple<double,double> dxplus = derivation1x_plus(i,j,LevelSet1.Psi,LevelSet1.dx,LevelSet1.num_xmesh);
//					tuple<double,double> dxminus = derivation1x_minus(i,j,LevelSet1.Psi,LevelSet1.dx,LevelSet1.num_xmesh);
//					tuple<double,double> dyplus = derivation1y_plus(i,j,LevelSet1.Psi,LevelSet1.dy,LevelSet1.num_ymesh);
//					tuple<double,double> dyminus = derivation1y_minus(i,j,LevelSet1.Psi,LevelSet1.dy,LevelSet1.num_ymesh);
//
//					double derv_xplus, derv_xminus, derv_yplus, derv_yminus;
//					double delta1, delta2, delta3, delta4;
//
//					// D+-xy .   i.e. difference value
//					derv_xplus = get<0>(dxplus) ;
//					derv_xminus = get<0>(dxminus) ;
//					derv_yplus = get<0>(dyplus) ;
//					derv_yminus = get<0>(dyminus) ;
//					// delta x+,delta x-,delta y+,delta y-.
//					delta1 = get<1>(dxplus) ;
//					delta2 = get<1>(dxminus) ;
//					delta3 = get<1>(dyplus) ;
//					delta4 = get<1>(dyminus) ;
//
//					dt_ij = cfl*min(min(delta1, delta2),min(delta3,delta4));
//
//					//Psi_temp1[i][j] = Psi1[i][j] - dt_ij*sign_func(Psi1[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,Psi1[i][j])-1.0);
//					LevelSet1.Psi[i][j] = LevelSet1.Psi[i][j] - dt_ij*sign_func(LevelSet1.Psi[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,LevelSet1.Psi[i][j])-1.0);
//
//					levelset_output1 << LevelSet1.Psi[i][j] << " ";
//					error = sign_distance_level_set.Psi[i][j] - LevelSet1.Psi[i][j];
//					levelset_output3 << error << " ";
//
//					// Compute L1 and L2 error.
//					if (l0_error<abs(error))
//					{
//						l0_error = abs(error);
//					}
//					l1_error += abs(error);
//
//				} //end 'for' with j
//
//				levelset_output1 << " "<<"\n";
//				levelset_output3 << " "<<"\n";
//			} //end 'for' with i
//			break;
//
//
//
//
//		case 1 : 
//			for (int i = 0; i < LevelSet1.num_xmesh; i++)
//			{
//				for (int j = LevelSet1.num_ymesh-1; j > -1; j--)
//				{
//					tuple<double,double> dxplus = derivation1x_plus(i,j,LevelSet1.Psi,LevelSet1.dx,LevelSet1.num_xmesh);
//					tuple<double,double> dxminus = derivation1x_minus(i,j,LevelSet1.Psi,LevelSet1.dx,LevelSet1.num_xmesh);
//					tuple<double,double> dyplus = derivation1y_plus(i,j,LevelSet1.Psi,LevelSet1.dy,LevelSet1.num_ymesh);
//					tuple<double,double> dyminus = derivation1y_minus(i,j,LevelSet1.Psi,LevelSet1.dy,LevelSet1.num_ymesh);
//
//					double derv_xplus, derv_xminus, derv_yplus, derv_yminus;
//					double delta1, delta2, delta3, delta4;
//
//					// D+-xy .   i.e. difference value
//					derv_xplus = get<0>(dxplus) ;
//					derv_xminus = get<0>(dxminus) ;
//					derv_yplus = get<0>(dyplus) ;
//					derv_yminus = get<0>(dyminus) ;
//					// delta x+,delta x-,delta y+,delta y-.
//					delta1 = get<1>(dxplus) ;
//					delta2 = get<1>(dxminus) ;
//					delta3 = get<1>(dyplus) ;
//					delta4 = get<1>(dyminus) ;
//
//					dt_ij = cfl*min(min(delta1, delta2),min(delta3,delta4));
//
//					//Psi_temp1[i][j] = LevelSet1.Psi[i][j] - dt_ij*sign_func(LevelSet1.Psi[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,LevelSet1.Psi[i][j])-1.0);
//					LevelSet1.Psi[i][j] = LevelSet1.Psi[i][j] - dt_ij*sign_func(LevelSet1.Psi[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,LevelSet1.Psi[i][j])-1.0);
//					/*levelset_output1 << LevelSet1.Psi[i][j] << " ";
//					levelset_output3 << sign_distance_level_set[i][j] - LevelSet1.Psi[i][j] << " ";*/
//				} //end 'for' with j
//
//				for (int j = 0; j  < LevelSet1.num_ymesh; j ++)
//				{
//					levelset_output1 << LevelSet1.Psi[i][j] << " ";
//					error = sign_distance_level_set.Psi[i][j] - LevelSet1.Psi[i][j];
//					levelset_output3 << error << " ";
//
//					// Compute L1 and L2 error.
//					if (l0_error<abs(error))
//					{
//						l0_error = abs(error);
//					}
//					l1_error += abs(error);
//				}
//
//				levelset_output1 << " "<<"\n";
//				levelset_output3 << " "<<"\n";
//			} //end 'for' with i
//			break;
//
//
//
//
//		case 2 : 
//			for (int i = LevelSet1.num_xmesh-1; i>-1 ; i--)
//			{
//				for (int j = 0; j < LevelSet1.num_ymesh; j++)
//				{
//					tuple<double,double> dxplus = derivation1x_plus(i,j,LevelSet1.Psi,LevelSet1.dx,LevelSet1.num_xmesh);
//					tuple<double,double> dxminus = derivation1x_minus(i,j,LevelSet1.Psi,LevelSet1.dx,LevelSet1.num_xmesh);
//					tuple<double,double> dyplus = derivation1y_plus(i,j,LevelSet1.Psi,LevelSet1.dy,LevelSet1.num_ymesh);
//					tuple<double,double> dyminus = derivation1y_minus(i,j,LevelSet1.Psi,LevelSet1.dy,LevelSet1.num_ymesh);
//
//					double derv_xplus, derv_xminus, derv_yplus, derv_yminus;
//					double delta1, delta2, delta3, delta4;
//
//					// D+-xy .   i.e. difference value
//					derv_xplus = get<0>(dxplus) ;
//					derv_xminus = get<0>(dxminus) ;
//					derv_yplus = get<0>(dyplus) ;
//					derv_yminus = get<0>(dyminus) ;
//					// delta x+,delta x-,delta y+,delta y-.
//					delta1 = get<1>(dxplus) ;
//					delta2 = get<1>(dxminus) ;
//					delta3 = get<1>(dyplus) ;
//					delta4 = get<1>(dyminus) ;
//
//					dt_ij = cfl*min(min(delta1, delta2),min(delta3,delta4));
//
//					//Psi_temp1[i][j] = LevelSet1.Psi[i][j] - dt_ij*sign_func(LevelSet1.Psi[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,LevelSet1.Psi[i][j])-1.0);
//					LevelSet1.Psi[i][j] = LevelSet1.Psi[i][j] - dt_ij*sign_func(LevelSet1.Psi[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,LevelSet1.Psi[i][j])-1.0);
//					//levelset_output1 << LevelSet1.Psi[i][j] << " ";
//					//levelset_output3 << sign_distance_level_set[i][j] - LevelSet1.Psi[i][j] << " ";
//				} //end 'for' with j
//				//levelset_output1 << " "<<"\n";
//				//levelset_output3 << " "<<"\n";
//			} //end 'for' with i
//
//			for (int i = 0; i < LevelSet1.num_xmesh; i++)
//			{
//				for (int j = 0; j < LevelSet1.num_ymesh; j++)
//				{
//					levelset_output1 << LevelSet1.Psi[i][j] << " ";
//					error = sign_distance_level_set.Psi[i][j] - LevelSet1.Psi[i][j];
//					levelset_output3 << error << " ";
//
//					// Compute L1 and L2 error.
//					if (l0_error<abs(error))
//					{
//						l0_error = abs(error);
//					}
//					l1_error += abs(error);
//				}
//				levelset_output1 << " "<<"\n";
//				levelset_output3 << " "<<"\n";
//			}
//			break;
//
//
//
//		case 3 : 
//			for (int i = LevelSet1.num_xmesh-1; i>-1 ; i--)
//			{
//				for (int j = LevelSet1.num_ymesh-1; j > -1; j--)
//				{
//
//					tuple<double,double> dxplus = derivation1x_plus(i,j,LevelSet1.Psi,LevelSet1.dx,LevelSet1.num_xmesh);
//					tuple<double,double> dxminus = derivation1x_minus(i,j,LevelSet1.Psi,LevelSet1.dx,LevelSet1.num_xmesh);
//					tuple<double,double> dyplus = derivation1y_plus(i,j,LevelSet1.Psi,LevelSet1.dy,LevelSet1.num_ymesh);
//					tuple<double,double> dyminus = derivation1y_minus(i,j,LevelSet1.Psi,LevelSet1.dy,LevelSet1.num_ymesh);
//
//					double derv_xplus, derv_xminus, derv_yplus, derv_yminus;
//					double delta1, delta2, delta3, delta4;
//
//					// D+-xy .   i.e. difference value
//					derv_xplus = get<0>(dxplus) ;
//					derv_xminus = get<0>(dxminus) ;
//					derv_yplus = get<0>(dyplus) ;
//					derv_yminus = get<0>(dyminus) ;
//					// delta x+,delta x-,delta y+,delta y-.
//					delta1 = get<1>(dxplus) ;
//					delta2 = get<1>(dxminus) ;
//					delta3 = get<1>(dyplus) ;
//					delta4 = get<1>(dyminus) ;
//
//					dt_ij = cfl*min(min(delta1, delta2),min(delta3,delta4));
//
//					//Psi_temp1[i][j] = LevelSet1.Psi[i][j] - dt_ij*sign_func(LevelSet1.Psi[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,LevelSet1.Psi[i][j])-1.0);
//					LevelSet1.Psi[i][j] = LevelSet1.Psi[i][j] - dt_ij*sign_func(LevelSet1.Psi[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,LevelSet1.Psi[i][j])-1.0);
//
//					//levelset_output1 << LevelSet1.Psi[i][j] << " ";
//					//levelset_output3 << sign_distance_level_set[i][j] - LevelSet1.Psi[i][j] << " ";
//				} //end 'for' with j
//				//levelset_output1 << " "<<"\n";
//				//levelset_output3 << " "<<"\n";
//			} //end 'for' with i
//
//			for (int i = 0; i < LevelSet1.num_xmesh; i++)
//			{
//				for (int j = 0; j < LevelSet1.num_ymesh; j++)
//				{
//					levelset_output1 << LevelSet1.Psi[i][j] << " ";
//					error = sign_distance_level_set.Psi[i][j] - LevelSet1.Psi[i][j];
//					levelset_output3 << error << " ";
//
//					// Compute L1 and L2 error.
//					if (l0_error<abs(error))
//					{
//						l0_error = abs(error);
//					}
//					l1_error += abs(error);
//				}
//				levelset_output1 << " "<<"\n";
//				levelset_output3 << " "<<"\n";
//			}
//			break;
//
//		} // end 'switch'
//
//		levelset_output1.close();
//
//
//		L1[t] = l1_error*(LevelSet1.dx*LevelSet1.dy);
//		L0[t] = l0_error;
//
//		ostringstream file_name4;
//		file_name4 <<  "D:/Gauss Seidel/L0 "<<t<<".dat";
//		ofstream L0_output4(file_name4.str());
//		assert(L0_output4.is_open());
//		L0_output4 << L0[t] << " ";
//		L0_output4.close();
//
//		ostringstream file_name5;
//		file_name5 <<  "D:/Gauss Seidel/L1 "<<t<<".dat";
//		ofstream L1_output(file_name5.str());
//		assert(L1_output.is_open());
//		L1_output << L1[t] << " ";
//		L1_output.close();
//
//	}//end 'for' with t
//}


////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////
////                       Gauss-Seidel iteration
////
////                             END
////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////