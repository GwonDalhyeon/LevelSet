#include "reinitialization_Function.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////
////                       TVD Runge-Kutta
////
////                             START
////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
void TVD_RK(LEVELSET InputLevelSet, double total_time ,double delta_t)
{
	double cfl = 0.45;
	double dt_ij;

	LEVELSET LevelSet1 = InputLevelSet;
	LEVELSET LevelSet2 = InputLevelSet;
	LEVELSET LevelSetTemp1 = InputLevelSet;
	LEVELSET LevelSetTemp2 = InputLevelSet;

	double dt = delta_t; // time step.
	double t_total = total_time;
	int num_step = ceil(t_total/dt);
	cout<<"Start TVD Runge-Kutta"<<" \n";
	for (int t = 0; t < num_step; t++)
	{
		//cout<<"TVD Runge-Kutta time step " <<t<<"\n";

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
	cout<<"End TVD Runge-Kutta"<<" \n";
}



void TVD_RK(LEVELSET InputLevelSet, double total_time ,double delta_t, class LEVELSET sign_distance_level_set)
{
	double l0_error;
	double l1_error;
	double error;
	double cfl = 0.45;
	double dt_ij;


	LEVELSET LevelSet1= InputLevelSet;
	LEVELSET LevelSet2=InputLevelSet;
	LEVELSET LevelSetTemp1=InputLevelSet;
	LEVELSET LevelSetTemp2=InputLevelSet;

	double dt = delta_t; // time step.
	double t_total = total_time;
	int num_step = ceil(t_total/dt);

	ostringstream file_name4;
	file_name4 <<  "D:/TVD RK/L0.dat";
	ofstream L0_output4(file_name4.str());
	assert(L0_output4.is_open());

	ostringstream file_name5;
	file_name5 <<  "D:/TVD RK/L1.dat";
	ofstream L1_output(file_name5.str());
	assert(L1_output.is_open());

	cout<<"Start TVD Runge-Kutta"<<" \n";
	for (int t = 0; t < num_step; t++)
	{
		//cout<<"TVD Runge-Kutta time step " <<t<<"\n";

		l1_error = 0;
		l0_error = 0;
		error = 0;
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

		ostringstream file_name3;
		file_name3 <<  "D:/TVD RK/Error "<<t<<".dat";
		ofstream levelset_output3(file_name3.str());
		assert(levelset_output3.is_open());

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
				error =sign_distance_level_set.Psi[i][j] - LevelSet2.Psi[i][j];
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

		L0_output4 <<l0_error  << " ";
		L1_output <<l1_error*(LevelSet1.dx*LevelSet1.dy)<< " ";


		// Update Psi2 . Psi1.
		LevelSet1 = LevelSet2;
	} // end 'for' with t
	L0_output4.close();
	L1_output.close();

	cout<<"End TVD Runge-Kutta"<<" \n";

}

void TVD_RK_subcell(LEVELSET InputLevelSet, double total_time ,double delta_t)
{
	double cfl = 0.45;
	double dt_ij;

	LEVELSET LevelSet1 = InputLevelSet;
	LEVELSET LevelSet2 = InputLevelSet;
	LEVELSET LevelSetTemp1 = InputLevelSet;
	LEVELSET LevelSetTemp2 = InputLevelSet;

	double dt = delta_t; // time step.
	double t_total = total_time;
	int num_step = ceil(t_total/dt);
	cout<<"Start TVD Runge-Kutta"<<" \n";
	for (int t = 0; t < num_step; t++)
	{
		//cout<<"TVD Runge-Kutta time step " <<t<<"\n";

		//ostringstream file_name1;
		//file_name1 <<  "D:/TVD RK/reinitialization set temp "<<t<<".dat";
		//ofstream levelset_output1(file_name1.str());
		//assert(levelset_output1.is_open());

		// Update Psi_temp1.
		for (int i = 0; i <LevelSet1.num_xmesh; i++)
		{
			for (int j = 0; j < LevelSet1.num_ymesh; j++)
			{

				tuple<double,double> dxplus;
				tuple<double,double> dxminus;
				tuple<double,double> dyplus;
				tuple<double,double> dyminus;
				if (i<LevelSet1.num_xmesh-1 && LevelSet1.Psi[i][j]*LevelSet1.Psi[i+1][j]<0)
				{
					dxplus= derivation1x_plus_subcell(i,j,LevelSet1.Psi,LevelSet1.dx,LevelSet1.num_xmesh);
				}
				else
				{
					dxplus = derivation1x_plus(i,j,LevelSet1.Psi,LevelSet1.dx,LevelSet1.num_xmesh);
				}

				if (i>0 && LevelSet1.Psi[i-1][j]*LevelSet1.Psi[i][j]<0)
				{
					dxminus = derivation1x_minus_subcell(i,j,LevelSet1.Psi,LevelSet1.dx,LevelSet1.num_xmesh);
				}
				else
				{
					dxminus = derivation1x_minus(i,j,LevelSet1.Psi,LevelSet1.dx,LevelSet1.num_xmesh);
				}

				if (j<LevelSet1.num_ymesh-1 && LevelSet1.Psi[i][j]*LevelSet1.Psi[i][j+1]<0)
				{
					dyplus = derivation1y_plus_subcell(i,j,LevelSet1.Psi,LevelSet1.dy,LevelSet1.num_ymesh);
				}
				else
				{
					dyplus = derivation1y_plus(i,j,LevelSet1.Psi,LevelSet1.dy,LevelSet1.num_ymesh);
				}

				if (j>0 && LevelSet1.Psi[i][j-1]*LevelSet1.Psi[i][j]<0)
				{
					dyminus = derivation1y_minus_subcell(i,j,LevelSet1.Psi,LevelSet1.dy,LevelSet1.num_ymesh);
				}
				else
				{
					dyminus = derivation1y_minus(i,j,LevelSet1.Psi,LevelSet1.dy,LevelSet1.num_ymesh);
				}


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
		for (int i = 0; i <LevelSetTemp1.num_xmesh; i++)
		{
			for (int j = 0; j < LevelSetTemp1.num_ymesh; j++)
			{

				tuple<double,double> dxplus;
				tuple<double,double> dxminus;
				tuple<double,double> dyplus;
				tuple<double,double> dyminus;
				if (i<LevelSetTemp1.num_xmesh-1 && LevelSetTemp1.Psi[i][j]*LevelSetTemp1.Psi[i+1][j]<0)
				{
					dxplus= derivation1x_plus_subcell(i,j,LevelSetTemp1.Psi,LevelSetTemp1.dx,LevelSetTemp1.num_xmesh);
				}
				else
				{
					dxplus = derivation1x_plus(i,j,LevelSetTemp1.Psi,LevelSetTemp1.dx,LevelSetTemp1.num_xmesh);
				}

				if (i>0 && LevelSetTemp1.Psi[i-1][j]*LevelSetTemp1.Psi[i][j]<0)
				{
					dxminus = derivation1x_minus_subcell(i,j,LevelSetTemp1.Psi,LevelSetTemp1.dx,LevelSetTemp1.num_xmesh);
				}
				else
				{
					dxminus = derivation1x_minus(i,j,LevelSetTemp1.Psi,LevelSetTemp1.dx,LevelSetTemp1.num_xmesh);
				}

				if (j<LevelSetTemp1.num_ymesh-1 && LevelSetTemp1.Psi[i][j]*LevelSetTemp1.Psi[i][j+1]<0)
				{
					dyplus = derivation1y_plus_subcell(i,j,LevelSetTemp1.Psi,LevelSetTemp1.dy,LevelSetTemp1.num_ymesh);
				}
				else
				{
					dyplus = derivation1y_plus(i,j,LevelSetTemp1.Psi,LevelSetTemp1.dy,LevelSetTemp1.num_ymesh);
				}

				if (j>0 && LevelSetTemp1.Psi[i][j-1]*LevelSetTemp1.Psi[i][j]<0)
				{
					dyminus = derivation1y_minus_subcell(i,j,LevelSetTemp1.Psi,LevelSetTemp1.dy,LevelSetTemp1.num_ymesh);
				}
				else
				{
					dyminus = derivation1y_minus(i,j,LevelSetTemp1.Psi,LevelSetTemp1.dy,LevelSetTemp1.num_ymesh);
				}

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
	cout<<"End TVD Runge-Kutta"<<" \n";

}
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
void FE(LEVELSET InputLevelSet, double total_time ,double delta_t)
{
	double cfl = 0.45;
	double dt_ij;


	LEVELSET LevelSet1=InputLevelSet;
	LEVELSET LevelSet2=InputLevelSet;

	double dt = delta_t; // time step.
	double t_total = total_time;
	int num_step = ceil(t_total/dt);

	cout<<"Start  Forward Euler"<<" \n";
	for (int t = 0; t < num_step; t++)
	{
		//cout<<"Forward Euler time step " <<t<<"\n";

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
		LevelSet1 = LevelSet2;
	}
	cout<<"End  Forward Euler"<<" \n";

}




void FE(LEVELSET InputLevelSet, double total_time ,double delta_t, class LEVELSET sign_distance_level_set)
{
	double l0_error;
	double l1_error;
	double error;
	double cfl = 0.45;
	double dt_ij;


	LEVELSET LevelSet1=InputLevelSet;
	LEVELSET LevelSet2=InputLevelSet;


	double dt = delta_t; // time step.
	double t_total = total_time;
	int num_step = ceil(t_total/dt);



	ostringstream file_name4;
	file_name4 <<  "D:/FE/L0.dat";
	ofstream L0_output4(file_name4.str());
	assert(L0_output4.is_open());


	ostringstream file_name5;
	file_name5 <<  "D:/FE/L1.dat";
	ofstream L1_output(file_name5.str());
	assert(L1_output.is_open());



	cout<<"Start  Forward Euler"<<" \n";
	for (int t = 0; t < num_step; t++)
	{
		//cout<<"Forward Euler time step " <<t<<"\n";

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

		L0_output4 << l0_error << " ";
		L1_output << l1_error*(LevelSet1.dx*LevelSet1.dy) << " ";

		// update levelset2 -> levelset1.
		LevelSet1 = LevelSet2;
	}
	cout<<"End  Forward Euler"<<" \n";

	L0_output4.close();
	L1_output.close();

}

void FE_subcell(LEVELSET InputLevelSet, double total_time ,double delta_t)
{
	double cfl = 0.45;
	double dt_ij;


	LEVELSET LevelSet1=InputLevelSet;
	LEVELSET LevelSet2=InputLevelSet;

	double dt = delta_t; // time step.
	double t_total = total_time;
	int num_step = ceil(t_total/dt);

	cout<<"Start  Forward Euler"<<" \n";
	for (int t = 0; t < num_step; t++)
	{
		//cout<<"Forward Euler time step " <<t<<"\n";

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
				tuple<double,double> dxplus;
				tuple<double,double> dxminus;
				tuple<double,double> dyplus;
				tuple<double,double> dyminus;
				if (i<LevelSet1.num_xmesh-1 && LevelSet1.Psi[i][j]*LevelSet1.Psi[i+1][j]<0)
				{
					dxplus= derivation1x_plus_subcell(i,j,LevelSet1.Psi,LevelSet1.dx,LevelSet1.num_xmesh);
				}
				else
				{
					dxplus = derivation1x_plus(i,j,LevelSet1.Psi,LevelSet1.dx,LevelSet1.num_xmesh);
				}

				if (i>0 && LevelSet1.Psi[i-1][j]*LevelSet1.Psi[i][j]<0)
				{
					dxminus = derivation1x_minus_subcell(i,j,LevelSet1.Psi,LevelSet1.dx,LevelSet1.num_xmesh);
				}
				else
				{
					dxminus = derivation1x_minus(i,j,LevelSet1.Psi,LevelSet1.dx,LevelSet1.num_xmesh);
				}

				if (j<LevelSet1.num_ymesh-1 && LevelSet1.Psi[i][j]*LevelSet1.Psi[i][j+1]<0)
				{
					dyplus = derivation1y_plus_subcell(i,j,LevelSet1.Psi,LevelSet1.dy,LevelSet1.num_ymesh);
				}
				else
				{
					dyplus = derivation1y_plus(i,j,LevelSet1.Psi,LevelSet1.dy,LevelSet1.num_ymesh);
				}

				if (j>0 && LevelSet1.Psi[i][j-1]*LevelSet1.Psi[i][j]<0)
				{
					dyminus = derivation1y_minus_subcell(i,j,LevelSet1.Psi,LevelSet1.dy,LevelSet1.num_ymesh);
				}
				else
				{
					dyminus = derivation1y_minus(i,j,LevelSet1.Psi,LevelSet1.dy,LevelSet1.num_ymesh);
				}


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
		LevelSet1 = LevelSet2;
	}
	cout<<"End  Forward Euler"<<" \n";


	


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
void GS(LEVELSET InputLevelSet, double total_time ,double delta_t)
{
	double cfl = 0.45;
	double dt_ij;

	LEVELSET LevelSet1=InputLevelSet;

	double dt = delta_t; // time step.
	double t_total = total_time;
	int num_step = ceil(t_total/dt);

	cout<<"Start  Gauss-Seidel"<<" \n";
	for (int t = 0; t < num_step; t++)
	{
		//cout<<"Gauss-Seidel iteration time step " <<t<<"\n";

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

		// Update Psi1, itself, using four raster-scan visiting.
		switch (raster_visiting)
		{

		case 0 : 
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
					LevelSet1.Psi[i][j] = LevelSet1.Psi[i][j] - dt_ij*sign_func(LevelSet1.Psi[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,LevelSet1.Psi[i][j])-1.0);

					levelset_output1 << LevelSet1.Psi[i][j] << " ";
				} //end 'for' with j

				levelset_output1 << " "<<"\n";
			} //end 'for' with i
			break;




		case 1 : 
			for (int i = 0; i < LevelSet1.num_xmesh; i++)
			{
				for (int j = LevelSet1.num_ymesh-1; j > -1; j--)
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

					//Psi_temp1[i][j] = LevelSet1.Psi[i][j] - dt_ij*sign_func(LevelSet1.Psi[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,LevelSet1.Psi[i][j])-1.0);
					LevelSet1.Psi[i][j] = LevelSet1.Psi[i][j] - dt_ij*sign_func(LevelSet1.Psi[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,LevelSet1.Psi[i][j])-1.0);
					/*levelset_output1 << LevelSet1.Psi[i][j] << " ";
					levelset_output3 << sign_distance_level_set[i][j] - LevelSet1.Psi[i][j] << " ";*/
				} //end 'for' with j

				for (int j = 0; j  < LevelSet1.num_ymesh; j ++)
				{
					levelset_output1 << LevelSet1.Psi[i][j] << " ";
				}

				levelset_output1 << " "<<"\n";
			} //end 'for' with i
			break;




		case 2 : 
			for (int i = LevelSet1.num_xmesh-1; i>-1 ; i--)
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

					//Psi_temp1[i][j] = LevelSet1.Psi[i][j] - dt_ij*sign_func(LevelSet1.Psi[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,LevelSet1.Psi[i][j])-1.0);
					LevelSet1.Psi[i][j] = LevelSet1.Psi[i][j] - dt_ij*sign_func(LevelSet1.Psi[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,LevelSet1.Psi[i][j])-1.0);
					//levelset_output1 << LevelSet1.Psi[i][j] << " ";
					//levelset_output3 << sign_distance_level_set[i][j] - LevelSet1.Psi[i][j] << " ";
				} //end 'for' with j
				//levelset_output1 << " "<<"\n";
				//levelset_output3 << " "<<"\n";
			} //end 'for' with i

			for (int i = 0; i < LevelSet1.num_xmesh; i++)
			{
				for (int j = 0; j < LevelSet1.num_ymesh; j++)
				{
					levelset_output1 << LevelSet1.Psi[i][j] << " ";
				}
				levelset_output1 << " "<<"\n";
			}
			break;



		case 3 : 
			for (int i = LevelSet1.num_xmesh-1; i>-1 ; i--)
			{
				for (int j = LevelSet1.num_ymesh-1; j > -1; j--)
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

					//Psi_temp1[i][j] = LevelSet1.Psi[i][j] - dt_ij*sign_func(LevelSet1.Psi[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,LevelSet1.Psi[i][j])-1.0);
					LevelSet1.Psi[i][j] = LevelSet1.Psi[i][j] - dt_ij*sign_func(LevelSet1.Psi[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,LevelSet1.Psi[i][j])-1.0);

					//levelset_output1 << LevelSet1.Psi[i][j] << " ";
					//levelset_output3 << sign_distance_level_set[i][j] - LevelSet1.Psi[i][j] << " ";
				} //end 'for' with j
				//levelset_output1 << " "<<"\n";
				//levelset_output3 << " "<<"\n";
			} //end 'for' with i

			for (int i = 0; i < LevelSet1.num_xmesh; i++)
			{
				for (int j = 0; j < LevelSet1.num_ymesh; j++)
				{
					levelset_output1 << LevelSet1.Psi[i][j] << " ";
				}
				levelset_output1 << " "<<"\n";
			}
			break;

		} // end 'switch'

		levelset_output1.close();

	}//end 'for' with t
	cout<<"End  Gauss-Seidel"<<" \n";

}




void GS(LEVELSET InputLevelSet, double total_time ,double delta_t, class LEVELSET sign_distance_level_set)
{
	double l0_error;
	double l1_error;
	double error;
	double cfl = 0.45;
	double dt_ij;

	LEVELSET LevelSet1=InputLevelSet;

	double dt = delta_t; // time step.
	double t_total = total_time;
	int num_step = ceil(t_total/dt);

	ostringstream file_name4;
	file_name4 <<  "D:/Gauss Seidel/L0.dat";
	ofstream L0_output4(file_name4.str());
	assert(L0_output4.is_open());


	ostringstream file_name5;
	file_name5 <<  "D:/Gauss Seidel/L1.dat";
	ofstream L1_output(file_name5.str());
	assert(L1_output.is_open());


	cout<<"Start  Gauss-Seidel"<<" \n";
	for (int t = 0; t < num_step; t++)
	{
		//cout<<"Gauss-Seidel iteration time step " <<t<<"\n";

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
					LevelSet1.Psi[i][j] = LevelSet1.Psi[i][j] - dt_ij*sign_func(LevelSet1.Psi[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,LevelSet1.Psi[i][j])-1.0);

					levelset_output1 << LevelSet1.Psi[i][j] << " ";
					error = sign_distance_level_set.Psi[i][j] - LevelSet1.Psi[i][j];
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
			for (int i = 0; i < LevelSet1.num_xmesh; i++)
			{
				for (int j = LevelSet1.num_ymesh-1; j > -1; j--)
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

					//Psi_temp1[i][j] = LevelSet1.Psi[i][j] - dt_ij*sign_func(LevelSet1.Psi[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,LevelSet1.Psi[i][j])-1.0);
					LevelSet1.Psi[i][j] = LevelSet1.Psi[i][j] - dt_ij*sign_func(LevelSet1.Psi[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,LevelSet1.Psi[i][j])-1.0);
					/*levelset_output1 << LevelSet1.Psi[i][j] << " ";
					levelset_output3 << sign_distance_level_set[i][j] - LevelSet1.Psi[i][j] << " ";*/
				} //end 'for' with j

				for (int j = 0; j  < LevelSet1.num_ymesh; j ++)
				{
					levelset_output1 << LevelSet1.Psi[i][j] << " ";
					error = sign_distance_level_set.Psi[i][j] - LevelSet1.Psi[i][j];
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
			for (int i = LevelSet1.num_xmesh-1; i>-1 ; i--)
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

					//Psi_temp1[i][j] = LevelSet1.Psi[i][j] - dt_ij*sign_func(LevelSet1.Psi[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,LevelSet1.Psi[i][j])-1.0);
					LevelSet1.Psi[i][j] = LevelSet1.Psi[i][j] - dt_ij*sign_func(LevelSet1.Psi[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,LevelSet1.Psi[i][j])-1.0);
					//levelset_output1 << LevelSet1.Psi[i][j] << " ";
					//levelset_output3 << sign_distance_level_set[i][j] - LevelSet1.Psi[i][j] << " ";
				} //end 'for' with j
				//levelset_output1 << " "<<"\n";
				//levelset_output3 << " "<<"\n";
			} //end 'for' with i

			for (int i = 0; i < LevelSet1.num_xmesh; i++)
			{
				for (int j = 0; j < LevelSet1.num_ymesh; j++)
				{
					levelset_output1 << LevelSet1.Psi[i][j] << " ";
					error = sign_distance_level_set.Psi[i][j] - LevelSet1.Psi[i][j];
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
			for (int i = LevelSet1.num_xmesh-1; i>-1 ; i--)
			{
				for (int j = LevelSet1.num_ymesh-1; j > -1; j--)
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

					//Psi_temp1[i][j] = LevelSet1.Psi[i][j] - dt_ij*sign_func(LevelSet1.Psi[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,LevelSet1.Psi[i][j])-1.0);
					LevelSet1.Psi[i][j] = LevelSet1.Psi[i][j] - dt_ij*sign_func(LevelSet1.Psi[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,LevelSet1.Psi[i][j])-1.0);

					//levelset_output1 << LevelSet1.Psi[i][j] << " ";
					//levelset_output3 << sign_distance_level_set[i][j] - LevelSet1.Psi[i][j] << " ";
				} //end 'for' with j
				//levelset_output1 << " "<<"\n";
				//levelset_output3 << " "<<"\n";
			} //end 'for' with i

			for (int i = 0; i < LevelSet1.num_xmesh; i++)
			{
				for (int j = 0; j < LevelSet1.num_ymesh; j++)
				{
					levelset_output1 << LevelSet1.Psi[i][j] << " ";
					error = sign_distance_level_set.Psi[i][j] - LevelSet1.Psi[i][j];
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



		L0_output4 << l0_error << " ";

		L1_output << l1_error*(LevelSet1.dx*LevelSet1.dy)<< " ";


	}//end 'for' with t
	cout<<"End  Gauss-Seidel"<<" \n";

	L0_output4.close();
	L1_output.close();

}



void GS_subcell(LEVELSET InputLevelSet, double total_time ,double delta_t)
{
	double cfl = 0.45;
	double dt_ij;

	LEVELSET LevelSet1=InputLevelSet;

	double dt = delta_t; // time step.
	double t_total = total_time;
	int num_step = ceil(t_total/dt);

	cout<<"Start  Gauss-Seidel"<<" \n";
	for (int t = 0; t < num_step; t++)
	{
		//cout<<"Gauss-Seidel iteration time step " <<t<<"\n";

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

		// Update Psi1, itself, using four raster-scan visiting.
		switch (raster_visiting)
		{

		case 0 : 
			for (int i = 0; i < LevelSet1.num_xmesh; i++)
			{
				for (int j = 0; j < LevelSet1.num_ymesh; j++)
				{
					tuple<double,double> dxplus;
					tuple<double,double> dxminus;
					tuple<double,double> dyplus;
					tuple<double,double> dyminus;
					if (i<LevelSet1.num_xmesh-1 && LevelSet1.Psi[i][j]*LevelSet1.Psi[i+1][j]<0)
					{
						dxplus= derivation1x_plus_subcell(i,j,LevelSet1.Psi,LevelSet1.dx,LevelSet1.num_xmesh);
					}
					else
					{
						dxplus = derivation1x_plus(i,j,LevelSet1.Psi,LevelSet1.dx,LevelSet1.num_xmesh);
					}

					if (i>0 && LevelSet1.Psi[i-1][j]*LevelSet1.Psi[i][j]<0)
					{
						dxminus = derivation1x_minus_subcell(i,j,LevelSet1.Psi,LevelSet1.dx,LevelSet1.num_xmesh);
					}
					else
					{
						dxminus = derivation1x_minus(i,j,LevelSet1.Psi,LevelSet1.dx,LevelSet1.num_xmesh);
					}

					if (j<LevelSet1.num_ymesh-1 && LevelSet1.Psi[i][j]*LevelSet1.Psi[i][j+1]<0)
					{
						dyplus = derivation1y_plus_subcell(i,j,LevelSet1.Psi,LevelSet1.dy,LevelSet1.num_ymesh);
					}
					else
					{
						dyplus = derivation1y_plus(i,j,LevelSet1.Psi,LevelSet1.dy,LevelSet1.num_ymesh);
					}

					if (j>0 && LevelSet1.Psi[i][j-1]*LevelSet1.Psi[i][j]<0)
					{
						dyminus = derivation1y_minus_subcell(i,j,LevelSet1.Psi,LevelSet1.dy,LevelSet1.num_ymesh);
					}
					else
					{
						dyminus = derivation1y_minus(i,j,LevelSet1.Psi,LevelSet1.dy,LevelSet1.num_ymesh);
					}

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
					LevelSet1.Psi[i][j] = LevelSet1.Psi[i][j] - dt_ij*sign_func(LevelSet1.Psi[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,LevelSet1.Psi[i][j])-1.0);

					levelset_output1 << LevelSet1.Psi[i][j] << " ";
				} //end 'for' with j

				levelset_output1 << " "<<"\n";
			} //end 'for' with i
			break;




		case 1 : 
			for (int i = 0; i < LevelSet1.num_xmesh; i++)
			{
				for (int j = LevelSet1.num_ymesh-1; j > -1; j--)
				{
					tuple<double,double> dxplus;
					tuple<double,double> dxminus;
					tuple<double,double> dyplus;
					tuple<double,double> dyminus;
					if (i<LevelSet1.num_xmesh-1 && LevelSet1.Psi[i][j]*LevelSet1.Psi[i+1][j]<0)
					{
						dxplus= derivation1x_plus_subcell(i,j,LevelSet1.Psi,LevelSet1.dx,LevelSet1.num_xmesh);
					}
					else
					{
						dxplus = derivation1x_plus(i,j,LevelSet1.Psi,LevelSet1.dx,LevelSet1.num_xmesh);
					}

					if (i>0 && LevelSet1.Psi[i-1][j]*LevelSet1.Psi[i][j]<0)
					{
						dxminus = derivation1x_minus_subcell(i,j,LevelSet1.Psi,LevelSet1.dx,LevelSet1.num_xmesh);
					}
					else
					{
						dxminus = derivation1x_minus(i,j,LevelSet1.Psi,LevelSet1.dx,LevelSet1.num_xmesh);
					}

					if (j<LevelSet1.num_ymesh-1 && LevelSet1.Psi[i][j]*LevelSet1.Psi[i][j+1]<0)
					{
						dyplus = derivation1y_plus_subcell(i,j,LevelSet1.Psi,LevelSet1.dy,LevelSet1.num_ymesh);
					}
					else
					{
						dyplus = derivation1y_plus(i,j,LevelSet1.Psi,LevelSet1.dy,LevelSet1.num_ymesh);
					}

					if (j>0 && LevelSet1.Psi[i][j-1]*LevelSet1.Psi[i][j]<0)
					{
						dyminus = derivation1y_minus_subcell(i,j,LevelSet1.Psi,LevelSet1.dy,LevelSet1.num_ymesh);
					}
					else
					{
						dyminus = derivation1y_minus(i,j,LevelSet1.Psi,LevelSet1.dy,LevelSet1.num_ymesh);
					}

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

					//Psi_temp1[i][j] = LevelSet1.Psi[i][j] - dt_ij*sign_func(LevelSet1.Psi[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,LevelSet1.Psi[i][j])-1.0);
					LevelSet1.Psi[i][j] = LevelSet1.Psi[i][j] - dt_ij*sign_func(LevelSet1.Psi[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,LevelSet1.Psi[i][j])-1.0);
					/*levelset_output1 << LevelSet1.Psi[i][j] << " ";
					levelset_output3 << sign_distance_level_set[i][j] - LevelSet1.Psi[i][j] << " ";*/
				} //end 'for' with j

				for (int j = 0; j  < LevelSet1.num_ymesh; j ++)
				{
					levelset_output1 << LevelSet1.Psi[i][j] << " ";
				}

				levelset_output1 << " "<<"\n";
			} //end 'for' with i
			break;




		case 2 : 
			for (int i = LevelSet1.num_xmesh-1; i>-1 ; i--)
			{
				for (int j = 0; j < LevelSet1.num_ymesh; j++)
				{
					tuple<double,double> dxplus;
					tuple<double,double> dxminus;
					tuple<double,double> dyplus;
					tuple<double,double> dyminus;
					if (i<LevelSet1.num_xmesh-1 && LevelSet1.Psi[i][j]*LevelSet1.Psi[i+1][j]<0)
					{
						dxplus= derivation1x_plus_subcell(i,j,LevelSet1.Psi,LevelSet1.dx,LevelSet1.num_xmesh);
					}
					else
					{
						dxplus = derivation1x_plus(i,j,LevelSet1.Psi,LevelSet1.dx,LevelSet1.num_xmesh);
					}

					if (i>0 && LevelSet1.Psi[i-1][j]*LevelSet1.Psi[i][j]<0)
					{
						dxminus = derivation1x_minus_subcell(i,j,LevelSet1.Psi,LevelSet1.dx,LevelSet1.num_xmesh);
					}
					else
					{
						dxminus = derivation1x_minus(i,j,LevelSet1.Psi,LevelSet1.dx,LevelSet1.num_xmesh);
					}

					if (j<LevelSet1.num_ymesh-1 && LevelSet1.Psi[i][j]*LevelSet1.Psi[i][j+1]<0)
					{
						dyplus = derivation1y_plus_subcell(i,j,LevelSet1.Psi,LevelSet1.dy,LevelSet1.num_ymesh);
					}
					else
					{
						dyplus = derivation1y_plus(i,j,LevelSet1.Psi,LevelSet1.dy,LevelSet1.num_ymesh);
					}

					if (j>0 && LevelSet1.Psi[i][j-1]*LevelSet1.Psi[i][j]<0)
					{
						dyminus = derivation1y_minus_subcell(i,j,LevelSet1.Psi,LevelSet1.dy,LevelSet1.num_ymesh);
					}
					else
					{
						dyminus = derivation1y_minus(i,j,LevelSet1.Psi,LevelSet1.dy,LevelSet1.num_ymesh);
					}

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

					//Psi_temp1[i][j] = LevelSet1.Psi[i][j] - dt_ij*sign_func(LevelSet1.Psi[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,LevelSet1.Psi[i][j])-1.0);
					LevelSet1.Psi[i][j] = LevelSet1.Psi[i][j] - dt_ij*sign_func(LevelSet1.Psi[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,LevelSet1.Psi[i][j])-1.0);
					//levelset_output1 << LevelSet1.Psi[i][j] << " ";
					//levelset_output3 << sign_distance_level_set[i][j] - LevelSet1.Psi[i][j] << " ";
				} //end 'for' with j
				//levelset_output1 << " "<<"\n";
				//levelset_output3 << " "<<"\n";
			} //end 'for' with i

			for (int i = 0; i < LevelSet1.num_xmesh; i++)
			{
				for (int j = 0; j < LevelSet1.num_ymesh; j++)
				{
					levelset_output1 << LevelSet1.Psi[i][j] << " ";
				}
				levelset_output1 << " "<<"\n";
			}
			break;



		case 3 : 
			for (int i = LevelSet1.num_xmesh-1; i>-1 ; i--)
			{
				for (int j = LevelSet1.num_ymesh-1; j > -1; j--)
				{

					tuple<double,double> dxplus;
					tuple<double,double> dxminus;
					tuple<double,double> dyplus;
					tuple<double,double> dyminus;
					if (i<LevelSet1.num_xmesh-1 && LevelSet1.Psi[i][j]*LevelSet1.Psi[i+1][j]<0)
					{
						dxplus= derivation1x_plus_subcell(i,j,LevelSet1.Psi,LevelSet1.dx,LevelSet1.num_xmesh);
					}
					else
					{
						dxplus = derivation1x_plus(i,j,LevelSet1.Psi,LevelSet1.dx,LevelSet1.num_xmesh);
					}

					if (i>0 && LevelSet1.Psi[i-1][j]*LevelSet1.Psi[i][j]<0)
					{
						dxminus = derivation1x_minus_subcell(i,j,LevelSet1.Psi,LevelSet1.dx,LevelSet1.num_xmesh);
					}
					else
					{
						dxminus = derivation1x_minus(i,j,LevelSet1.Psi,LevelSet1.dx,LevelSet1.num_xmesh);
					}

					if (j<LevelSet1.num_ymesh-1 && LevelSet1.Psi[i][j]*LevelSet1.Psi[i][j+1]<0)
					{
						dyplus = derivation1y_plus_subcell(i,j,LevelSet1.Psi,LevelSet1.dy,LevelSet1.num_ymesh);
					}
					else
					{
						dyplus = derivation1y_plus(i,j,LevelSet1.Psi,LevelSet1.dy,LevelSet1.num_ymesh);
					}

					if (j>0 && LevelSet1.Psi[i][j-1]*LevelSet1.Psi[i][j]<0)
					{
						dyminus = derivation1y_minus_subcell(i,j,LevelSet1.Psi,LevelSet1.dy,LevelSet1.num_ymesh);
					}
					else
					{
						dyminus = derivation1y_minus(i,j,LevelSet1.Psi,LevelSet1.dy,LevelSet1.num_ymesh);
					}

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

					//Psi_temp1[i][j] = LevelSet1.Psi[i][j] - dt_ij*sign_func(LevelSet1.Psi[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,LevelSet1.Psi[i][j])-1.0);
					LevelSet1.Psi[i][j] = LevelSet1.Psi[i][j] - dt_ij*sign_func(LevelSet1.Psi[i][j]) * (HG(derv_xplus,derv_xminus,derv_yplus,derv_yminus,LevelSet1.Psi[i][j])-1.0);

					//levelset_output1 << LevelSet1.Psi[i][j] << " ";
					//levelset_output3 << sign_distance_level_set[i][j] - LevelSet1.Psi[i][j] << " ";
				} //end 'for' with j
				//levelset_output1 << " "<<"\n";
				//levelset_output3 << " "<<"\n";
			} //end 'for' with i

			for (int i = 0; i < LevelSet1.num_xmesh; i++)
			{
				for (int j = 0; j < LevelSet1.num_ymesh; j++)
				{
					levelset_output1 << LevelSet1.Psi[i][j] << " ";
				}
				levelset_output1 << " "<<"\n";
			}
			break;

		} // end 'switch'

		levelset_output1.close();

	}//end 'for' with t
	cout<<"End  Gauss-Seidel"<<" \n";

}

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////
////                       Gauss-Seidel iteration
////
////                             END
////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////


