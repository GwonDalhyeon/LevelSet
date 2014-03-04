
#include "reinitialzation_scheme.h"

using namespace std;


int main()
{
	
	int num_mesh, num_xmesh, num_ymesh;  // Decide the number of mesh , not size. And times.
	double side_length, side_xlength, side_ylength; // length foreach axis. 
	// codinate vectors for each axis, and TImes.

	side_length = 4.0;
	num_mesh = 151;

	num_xmesh = num_mesh;
	num_ymesh = num_mesh;
	side_xlength = side_length;
	side_ylength = side_length;


	// Define initial level set which is not a signed distance function

	class LEVELSET initial_level_set(side_xlength, side_ylength, num_xmesh, num_ymesh);
	class LEVELSET sign_distance_level_set(side_xlength, side_ylength, num_xmesh,num_ymesh);

	ofstream initial_levelset_output("initial level set.dat"); // Write data file.
	assert(initial_levelset_output.is_open());
	ofstream distance_levelset_output("sign distance level set.dat");
	assert(distance_levelset_output.is_open());
	double a=0.7;
	double r=1.0;
	for (int i = 0; i < initial_level_set.num_xmesh; i++) // Loop for X Axis
	{
		for (int j = 0; j < initial_level_set.num_ymesh; j++) // Loop for Y axis
		{
			//// A circle with center at the origen and radius 1.
			initial_level_set.Psi[i][j] = (sqrt(initial_level_set.X[i]*initial_level_set.X[i]+initial_level_set.Y[j]*initial_level_set.Y[j])-1.0)
				*((initial_level_set.X[i]-1.0)*(initial_level_set.X[i]-1.0)+(initial_level_set.Y[j]-1.0)*(initial_level_set.Y[j]-1.0)+0.1);
			sign_distance_level_set.Psi[i][j] = sqrt(sign_distance_level_set.X[i]*sign_distance_level_set.X[i]+sign_distance_level_set.Y[j]*sign_distance_level_set.Y[j])-1.0;


			//// Two circles of radius r are placed at (+-a,0) on the plane. Let 0<a<r, sh that the two circles intersect each other.
			/*double temp1 = (a-initial_level_set.X[i])/sqrt((a-initial_level_set.X[i])*(a-initial_level_set.X[i])+initial_level_set.Y[j]*initial_level_set.Y[j]);
			double temp2 = (a+initial_level_set.X[i])/sqrt((a+initial_level_set.X[i])*(a+initial_level_set.X[i])+initial_level_set.Y[j]*initial_level_set.Y[j]);
			if (temp1>= a/r && temp2 >=a/r)
			{
				double temp3 = min(  initial_level_set.X[i]*initial_level_set.X[i] + (initial_level_set.Y[j]+sqrt(r*r-a*a))*(initial_level_set.Y[j]+sqrt(r*r-a*a)),
					initial_level_set. X[i]*initial_level_set.X[i] + (initial_level_set.Y[j]-sqrt(r*r-a*a))*(initial_level_set.Y[j]-sqrt(r*r-a*a)));
				initial_level_set.Psi[i][j] =sqrt(temp3) * ((initial_level_set.X[i]-1.0)*(initial_level_set.X[i]-1.0)+(initial_level_set.Y[j]-1.0)*(initial_level_set.Y[j]-1.0)+0.1);
				sign_distance_level_set.Psi[i][j] =sqrt(temp3);
			}
			else
			{
				double temp3 = min(  sqrt((initial_level_set.X[i]+a)*(initial_level_set.X[i]+a) + initial_level_set.Y[j]*initial_level_set.Y[j]) - r ,
					sqrt((initial_level_set.X[i]-a)*(initial_level_set.X[i]-a) + initial_level_set.Y[j]*initial_level_set.Y[j]) - r);
				initial_level_set.Psi[i][j] =temp3 * ((initial_level_set.X[i]-1.0)*(initial_level_set.X[i]-1.0)+(initial_level_set.Y[j]-1.0)*(initial_level_set.Y[j]-1.0)+0.1); 
				sign_distance_level_set.Psi[i][j] =temp3 ;
			}*/


			initial_levelset_output << initial_level_set.Psi[i][j] << " ";
			distance_levelset_output << sign_distance_level_set.Psi[i][j] << " ";
		}
		initial_levelset_output << " "<<"\n";
		distance_levelset_output << " "<<"\n";
	}
	initial_levelset_output.close();
	distance_levelset_output.close();


	//////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	//                       reinitialization
	//
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	
	double dt, t_total;
	dt = 0.0005;
	t_total = 0.2;
	//t_total = 0.001;

	////Without error
	//TVD_RK(initial_level_set,t_total,dt);
	//FE(initial_level_set,t_total,dt);
	//GS(initial_level_set,t_total,dt);

	//// With error
	TVD_RK(initial_level_set,t_total,dt, sign_distance_level_set);
	FE(initial_level_set,t_total,dt, sign_distance_level_set);
	GS(initial_level_set,t_total,dt, sign_distance_level_set);
}