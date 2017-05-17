#include <cmath>
#include "tbb/tbb.h"
#include "Multigrid3D.h"
#include "array.h"
#include "fluid_buffer3D.h"
#include "Smoke_solver3D.h"
#include "SDFGen.h"
#include <iostream>

SDFGenerator fan_blade;
float g_time=0.0f;
float g_dt = 0.01f;
float g_u_in;
#include <ctime>
using namespace ivock;
double diffclock(clock_t clock1,clock_t clock2)
{
	double diffticks=clock1-clock2;
	double diffms=(diffticks)/(CLOCKS_PER_SEC/1000);
	return diffms/1000.0;
}
float sphere_phi(const Vec3f &pos, const Vec3f &center, float r)
{

	
	return dist(pos, center) - r;
}
Vec3f sphere_vel(const Vec3f &pos, const Vec3f &center, float r)
{
	Vec3f vel(0.0f);
	
	return vel;
}
float boundary_phi(const Vec3f & pos)
{

	return sphere_phi(pos, Vec3f(3.0,5.0,5.0),1.0);
}
Vec3f boundary_vel(const Vec3f & pos)
{
	return sphere_vel(pos, Vec3f(3.0,5.0,5.0),1.0);
}
int main(int argc, char** argv) {


	if(argc<7)
	{
		printf("please specify output path, dt, uniform velocity, nx, ny, nz, for example: Smoke_solver.exe C:/ 0.03 5.0 128 64 64\n");
		return 0;
	}
	else{
		int advection_type;
		float vort_confine_str=0;
		int nx, ny, nz;
		float ratio;
		char file_path[256];
		int n=sprintf(file_path,"%s",argv[1]);
		sscanf(argv[2], "%f", &g_dt);
		sscanf(argv[3], "%f", &g_u_in);
		sscanf(argv[4], "%d", &nx);
		sscanf(argv[5], "%d", &ny);
		sscanf(argv[6], "%d", &nz);
		ratio = 0.25f; // please modify this value to modify subdomain h 

		//!!!
		//int nx=128,ny=64,nz=64;
		float L = 20;
		float g_h = L/(float)nx;
		SmokeSolver3D g_smokeSolver;
		g_smokeSolver.init(nx,ny,nz,L, ratio);
		g_smokeSolver.initSubdomain(Vec3f(1.5, 3.5, 3.5),Vec3f(4.5,6.5,6.5),g_smokeSolver._hx*ratio);



		// an example here to initialize sdf from an obj, not user friendly at this moment.
		//80 40 40  dt 0.03
		//fan_blade.dx=(0.25f*g_smokeSolver._hx);
		//fan_blade.objToSDF("disk.obj",Vec3f(-2.5,-1.0,-2.5),Vec3f(2.5,1.0,2.5));

		//printf("sdf done\n");



		buffer3Dc bc_des;
		bc_des.init(nx,ny,nz);
		
		bc_des.setZero();
		for (int k=0;k<nz;k++)for(int j=0;j<ny;j++)for(int i=0;i<nx;i++)
		{
			//0:fluid;1:air;2:solid
			if(i<=1) bc_des(i,j,k) = 1;
			if(j<=1) bc_des(i,j,k) = 1;
			if(k<=1) bc_des(i,j,k) = 1;

			if(i>=nx-2) bc_des(i,j,k) = 1;
			if(j>=ny-2) bc_des(i,j,k) = 1;
			if(k>=nz-2) bc_des(i,j,k) = 1;

		}
		g_smokeSolver.set_boundary(bc_des);//in flow condition
		g_smokeSolver.set_boundary_phi(boundary_phi);
		g_smokeSolver.set_boundary_vel(boundary_vel);
		
		g_smokeSolver.setEmitter(Vec3f(5,1,5), 1.0, 5000);//should do it in the way you want, or add multiple emitters



		clock_t start = clock();
		for (int frame = 0; frame<600;frame++)
		{
			{

				g_smokeSolver.set_viscosity(vort_confine_str);
				g_smokeSolver.set_U_in(g_u_in);
				g_smokeSolver.set_boundary_phi(boundary_phi);
				g_smokeSolver.set_boundary_vel(boundary_vel);
				g_smokeSolver.time_step(g_dt,0,frame);
				g_time += g_dt;
			}
			
			g_smokeSolver.write_tracers(file_path,frame);
			printf("frame %d done\n",frame);
		}
		clock_t end = clock();
		cout << diffclock(end,start)/6.0<<endl;




		return 0;
	}
	
}
