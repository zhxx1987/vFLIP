#ifndef _smoke_solver3D_
#define _smoke_solver3D_
#include "array.h"
#include "tbb/tbb.h"
#include <vector>
#include "Multigrid3D.h"
#include "fluid_buffer3D.h"
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <cstdio>
#include "fluid_particle.h"
#include "vec.h"
#include "pcg_solver.h"
#include "array3.h"
#include "sub_domain.h"
#include "GeometricLevelGen.h"

#include "AlgebraicMultigrid.h"
using namespace std;
namespace ivock{

	class SmokeSolver3D
	{
	public:
		SmokeSolver3D(){}
		~SmokeSolver3D(){Finalize();}
		SubDomainSolver sub_domain;
		Vec3f emitter_pos;
		float emitter_r;
		uint  emitter_n;
		vector<Vec4f> tracers;
		Vec3f origin;
		float frand(float a, float b)
		{
			return a + (b-a)*((float)(rand()%RAND_MAX)/(float)RAND_MAX);
		}
		void setEmitter(Vec3f & pos, float r, uint n)
		{
			emitter_pos = pos;
			emitter_r = r;
			emitter_n = n;
		}
		void emit_tracers()
		{
			cout<<"emitting tracers:"<<endl;
			vector<Vec4f> tracers_temp;
			tracers_temp.resize(0);
			for (uint i=0;i<tracers.size();i++)
			{
				if (tracers[i][0]>2*_hx/* + origin[0]*/ && 
					tracers[i][1]>2*_hx/* + origin[1]*/ &&
					tracers[i][2]>2*_hx/* + origin[2]*/ &&
					tracers[i][0]<_lx-2*_hx/* + origin[0]*/ &&
					tracers[i][1]<_ly-2*_hx/* + origin[1]*/ &&
					tracers[i][2]<_lz-2*_hx/* + origin[2]*/ &&
					tracers[i][3]>0.01)
				{
					Vec3f pos(tracers[i][0], tracers[i][1],tracers[i][2]);
					//if(get_solid_phi(pos)>0.1f*_hx)
						tracers_temp.push_back(tracers[i]);
				}
			}
			tracers.swap(tracers_temp);
			uint num = 0;
			while(num<emitter_n)
			{
				float r = emitter_r;
				float x = frand(-r-_hx, r+_hx);
				float y = frand(-r-_hx, r+_hx);
				float z = frand(-r-_hx, r+_hx);

				if (x*x + z*z + y*y <= r*r)
				{
					tracers.push_back(Vec4f(emitter_pos[0]+x,
						emitter_pos[1]+y,
						emitter_pos[2]+z,
						1.0)
						);
					num++;
				}
			}

			for (int k=0;k<sub_domain._nz;k++)for(int j=0;j<sub_domain._ny;j++)for(int i=0;i<sub_domain._nx;i++)
			{
				for (int p=0;p<4;p++)
				{
					Vec3f pos = sub_domain.bmin + Vec3f(i,j,k)*sub_domain._hx + Vec3f(frand(0,1),frand(0,1),frand(0,1))*sub_domain._hx;
					if (get_solid_phi(pos)>0&&get_solid_phi(pos)<sub_domain._hx)
					{
						tracers.push_back(Vec4f(pos[0],pos[1],pos[2],1.0));
					}
				}
			}


			//cout<<"emitting tracers done:"<<tracers.size()<<" tracers"<<endl;
		}
		Vec3f traceRK1(Vec3f & pos, float dt);
		void advect_tracers(float dt)
		{
			tbb::parallel_for((size_t)0,
				(size_t)tracers.size(),
				(size_t)1,
				[&](size_t i)
			{
				Vec3f pos = Vec3f(tracers[i][0],tracers[i][1],tracers[i][2]);
				pos = traceRK3(pos,dt);
				///*Vec3f vel = get_velocity(pos);
				//Vec3f mpos = pos + 0.5f * dt * vel;
				//pos = pos + dt*get_velocity(mpos);*/
				float phi_value = get_solid_phi(pos);
				if (phi_value<=0)
				{
					Vec3f normal = get_solid_normal(pos);//_nodal_solid_phi.sample_grad(pos[0],pos[1],pos[2]);
					//normalize(normal);
					pos -= phi_value*normal;

				}
				float d = tracers[i][3];
				d = d / (1.0 + _dens_decay * dt);
				tracers[i] = Vec4f(pos[0],pos[1],pos[2], d);
			});
		}
		void write_tracers(char * file_path, int frame)
		{
			char file_name[256];
			sprintf(file_name,"%s/Particle_data%04d.bin", file_path,frame);
			float* data;
			data = new float[4*tracers.size()];

			tbb::parallel_for((size_t)0,
				(size_t)tracers.size(),
				(size_t)1,
				[&](size_t i)
			{
				data[i*4+0] = tracers[i][0]/_lx - 0.5;
				data[i*4+1] = tracers[i][1]/_lx;
				data[i*4+2] = tracers[i][2]/_lx - 0.5;
				data[i*4+3] = tracers[i][3];
			});


			
			FILE *data_file = fopen(file_name,"wb");
			fwrite(data,sizeof(float)*4,tracers.size(),data_file);
			fclose(data_file);
			delete []data;







			


		}
		uint _nx, _ny, _nz, _n;
		float _hx, _hy, _hz;
		float _lx, _ly, _lz;
		float _temp_decay, _dens_decay;
		float _alpha, _beta;
		float _smoke_heat, _smoke_dens, _smoke_fuel;
		float _smoke_radius;
		float _U_in;
		Vec3f _smoke_center;
		//buffers:
		buffer3Df _un, _vn, _wn, _utemp,_vtemp,_wtemp,_unp1,_vnp1,_wnp1;
		buffer3Df _u_weights, _v_weights, _w_weights;
		buffer3Df _u_solid, _v_solid, _w_solid;
		buffer3Df _nodal_solid_phi;
		Array3f u_extrap,v_extrap,w_extrap;
		Array3c u_valid, v_valid, w_valid;
		buffer3Dc _b_desc, _h_desc;;
		float _cfl;
		float _vort_confine_coef;
		float nu;
		std::vector< std::vector<int> > particle_hash;
		buffer3Df _rho, _Temperature, _temp;
		//Solver data
		PCGSolver<double> solver;
		SparseMatrixd matrix;
		std::vector<double> rhs;
		std::vector<double> pressure;


		void setSmokeDens(float _v) {_smoke_dens = _v;}
		void setSmokeHeat(float _v) {_smoke_heat = _v;}
		void setSmokeHeatDecay(float _v) {_temp_decay = _v;}
		void setSmokeDensDecay(float _v) {_dens_decay = _v;}
		void setSmokeRadius(float _v) {_smoke_radius = _v;}
		void setSmokeCenter(Vec3f &_v) {_smoke_center = _v;}
		void setDensBuoy(float _v){_alpha = _v;}
		void setTempBuoy(float _v){_beta = _v;}





		float PSE_kernel(Vec3f &xi, Vec3f & xj, float tau)
		{
			//return exp(-dist2(xi,xj)*0.25f/tau);
			float frac = sqrt(4.0*3.1415926*tau);
			return 1.0f/frac*frac*frac*exp(-dist2(xi,xj)*0.25/tau);
		}
		void set_viscosity(float _nu) {nu = _nu;}
		void seed_ghost();
		void diffuse(float dt);
		void grid_diffuse(float dt);
		std::vector<particle> FLIP;
		void getDelta(buffer3Df &field_old, buffer3Df &field_new, buffer3Df &fieldtemp);
		void compute_delta_u();
		float H(float r)
		{
			float res = 0;
			if(r>=-1 && r<0) res = 1+r;
			if(r>=0 && r<1) res = 1-r;
			return res;
		}
		float compute_weight(float gx,float gy,float gz,
			float px,float py,float pz)
		{
			//k(x,y,z) = H(dx/hx)H(dy/hx)H(dz/hx)
			//H(r) = 1-r 0<=r<1  1+r -1<=r<0 0 else;
			float dx = px - gx;
			float dy = py - gy;
			float dz = pz - gz;
			return H(dx/_hx)*H(dy/_hx)*H(dz/_hx);
		}
		void particle_from_grid();
		void particle_to_grid(buffer3Df & field,buffer3Df &coef, int component);
		void advect_FLIP(float dt);
		void hash_FLIP();
		void reorder_FLIP();
		void seedHalo();
		Vec3f get_solid_normal(Vec3f & pos);
		Vec3f get_solid_vel(Vec3f & pos);
		inline Vec3f get_velocity(Vec3f & pos);
		inline Vec3f SmokeSolver3D::get_velocity_self(Vec3f & pos);
		float get_solid_phi(Vec3f & pos);
		Vec3f get_velocity_diff(Vec3f & pos);
		Vec3f traceRK3(Vec3f & pos, float dt);
		Vec3f traceRK2(Vec3f & pos, float dt);
		Vec3f trace(float dt, Vec3f & pos);
		void getCFL(float dt );
		void extrapolate(Array3f & grid, buffer3Df & u, Array3c & valid);
		void constrain_velocity();
		void set_vort_confine(float str) { _vort_confine_coef = str; }
		void set_U_in(float val) {_U_in = val;}

		void Finalize()
		{
			_un.free();
			_utemp.free();
			_unp1.free();
			_vn.free();
			_vtemp.free();
			_vnp1.free();
			_wn.free();
			_wtemp.free();
			_wnp1.free();

			

			_b_desc.free();
			_h_desc.free();
			_nodal_solid_phi.free();
			_u_weights.free();
			_v_weights.free();
			_w_weights.free();
			_u_solid.free();
			_v_solid.free();
			_w_solid.free();
			_rho.free();
			_Temperature.free();
			_temp.free();

		}
		//init
		void init(uint nx, uint ny, uint nz, float L, float ratio)
		{
			origin = Vec3f(0);
			_nx=nx;
			_ny=ny;
			_nz=nz;
			_hx = L/(float)_nx;
			_hz = _hy = _hx;
			_lx = L; _ly = _hy*(float)ny; _lz = _hz*(float)nz;
			_temp_decay = 0;
			_alpha = _beta = 0;
			_smoke_dens = _smoke_heat = 0;

			//FLIP_Solver.init(_nx,_ny,_nz,_hx);

			_un.init(_nx+1,_ny,_nz,_hx,0.5,0,0);
			_utemp.init(_nx+1,_ny,_nz,_hx,0.5,0,0);
			_unp1.init(_nx+1,_ny,_nz,_hx,0.5,0,0);
			_vn.init(_nx,_ny+1,_nz,_hx,0.0,0.5,0.0);
			_vtemp.init(_nx,_ny+1,_nz,_hx,0.0,0.5,0.0);
			_vnp1.init(_nx,_ny+1,_nz,_hx,0.0,0.5,0.0);
			_wn.init(_nx,_ny,_nz+1,_hx,0.0,0.0,0.5);
			_wtemp.init(_nx,_ny,_nz+1,_hx,0.0,0.0,0.5);
			_wnp1.init(_nx,_ny,_nz+1,_hx,0.0,0.0,0.5);
			_u_solid.init(_nx+1,_ny,_nz,_hx,0.5,0,0);
			_v_solid.init(_nx,_ny+1,_nz,_hx,0.0,0.5,0.0);
			_w_solid.init(_nx,_ny,_nz+1,_hx,0.0,0.0,0.5);
			
			
			_rho.init(_nx*2, _ny*2, _nz*2, 0.5*_hx, 0.0,0.0,0.0);
			_Temperature.init(_nx*2, _ny*2, _nz*2, 0.5*_hx, 0.0,0.0,0.0);
			_temp.init(_nx*2, _ny*2, _nz*2, 0.5*_hx, 0.0,0.0,0.0);
			
			u_valid.resize(_nx+1,_ny,_nz);
			v_valid.resize(_nx,_ny+1,_nz);
			w_valid.resize(_nx,_ny,_nz+1);
			u_extrap.resize(_nx+1,_ny,_nz);
			v_extrap.resize(_nx,_ny+1,_nz);
			w_extrap.resize(_nx,_ny,_nz+1);

			_u_weights.init(_nx+1,_ny,_nz);
			_v_weights.init(_nx,_ny+1,_nz);
			_w_weights.init(_nx,_ny,_nz+1);

			_nodal_solid_phi.init(_nx+1,_ny+1,_nz+1,_hx,0.5,0.5,0.5);
			
			

			_b_desc.init(_nx,_ny,_nz);
			
			_h_desc.init(_nx,_ny,_nz);
			
			hash_FLIP();
			reorder_FLIP();
			printf("after initialize: %d\n",FLIP.size());
			
		}
		void initSubdomain(Vec3f &submin, Vec3f &submax, float h)
		{
			sub_domain.init(submin, submax, h);
		}
		void setSmoke(double temp_decay, double alpha,double beta,double smoke_heat,double smoke_dens){
			_temp_decay = temp_decay;
			_alpha = alpha;
			_beta = beta;
			_smoke_heat = smoke_heat;
			_smoke_dens = smoke_dens;
		}
		void advect(float dt);
		void advect_field(float dt, buffer3Df & field, buffer3Df &field_new);
		void advect_field_cubic(float dt, buffer3Df & field, buffer3Df &field_new);
		void advect_field_cubic_clamp(float dt, buffer3Df & field, buffer3Df &field_new);

		void gen_heat(float dt);
		void add_smokeforce(float dt);
		void set_boundary(buffer3Dc & b_desc) {_b_desc.copy(b_desc);}
		void set_heat(buffer3Dc & h_desc)
		{
			_h_desc.copy(h_desc);

		}
		void time_step(float dt,int adv_type, int frame);
		void pcg_projection(float dt);
		//void output(uint nx, uint ny, uint nz, int frame, char* file_path);
		void compute_face_weights();
		void set_boundary_phi(float (*phi)(const Vec3f&));
		void set_boundary_vel(Vec3f (*vel)(const Vec3f&));

	};


}



#endif