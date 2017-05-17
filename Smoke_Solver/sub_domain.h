#ifndef _SubDomainSolver_
#define _SubDomainSolver_
#include "array.h"
#include "tbb/tbb.h"
#include "Multigrid3D.h"
#include "fluid_buffer3D.h"
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <cstdio>
#include "vec.h"
#include "pcg_solver.h"
#include "array3.h"
#include "fluid_particle.h"
#include "GeometricLevelGen.h"

#include "AlgebraicMultigrid.h"
using namespace std;
namespace ivock{

	class SubDomainSolver
	{
	public:
		SubDomainSolver(){}
		~SubDomainSolver(){Finalize();}
		Vec3f emitter_pos;
		float emitter_r;
		uint  emitter_n;
		vector<Vec4f> tracers;
		Vec3f bmin, bmax;
		float frand(float a, float b)
		{
			return a + (b-a)*((float)(rand()%RAND_MAX)/(float)RAND_MAX);
		}
		
		uint _nx, _ny, _nz, _n;
		float _hx, _hy, _hz;
		float _lx, _ly, _lz;
		float _temp_decay;
		float _alpha, _beta;
		float _smoke_heat, _smoke_dens, _smoke_fuel;
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

		//Solver data
		PCGSolver<double> solver;
		SparseMatrixd matrix;
		std::vector<double> rhs;
		std::vector<double> pressure;

		std::vector<std::vector<int>> particle_hash;
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
			float px,float py,float pz, float h)
		{
			//k(x,y,z) = H(dx/hx)H(dy/hx)H(dz/hx)
			//H(r) = 1-r 0<=r<1  1+r -1<=r<0 0 else;
			float dx = px - gx;
			float dy = py - gy;
			float dz = pz - gz;
			return H(dx/h)*H(dy/h)*H(dz/h);
		}
		void particle_from_grid();
		void extrapolate()
		{
			extrapolate(u_extrap,_un,u_valid);
			extrapolate(v_extrap,_vn,v_valid);
			extrapolate(w_extrap,_wn,w_valid);
		}
		void determineDist(float & h, int & l, int i, int j, int k)
		{
			int dist = min(min(min(i, (int)_nx - i), min(j, (int)_ny - j)), min(k, (int)_nz - k));
			if(dist<=4)
			{
				l = 4;
				h = 4.0*_hx;
			}
			else
			{
			
				l = 1;
				h = _hx;
			}
		}
		void particle_to_grid(buffer3Df & field,buffer3Df &coef, int component);
		void particle_to_grid(){
			_utemp.setZero();
			_vtemp.setZero();
			_wtemp.setZero();
			particle_to_grid(_un,_unp1,0);
			particle_to_grid(_vn,_vnp1,1);
			particle_to_grid(_wn,_wnp1,2);
			_utemp.copy(_un);
			_vtemp.copy(_vn);
			_wtemp.copy(_wn);
		}
		void advect_FLIP(float dt);
		void hash_FLIP();
		void reorder_FLIP();

		Vec3f get_solid_vel(Vec3f & pos);
		inline Vec3f get_velocity(Vec3f & pos);
		float get_solid_phi(Vec3f & pos);
		Vec3f get_solid_normal(Vec3f &pos);
		Vec3f get_velocity_diff(Vec3f & pos);
		Vec3f traceRK3(Vec3f & pos, float dt);
		Vec3f traceRK2(Vec3f & pos, float dt);
		Vec3f trace(float dt, Vec3f & pos);
		void getCFL(float dt );
		void extrapolate(Array3f & grid, buffer3Df & u, Array3c & valid);
		void constrain_velocity();
		void set_vort_confine(float str) { _vort_confine_coef = str; }


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
		}


		void init(Vec3f &_bmin, Vec3f &_bmax, float dx)
		{
			bmin = _bmin; bmax = _bmax;
			float L = bmax[0] - bmin[0];
			int nx = floor(L/dx);
			int ny = floor((bmax[1]-bmin[1])/dx);
			int nz = floor((bmax[2]-bmin[2])/dx);
			printf("%d,%d,%d\n",nx,ny,nz);
			L = (float)nx * dx;
			bmax = bmin + Vec3f(nx,ny,nz)*dx;
			init(nx,ny,nz,L);
			printf("init done\n");
			hash_FLIP();
			printf("hash done\n");
			reorder_FLIP();
			printf("reorder FLIP done\n");
		}
		bool is_in(Vec3f &pos)
		{
			if (pos[0]>bmin[0]+2.0*_hx&&pos[0]<bmax[0]-2.0*_hx
				&&pos[1]>bmin[1]+2.0*_hx&&pos[1]<bmax[1]-2.0*_hx
				&&pos[2]>bmin[2]+2.0*_hx&&pos[2]<bmax[2]-2.0*_hx)
			{
				return true;
			}
			else
			{
				return false;
			}
		}
		bool is_in2(Vec3f &pos)
		{
			if (pos[0]>bmin[0]+5.0*_hx&&pos[0]<bmax[0]-5.0*_hx
				&&pos[1]>bmin[1]+5.0*_hx&&pos[1]<bmax[1]-5.0*_hx
				&&pos[2]>bmin[2]+5.0*_hx&&pos[2]<bmax[2]-5.0*_hx)
			{
				return true;
			}
			else
			{
				return false;
			}
		}
		bool is_halo(Vec3f & pos, float dx)
		{
			if (pos[0]>bmin[0]-dx&&pos[0]<bmax[0]+dx
				&&pos[1]>bmin[1]-dx&&pos[1]<bmax[1]+dx
				&&pos[2]>bmin[2]-dx&&pos[2]<bmax[2]+dx)
			{
				return true;
			}
			else
			{
				return false;
			}
		}
		
		//init
		void init(uint nx, uint ny, uint nz, float L)
		{
			
			_nx=nx;
			_ny=ny;
			_nz=nz;
			_hx = L/(float)_nx;
			_hz = _hy = _hx;
			_lx = L; _ly = _hy*(float)ny; _lz = _hz*(float)nz;
			_temp_decay = 0;
			_alpha = _beta = 0;
			_smoke_dens = _smoke_heat = 0;

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

		
		void set_boundary(buffer3Dc & b_desc) {_b_desc.copy(b_desc);}
		void set_heat(buffer3Dc & h_desc)
		{
			_h_desc.copy(h_desc);

		}
		void time_step(float dt,int adv_type);
		void pcg_projection(float dt);
		//void output(uint nx, uint ny, uint nz, int frame, char* file_path);
		void compute_face_weights();
		void set_boundary_phi(float (*phi)(const Vec3f&));
		void set_boundary_vel(Vec3f (*vel)(const Vec3f&));
		float get_boundary_coef(Vec3f & pos);
	};


}



#endif