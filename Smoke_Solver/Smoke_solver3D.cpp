#include "Smoke_solver3D.h"
#include <fstream> 
#include <math.h>
#include "levelset_util.h"
using namespace std;
namespace ivock{

	void SmokeSolver3D::advect( float dt )
	{
		_utemp.setZero();
		_vtemp.setZero();
		_wtemp.setZero();
		


		advect_field_cubic_clamp(dt,_un,_utemp);
		advect_field_cubic_clamp(dt,_vn,_vtemp);
		advect_field_cubic_clamp(dt,_wn,_wtemp);

		


		_un.copy(_utemp);
		_vn.copy(_vtemp);
		_wn.copy(_wtemp);

		


	}
	
	void SmokeSolver3D::getCFL(float dt )
	{
		float max_v = _hx;
		for (int k=0; k<_nz;k++)for (int j=0; j<_ny;j++)for (int i=0; i<_nx+1;i++)
		{
			if (fabs(_un(i,j,k))>max_v)
			{
				max_v = _un(i,j,k);
			}
		}
		for (int k=0; k<_nz;k++)for (int j=0; j<_ny+1;j++)for (int i=0; i<_nx;i++)
		{
			if (fabs(_vn(i,j,k))>max_v)
			{
				max_v = _vn(i,j,k);
			}
		}
		for (int k=0; k<_nz+1;k++)for (int j=0; j<_ny;j++)for (int i=0; i<_nx;i++)
		{
			if (fabs(_wn(i,j,k))>max_v)
			{
				max_v = _wn(i,j,k);
			}
		}
		_cfl = _hx/max_v;
		printf("max_v: %f, dt: %f, cfl radio:%f\n",max_v,dt,dt/_cfl);

	}

	Vec3f SmokeSolver3D::get_solid_vel(Vec3f & pos)
	{
		float u = _u_solid.sample_linear(pos[0],pos[1],pos[2]);
		float v = _v_solid.sample_linear(pos[0],pos[1],pos[2]);
		float w = _w_solid.sample_linear(pos[0],pos[1],pos[2]);
		Vec3f vel(u,v,w);
		if(sub_domain.is_in(pos)) vel = sub_domain.get_solid_vel(pos);
		return vel;
	}

	inline Vec3f SmokeSolver3D::get_velocity(Vec3f & pos)
	{
		//if(sub_domain.is_in2(pos)) return  sub_domain.get_velocity(pos);
		float u = _un.sample_linear(pos[0],pos[1],pos[2]);
		float v = _vn.sample_linear(pos[0],pos[1],pos[2]);
		float w = _wn.sample_linear(pos[0],pos[1],pos[2]);
		Vec3f vel(u,v,w);

		//!!! commented out
		Vec3f vel2;
		if(sub_domain.is_in2(pos))
		{
			vel2 = sub_domain.get_velocity(pos);
			float c = sub_domain.get_boundary_coef(pos);
			vel = c*vel2 + (1.0f-c)*vel;
		}
		return vel;
	}
	inline Vec3f SmokeSolver3D::get_velocity_self(Vec3f & pos)
	{
		
		float u = _un.sample_linear(pos[0],pos[1],pos[2]);
		float v = _vn.sample_linear(pos[0],pos[1],pos[2]);
		float w = _wn.sample_linear(pos[0],pos[1],pos[2]);
		//Vec3f vel(u,v,w);

		return Vec3f(u,v,w);
	}
	Vec3f SmokeSolver3D::get_velocity_diff(Vec3f & pos)
	{
		float u = _utemp.sample_linear(pos[0],pos[1],pos[2]);
		float v = _vtemp.sample_linear(pos[0],pos[1],pos[2]);
		float w = _wtemp.sample_linear(pos[0],pos[1],pos[2]);
		return Vec3f(u,v,w);
	}
	float SmokeSolver3D::get_solid_phi(Vec3f & pos)
	{
		float val =  _nodal_solid_phi.sample_linear(pos[0],pos[1],pos[2]);
		
		//!!! commented out
		if(sub_domain.is_in(pos)) val = sub_domain.get_solid_phi(pos);
		
		return val;
	}
	Vec3f SmokeSolver3D::get_solid_normal(Vec3f & pos)
	{
		Vec3f val =  _nodal_solid_phi.sample_grad(pos[0],pos[1],pos[2]);
		
		//!!! commented out
		if(sub_domain.is_in(pos)) val = sub_domain.get_solid_normal(pos);
		
		normalize(val);
		return val;
	}

	Vec3f SmokeSolver3D::traceRK1(Vec3f & pos, float dt)
	{
		
		Vec3f input = pos;
		Vec3f velocity1 = get_velocity(input);
		input = pos + dt * velocity1;
		return input;
	}

	Vec3f SmokeSolver3D::traceRK2(Vec3f & pos, float dt)
	{
		
		Vec3f input = pos;
		Vec3f velocity1 = get_velocity(input);
		Vec3f midp1 = input + ((float)(0.5*dt))*velocity1;
		input = pos + dt * get_velocity(midp1);
		return input;
	}
	Vec3f SmokeSolver3D::traceRK3(Vec3f & pos, float dt)
	{
		float c1 = 2.0/9.0*dt, c2 = 3.0/9.0 * dt, c3 = (1.0-c1-c2) * dt;
		Vec3f input = pos;
		Vec3f velocity1 = get_velocity(input);
		Vec3f midp1 = input + ((float)(0.5*dt))*velocity1;
		Vec3f velocity2 = get_velocity(midp1);
		Vec3f midp2 = input + ((float)(0.75*dt))*velocity2;
		Vec3f velocity3 = get_velocity(midp2);
		//velocity = get_velocity(input + 0.5f*dt*velocity);
		//input += dt*velocity;
		input = input + c1*velocity1 + c2*velocity2 + c3*velocity3;
		return input;
	}
	Vec3f SmokeSolver3D::trace(float dt, Vec3f & pos)
	{
		if (dt>0)
		{
			float T=dt;
			Vec3f opos=pos;
			float t = 0;
			float substep = _cfl;
			while(t < T) {

				if(t + substep > T)
					substep = T - t;
				opos = traceRK3(opos,-substep);
				t+=substep;
			}
			//opos[0] = min(max(0.0f,opos[0]),(float)(_nx-1)*_hx );
			//opos[1] = min(max(0.0f,opos[1]),(float)(_ny-1)*_hx );
			//opos[2] = min(max(0.0f,opos[2]),(float)(_nz-1)*_hx );
			return opos;
		}
		else
		{
			float T = -dt;
			Vec3f opos=pos;
			float t = 0;
			float substep = _cfl;
			while(t < T) {

				if(t + substep > T)
					substep = T - t;
				opos = traceRK3(opos,substep);
				t+=substep;
			}
			//opos[0] = min(max(0.0f,opos[0]),(float)(_nx-1)*_hx );
			//opos[1] = min(max(0.0f,opos[1]),(float)(_ny-1)*_hx );
			//opos[2] = min(max(0.0f,opos[2]),(float)(_nz-1)*_hx );
			return opos;
		}

	}

	void SmokeSolver3D::grid_diffuse(float dt)
	{
		_unp1.copy(_un);
		_vnp1.copy(_vn);
		_wnp1.copy(_wn);
		int compute_num = _un._nx*_un._ny*_un._nz;
		int slice = _un._nx*_un._ny;
		tbb::parallel_for(0, compute_num, 1, [&](int thread_idx)
		{
			int k = thread_idx/slice;
			int j = (thread_idx%slice)/_un._nx;
			int i = thread_idx%_un._nx;
			if(k<_un._nz-1&&k>0 && j<_un._ny-1 &&j>0 && i<_un._nx-2 && i>1)
			{
				if(_u_weights(i,j,k)==0)
				{
					_unp1(i,j,k) = _u_solid(i,j,k);
				}
			}
		});
		compute_num = _vn._nx*_vn._ny*_vn._nz;
		slice = _vn._nx*_vn._ny;
		tbb::parallel_for(0, compute_num, 1, [&](int thread_idx)
		{
			int k = thread_idx/slice;
			int j = (thread_idx%slice)/_vn._nx;
			int i = thread_idx%_vn._nx;
			if(k<_vn._nz-1 &&k>0 && j>1 && j<_vn._ny-2 && i<_vn._nx-1 && i>0 )
			{
				if(_v_weights(i,j,k)==0)
				{
					_vnp1(i,j,k) = _v_solid(i,j,k);
				}
			}
		});
		compute_num = _wn._nx*_wn._ny*_wn._nz;
		slice = _wn._nx*_wn._ny;
		tbb::parallel_for(0, compute_num, 1, [&](int thread_idx)
		{
			int k = thread_idx/slice;
			int j = (thread_idx%slice)/_wn._nx;
			int i = thread_idx%_wn._nx;
			if(k>1 && k<_wn._nz-2 && j<_wn._ny-1&&j>0 && i<_wn._nx-1&&i>0 )
			{
				if(_u_weights(i,j,k)==0)
				{
					_wnp1(i,j,k) = _w_solid(i,j,k);
				}
			}
		});



		compute_num = _un._nx*_un._ny*_un._nz;
		slice = _un._nx*_un._ny;
		tbb::parallel_for(0, compute_num, 1, [&](int thread_idx)
		{
			int k = thread_idx/slice;
			int j = (thread_idx%slice)/_un._nx;
			int i = thread_idx%_un._nx;
			if(k<_un._nz-1&&k>0 && j<_un._ny-1 &&j>0 && i<_un._nx-2 && i>1)
			{
				float dudx2 = _unp1.at(i+1,j,k) - 2.0f*_unp1.at(i,j,k) + _unp1.at(i-1,j,k);
				float dudy2 = _unp1.at(i,j+1,k) - 2.0f*_unp1.at(i,j,k) + _unp1.at(i,j-1,k);
				float dudz2 = _unp1.at(i,j,k+1) - 2.0f*_unp1.at(i,j,k) + _unp1.at(i,j,k-1);
				float dvdxy = _vnp1.at(i,j+1,k) - _vnp1.at(i-1,j+1,k) - _vnp1.at(i,j,k) + _vnp1.at(i-1,j,k);
				float dwdzx = _wnp1.at(i,j,k+1) - _wnp1.at(i-1,j,k+1) - _wnp1.at(i,j,k) + _wnp1.at(i-1,j,k);

				float u_diffuse = dt*nu/sqr(_hx)*(2.0*dudx2 + dudy2 + dudz2 + dvdxy + dwdzx);
				_un(i,j,k) += u_diffuse;

			}
		});
		compute_num = _vn._nx*_vn._ny*_vn._nz;
		slice = _vn._nx*_vn._ny;
		tbb::parallel_for(0, compute_num, 1, [&](int thread_idx)
		{
			int k = thread_idx/slice;
			int j = (thread_idx%slice)/_vn._nx;
			int i = thread_idx%_vn._nx;
			if(k<_vn._nz-1 &&k>0 && j>1 && j<_vn._ny-2 && i<_vn._nx-1 && i>0 )
			{
				
				float dvdx2 = _vnp1.at(i+1,j,k) - 2.0f*_vnp1.at(i,j,k) + _vnp1.at(i-1,j,k);
				float dvdy2 = _vnp1.at(i,j+1,k) - 2.0f*_vnp1.at(i,j,k) + _vnp1.at(i,j-1,k);
				float dvdz2 = _vnp1.at(i,j,k+1) - 2.0f*_vnp1.at(i,j,k) + _vnp1.at(i,j,k-1);
				float dudyx = _unp1.at(i+1,j,k) - _unp1.at(i,j,k) - _unp1.at(i+1,j-1,k) + _unp1.at(i,j-1,k);
				float dwdyz = _wnp1.at(i,j,k+1) - _wnp1.at(i,j,k+1) - _wnp1.at(i,j-1,k+1) + _wnp1.at(i,j-1,k);

				float v_diffuse = dt*nu/sqr(_hx)*(dvdx2 + 2.0*dvdy2 + dvdz2 + dudyx + dwdyz);
				_vn(i,j,k) += v_diffuse;
			}
		});
		compute_num = _wn._nx*_wn._ny*_wn._nz;
		slice = _wn._nx*_wn._ny;
		tbb::parallel_for(0, compute_num, 1, [&](int thread_idx)
		{
			int k = thread_idx/slice;
			int j = (thread_idx%slice)/_wn._nx;
			int i = thread_idx%_wn._nx;
			if(k>1 && k<_wn._nz-2 && j<_wn._ny-1&&j>0 && i<_wn._nx-1&&i>0 )
			{
				float dwdx2 = _wnp1.at(i+1,j,k) - 2.0f*_wnp1.at(i,j,k) + _wnp1.at(i-1,j,k);
				float dwdy2 = _wnp1.at(i,j+1,k) - 2.0f*_wnp1.at(i,j,k) + _wnp1.at(i,j-1,k);
				float dwdz2 = _wnp1.at(i,j,k+1) - 2.0f*_wnp1.at(i,j,k) + _wnp1.at(i,j,k-1);
				float dudzx = _unp1.at(i+1,j,k) - _unp1.at(i,j,k) - _unp1.at(i+1,j,k-1) + _unp1.at(i,j,k-1);
				float dvdzy = _vnp1.at(i,j+1,k) - _vnp1.at(i,j,k) - _vnp1.at(i,j+1,k-1) + _vnp1.at(i,j,k-1);

				float w_diffuse = dt*nu/sqr(_hx)*(dwdx2 + dwdy2 + 2.0*dwdz2 + dudzx + dvdzy);
				_wn(i,j,k) += w_diffuse;
			}
		});
	}

	void SmokeSolver3D::time_step( float dt , int adv_type, int frame)
	{
		//printf("entering time step\n");
		getCFL(dt);
		float t = 0;
		
		


		
		emit_tracers();
		while(t < dt) {
			float substep = 100*_cfl;   
			if(t + substep > dt)
				substep = dt - t;
			
			
			advect_tracers(substep);
			_utemp.setZero();
			_vtemp.setZero();
			_wtemp.setZero();


			
			//FLIP.reserve(FLIP.size()+sub_domain.FLIP.size());
			

			diffuse(substep);
			for (int i=0;i<sub_domain.FLIP.size();i++)
			{
				FLIP.push_back(sub_domain.FLIP[i]);
			}
			
			printf("after seed halo %d\n",FLIP.size());
			//!!! replaced to
			//_unp1.copy(_un);
			//_vnp1.copy(_vn);
			//_wnp1.copy(_wn);

			advect_FLIP(substep);
			hash_FLIP();

			/* not really useful if only dealing with flow past obstacles;

			printf("entering advection\n");
			_temp.setZero();
			advect_field(substep, _rho, _temp);
			_rho.copy(_temp);
			_temp.setZero();
			
			advect_field(substep, _Temperature, _temp);
			_Temperature.copy(_temp);
			_temp.setZero();

			printf("smoke advection done\n");
			gen_heat(substep);
			add_smokeforce(substep);
			*/
			
			particle_to_grid(_un,_unp1,0);
			particle_to_grid(_vn,_vnp1,1);
			particle_to_grid(_wn,_wnp1,2);
			_utemp.copy(_un);
			_vtemp.copy(_vn);
			_wtemp.copy(_wn);

			

			pcg_projection(1.0);//you can put dt here, but actually any constant can be used... :)
			
			//extrapolate(u_extrap,_un,u_valid);
			//extrapolate(v_extrap,_vn,v_valid);
			//extrapolate(w_extrap,_wn,w_valid);
			//constrain_velocity();
			compute_delta_u();
			particle_from_grid();
			printf("particle from grid\n");
			
			
			//!!! commented out
			sub_domain.FLIP.resize(0);
			sub_domain.FLIP.reserve(FLIP.size());
			for (int i=0;i<FLIP.size();i++)
			{
				if (sub_domain.is_halo(FLIP[i].pos, 0.5*_hx))
				{
					sub_domain.FLIP.push_back(FLIP[i]);
				}
			}
			printf("before hash\n");
			sub_domain.hash_FLIP();
			printf("hash\n");
			sub_domain.particle_to_grid();
			printf("particle to grid\n");
			sub_domain.pcg_projection(1.0);
			sub_domain.compute_delta_u();
			sub_domain.particle_from_grid();
			sub_domain.reorder_FLIP();
			reorder_FLIP();
			printf("%d\n",FLIP.size());


			t+=substep;
		}
	}
	void SmokeSolver3D::constrain_velocity()
	{
		_unp1.setZero();
		_vnp1.setZero();
		_wnp1.setZero();
		_unp1.copy(_un);
		_vnp1.copy(_vn);
		_wnp1.copy(_wn);
		int compute_num; int slice;
		compute_num = _un._nx*_un._ny*_un._nz;
		slice = _un._nx*_un._ny;
		tbb::parallel_for(0, compute_num, 1, [&](int thread_idx)
		{
			int k = thread_idx/slice;
			int j = (thread_idx%slice)/_un._nx;
			int i = thread_idx%_un._nx;
			if(k<_un._nz-1&&k>0 && j<_un._ny-1 &&j>0 && i<_un._nx-2 && i>1)
			{
				if (_u_weights(i,j,k)==0)
				{
					Vec3f pos = Vec3f(i,j,k) * _hx 
						- Vec3f(_un._ox,_un._oy,_un._oz) * _hx;
				
					Vec3f vel = get_velocity_self(pos);
					Vec3f normal = get_solid_normal(pos); 
					float perp_component = dot(vel, normal);
					vel -= perp_component*normal;
					Vec3f v_solid = get_solid_vel(pos);
					vel += normal * dot(v_solid,normal);
					_unp1(i,j,k) = vel[0];
				
				}
				

			}
		});

		
		compute_num = _vn._nx*_vn._ny*_vn._nz;
		slice = _vn._nx*_vn._ny;
		tbb::parallel_for(0, compute_num, 1, [&](int thread_idx)
		{
			int k = thread_idx/slice;
			int j = (thread_idx%slice)/_vn._nx;
			int i = thread_idx%_vn._nx;
			if(k<_vn._nz-1 &&k>0 && j>1 && j<_vn._ny-2 && i<_vn._nx-1 && i>0 )
			{
				if (_v_weights(i,j,k)==0)
				{
					Vec3f pos = Vec3f(i,j,k) * _hx 
						- Vec3f(_vn._ox,_vn._oy,_vn._oz) * _hx;

					Vec3f vel = get_velocity_self(pos);
					Vec3f normal = get_solid_normal(pos); 
					float perp_component = dot(vel, normal);
					vel -= perp_component*normal;
					Vec3f v_solid = get_solid_vel(pos);
					vel += normal * dot(v_solid,normal);
					_vnp1(i,j,k) = vel[1];

				}

			}
		});

		
		compute_num = _wn._nx*_wn._ny*_wn._nz;
		slice = _wn._nx*_wn._ny;
		tbb::parallel_for(0, compute_num, 1, [&](int thread_idx)
		{
			int k = thread_idx/slice;
			int j = (thread_idx%slice)/_wn._nx;
			int i = thread_idx%_wn._nx;
			if(k>1 && k<_wn._nz-2 && j<_wn._ny-1&&j>0 && i<_wn._nx-1&&i>0 )
			{
				if (_w_weights(i,j,k)==0)
				{
					Vec3f pos = Vec3f(i,j,k) * _hx 
						- Vec3f(_wn._ox,_wn._oy,_wn._oz) * _hx;

					Vec3f vel = get_velocity_self(pos);
					Vec3f normal = get_solid_normal(pos); 
					float perp_component = dot(vel, normal);
					vel -= perp_component*normal;
					Vec3f v_solid = get_solid_vel(pos);
					vel += normal * dot(v_solid,normal);
					_wnp1(i,j,k) = vel[2];

				}
			}
		});



		compute_num = _un._nx*_un._ny*_un._nz;
		slice = _un._nx*_un._ny;
		tbb::parallel_for(0, compute_num, 1, [&](int thread_idx)
		{
			int k = thread_idx/slice;
			int j = (thread_idx%slice)/_un._nx;
			int i = thread_idx%_un._nx;
			if(k<_un._nz-1&&k>0 && j<_un._ny-1 &&j>0 && i<_un._nx-2 && i>1)
			{
				if (_u_weights(i,j,k)==0)
				{ 
					_un(i,j,k) = _unp1(i,j,k);

				}


			}
		});
		compute_num = _vn._nx*_vn._ny*_vn._nz;
		slice = _vn._nx*_vn._ny;
		tbb::parallel_for(0, compute_num, 1, [&](int thread_idx)
		{
			int k = thread_idx/slice;
			int j = (thread_idx%slice)/_vn._nx;
			int i = thread_idx%_vn._nx;
			if(k<_vn._nz-1 &&k>0 && j>1 && j<_vn._ny-2 && i<_vn._nx-1 && i>0 )
			{
				if (_v_weights(i,j,k)==0)
				{
					_vn(i,j,k) = _vnp1(i,j,k) ;

				}

			}
		});
		compute_num = _wn._nx*_wn._ny*_wn._nz;
		slice = _wn._nx*_wn._ny;
		tbb::parallel_for(0, compute_num, 1, [&](int thread_idx)
		{
			int k = thread_idx/slice;
			int j = (thread_idx%slice)/_wn._nx;
			int i = thread_idx%_wn._nx;
			if(k>1 && k<_wn._nz-2 && j<_wn._ny-1&&j>0 && i<_wn._nx-1&&i>0 )
			{
				if (_w_weights(i,j,k)==0)
				{
					_wn(i,j,k) = _wnp1(i,j,k);

				}
			}
		});
	}
	void SmokeSolver3D::advect_field(float dt, buffer3Df & field, buffer3Df &field_new)
	{
		int compute_elements = field._blockx*field._blocky*field._blockz;

		int slice = field._blockx*field._blocky;

		tbb::parallel_for(0, compute_elements, 1, [&](int thread_idx) {

			uint bk = thread_idx/slice;
			uint bj = (thread_idx%slice)/field._blockx;
			uint bi = thread_idx%(field._blockx);

			for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
			{
				uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
				if(i<field._nx && j<field._ny && k<field._nz)
				{

					float world_x = ((float)i-field._ox)*field._hx;
					float world_y = ((float)j-field._oy)*field._hy;
					float world_z = ((float)k-field._oz)*field._hz;
					Vec3f pos(world_x,world_y,world_z);
					Vec3f back_pos = traceRK3(pos,-dt);

					////Vec3f trace_pos = trace(dt, pos);

					//float u = _un.sample_linear(world_x,world_y,world_z);
					//float v = _vn.sample_linear(world_x,world_y,world_z);
					//float w = _wn.sample_linear(world_x,world_y,world_z);

					//float px = world_x - 0.5*dt * u, py = world_y - 0.5*dt *v, pz = world_z - 0.5*dt*w;
					//u = _un.sample_linear(px,py,pz);
					//v = _vn.sample_linear(px,py,pz);
					//w = _wn.sample_linear(px,py,pz);

					//px = world_x - dt*u; py = world_y - dt*v; pz = world_z - dt*w;

					float SLv = field.sample_linear(back_pos[0],
						back_pos[1],
						back_pos[2]);

					field_new(i,j,k) = SLv;

				}
			}
		});
	}
	void SmokeSolver3D::advect_field_cubic( float dt, buffer3Df & field, buffer3Df &field_new )
	{
		int compute_elements = field._blockx*field._blocky*field._blockz;

		int slice = field._blockx*field._blocky;

		tbb::parallel_for(0, compute_elements, 1, [&](int thread_idx) {

			uint bk = thread_idx/slice;
			uint bj = (thread_idx%slice)/field._blockx;
			uint bi = thread_idx%(field._blockx);

			for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
			{
				uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
				if(i<field._nx && j<field._ny && k<field._nz)
				{

					float world_x = ((float)i-field._ox)*_hx;
					float world_y = ((float)j-field._oy)*_hy;
					float world_z = ((float)k-field._oz)*_hz;
					Vec3f pos(world_x,world_y,world_z);
					Vec3f trace_pos = trace(dt, pos);

					//float u = _un.sample_linear(world_x,world_y,world_z);
					//float v = _vn.sample_linear(world_x,world_y,world_z);
					//float w = _wn.sample_linear(world_x,world_y,world_z);

					//float px = world_x - 0.5*dt * u, py = world_y - 0.5*dt *v, pz = world_z - 0.5*dt*w;
					//u = _un.sample_linear(px,py,pz);
					//v = _vn.sample_linear(px,py,pz);
					//w = _wn.sample_linear(px,py,pz);

					//px = world_x - dt*u; py = world_y - dt*v; pz = world_z - dt*w;

					float SLv = field.sample_cubic(trace_pos[0],trace_pos[1],trace_pos[2]);

					field_new(i,j,k) = SLv;



				}
			}
		});
	}

	void SmokeSolver3D::advect_field_cubic_clamp(float dt, buffer3Df & field, buffer3Df &field_new)
	{
		int compute_elements = field._blockx*field._blocky*field._blockz;

		int slice = field._blockx*field._blocky;

		tbb::parallel_for(0, compute_elements, 1, [&](int thread_idx) {

			uint bk = thread_idx/slice;
			uint bj = (thread_idx%slice)/field._blockx;
			uint bi = thread_idx%(field._blockx);

			for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
			{
				uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
				if(i<field._nx && j<field._ny && k<field._nz)
				{

					float world_x = ((float)i-field._ox)*_hx;
					float world_y = ((float)j-field._oy)*_hy;
					float world_z = ((float)k-field._oz)*_hz;
					Vec3f pos(world_x,world_y,world_z);
					//Vec3f trace_pos = trace(dt, pos);
					//
					float u = _un.sample_linear(world_x,world_y,world_z);
					float v = _vn.sample_linear(world_x,world_y,world_z);
					float w = _wn.sample_linear(world_x,world_y,world_z);

					float px = world_x - 0.5*dt * u, py = world_y - 0.5*dt *v, pz = world_z - 0.5*dt*w;
					u = _un.sample_linear(px,py,pz);
					v = _vn.sample_linear(px,py,pz);
					w = _wn.sample_linear(px,py,pz);

					px = world_x-dt*u, py = world_y-dt*v, pz = world_z-dt*w;
					float v0,v1,v2,v3,v4,v5,v6,v7;
					float SLvl = field.sample_cube_lerp(px,py,pz,
						v0,v1,v2,v3,v4,v5,v6,v7);

					float SLvc = field.sample_cubic(px,py,pz);

					float minv = min(min(min(min(min(min(min(v0,v1),v2),v3),v4),v5),v6),v7);
					float maxv = max(max(max(max(max(max(max(v0,v1),v2),v3),v4),v5),v6),v7);

					field_new(i,j,k) = SLvc;
					if( SLvc<=minv || SLvc >= maxv )
						field_new(i,j,k) = SLvl;

				}
			}
		});
	}

	//output to a pbrt scene....................
	//void SmokeSolver3D::output( uint nx, uint ny, uint nz, int frame, char* file_path)
	//{
	//	char file_name[256];
	//	sprintf(file_name,"%s/density_render.%04d.pbrt", file_path,frame);


	//	float *field = new float[nx*ny*nz];
	//	float *field_t=new float[nx*ny*nz];

	//	// normalize values
	//	//float max_dens = _rho(0,0,0);
	//	//for(int k=0;k<nz;k++)for(int j=0;j<ny;j++)for(int i=0;i<nx;i++) 
	//	//{
	//	//	if (_rho(0,0,0)>max_dens)
	//	//	{
	//	//		max_dens = _rho(i,j,k);
	//	//	}
	//	//}
	//	for(int k=0;k<nz;k++)for(int j=0;j<ny;j++)for(int i=0;i<nx;i++) 
	//	{
	//		if(_b_desc(i,j,k)==0){

	//			field[k*nx*ny + j*nx + i] = 0.5*max(_rho(i,j,k),0.0f)/10.0;
	//			field_t[k*nx*ny + j*nx + i] = _Tbf(i,j,k);
	//		}
	//		else
	//		{
	//			field[k*nx*ny + j*nx + i] =0;
	//			field_t[k*nx*ny + j*nx + i] = 0;

	//		}
	//	}


	//	std::fstream fout;
	//	fout.open(file_name, std::ios::out);

	//	int maxRes = (nx > ny) ? nx : ny;
	//	maxRes = (maxRes > nz) ? maxRes : nz;

	//	const float xSize = 1.0 / (float)maxRes * (float)nx;
	//	const float ySize = 1.0 / (float)maxRes * (float)ny;
	//	const float zSize = 1.0 / (float)maxRes * (float)nz;





	//	// dimensions
	//	fout<<"Volume \"volumegrid\" \n";
	//	fout<<" \"integer nx\""<<nx<<"\n";
	//	fout<<" \"integer ny\""<<ny<<"\n";
	//	fout<<" \"integer nz\""<<nz<<"\n";
	//	fout<<" \"point p0\" [ "<<-xSize/2.0<<" "<<0.0<<" "<<-zSize/2.0<<"] \"point p1\" "<<"[ "<<xSize/2.0<<" "<<ySize<<" "<<zSize/2.0<<" ]"<<"\n";
	//	fout<<" \"float density\" [ \n";
	//	for (int i = 0; i < nx*ny*nz; i++) 
	//		fout<<field[i]<<" ";
	//	fout<<"] \n";
	//	fout.close();

	//	//sprintf(file_name,"%s/density_render.%04d.bin", file_path,frame);
	//	//FILE *data_file = fopen(file_name,"wb");
	//	//fwrite(field,sizeof(float),nx*ny*nz,data_file);
	//	//fclose(data_file);

	//	//char file_name2[256];
	//	//sprintf(file_name2,"%s/density.%04d.bin", file_path,frame);
	//	//FILE *data_file = fopen(file_name2,"wb");
	//	//fwrite(field,sizeof(float),nx*ny*nz,data_file);
	//	//fclose(data_file);

	//	//char file_name3[256];
	//	//sprintf(file_name3,"%s/temperature.%04d.bin", file_path,frame);
	//	//data_file = fopen(file_name3,"wb");
	//	//fwrite(field_t,sizeof(float),nx*ny*nz,data_file);
	//	//fclose(data_file);




	//	delete[] field;
	//	delete[] field_t;
	//}
	void SmokeSolver3D::set_boundary_vel(Vec3f (*vel)(const Vec3f&))
	{
		_u_solid.setZero();
		_v_solid.setZero();
		_w_solid.setZero();
		int compute_num, slice;
		compute_num = _un._nx*_un._ny*_un._nz;
		slice = _un._nx*_un._ny;
		tbb::parallel_for(0, compute_num, 1, [&](int thread_idx)
		{
			int k = thread_idx/slice;
			int j = (thread_idx%slice)/_un._nx;
			int i = thread_idx%_un._nx;
			if(k<_un._nz&&k>=0 && j<_un._ny &&j>=0 && i<_un._nx && i>=0)
			{
				Vec3f pos = Vec3f(i,j,k)*_hx - Vec3f(_un._ox,_un._oy,_un._oz)*_hx;
				Vec3f velocity = vel(pos);
				_u_solid(i,j,k) = velocity[0];

			}
		});

		compute_num = _vn._nx*_vn._ny*_vn._nz;
		slice = _vn._nx*_vn._ny;
		tbb::parallel_for(0, compute_num, 1, [&](int thread_idx)
		{
			int k = thread_idx/slice;
			int j = (thread_idx%slice)/_vn._nx;
			int i = thread_idx%_vn._nx;
			if(k<_vn._nz &&k>=0 && j>=0 && j<_vn._ny && i<_vn._nx && i>=0 )
			{
				Vec3f pos = Vec3f(i,j,k)*_hx - Vec3f(_vn._ox,_vn._oy,_vn._oz)*_hx;
				Vec3f velocity = vel(pos);
				_v_solid(i,j,k) = velocity[1];

			}
		});

		compute_num = _wn._nx*_wn._ny*_wn._nz;
		slice = _wn._nx*_wn._ny;
		tbb::parallel_for(0, compute_num, 1, [&](int thread_idx)
		{
			int k = thread_idx/slice;
			int j = (thread_idx%slice)/_wn._nx;
			int i = thread_idx%_wn._nx;
			if(k>=0 && k<_wn._nz && j<_wn._ny&&j>=0 && i<_wn._nx&&i>=0 )
			{
				Vec3f pos = Vec3f(i,j,k)*_hx - Vec3f(_wn._ox,_wn._oy,_wn._oz)*_hx;
				Vec3f velocity = vel(pos);
				_w_solid(i,j,k) = velocity[2];
			}
		});

		//!!! commented out
		sub_domain.set_boundary_vel(vel);
	}


	void SmokeSolver3D::set_boundary_phi(float (*phi)(const Vec3f&))
	{
		int compute_num; int slice;
		compute_num = _nodal_solid_phi._nx*_nodal_solid_phi._ny*_nodal_solid_phi._nz;
		slice = _nodal_solid_phi._nx*_nodal_solid_phi._ny;
		tbb::parallel_for(0,compute_num,1,[&](int thread_idx)
		{
			int k = thread_idx/slice;
			int j = (thread_idx%slice)/_nodal_solid_phi._nx;
			int i = thread_idx%_nodal_solid_phi._nx;
			if(k<_nodal_solid_phi._nz&&j<_nodal_solid_phi._ny&&i<_nodal_solid_phi._nx)
			{

				Vec3f pos = Vec3f(i,j,k) * _hx - Vec3f(0.5,0.5,0.5)*_hx;
				_nodal_solid_phi(i,j,k) = phi(pos);
			}
		});

		//!!! commented out
		sub_domain.set_boundary_phi(phi);
		
		printf("set_boundary_phi done\n");
	}
	void SmokeSolver3D::compute_face_weights()
	{
		_u_weights.setZero();
		_v_weights.setZero();
		_w_weights.setZero();
		int compute_num; int slice;
		compute_num = _un._nx*_un._ny*_un._nz;
		slice = _un._nx*_un._ny;
		tbb::parallel_for(0,compute_num,1,[&](int thread_idx)
		{
			int k = thread_idx/slice;
			int j = (thread_idx%slice)/_un._nx;
			int i = thread_idx%_un._nx;
			if(i<_u_weights._nx&&j<_u_weights._ny&&k<_u_weights._nz)
			{

				_u_weights(i,j,k) = 1 - fraction_inside(_nodal_solid_phi(i,j,  k),

					_nodal_solid_phi(i,j+1,k),

					_nodal_solid_phi(i,j,  k+1),

					_nodal_solid_phi(i,j+1,k+1));

				_u_weights(i,j,k) = clamp(_u_weights(i,j,k),0.0f,1.0f);
			}
		});
		

		compute_num = _vn._nx*_vn._ny*_vn._nz;
		slice = _vn._nx*_vn._ny;
		tbb::parallel_for(0,compute_num,1,[&](int thread_idx)
		{
			int k = thread_idx/slice;
			int j = (thread_idx%slice)/_vn._nx;
			int i = thread_idx%_vn._nx;
			
			if(i<_v_weights._nx&&j<_v_weights._ny&&k<_v_weights._nz)
			{

				_v_weights(i,j,k) = 1 - fraction_inside(_nodal_solid_phi(i,  j,k),

					_nodal_solid_phi(i,  j,k+1),

					_nodal_solid_phi(i+1,j,k),

					_nodal_solid_phi(i+1,j,k+1));

				_v_weights(i,j,k) = clamp(_v_weights(i,j,k),0.0f,1.0f);
			}
		});

		compute_num = _wn._nx*_wn._ny*_wn._nz;
		slice = _wn._nx*_wn._ny;
		tbb::parallel_for(0,compute_num,1,[&](int thread_idx)
		{
			int k = thread_idx/slice;
			int j = (thread_idx%slice)/_wn._nx;
			int i = thread_idx%_wn._nx;
			if(i<_w_weights._nx&&j<_w_weights._ny&&k<_w_weights._nz)
			{
				_w_weights(i,j,k) = 1 - fraction_inside(_nodal_solid_phi(i,  j,  k),
					_nodal_solid_phi(i,  j+1,k),
					_nodal_solid_phi(i+1,j,  k),
					_nodal_solid_phi(i+1,j+1,k));
				_w_weights(i,j,k) = clamp(_w_weights(i,j,k),0.0f,1.0f);
			}
		});


	}
	void SmokeSolver3D::pcg_projection(float dt)
	{
		compute_face_weights();

		int ni = _nx;
		int nj = _ny;
		int nk = _nz;

		int system_size = ni*nj*nk;
		if(rhs.size() != system_size) {
			rhs.resize(system_size);
			pressure.resize(system_size);
			matrix.resize(system_size);
		}

		matrix.zero();
		rhs.assign(rhs.size(), 0);
		pressure.assign(pressure.size(), 0);
		//write boundary velocity;
		int compute_num = ni*nj*nk;
		int slice = ni*nj;
		tbb::parallel_for(0,compute_num,1,[&](int thread_idx)
		{
			int k = thread_idx/slice;
			int j = (thread_idx%slice)/ni;
			int i = thread_idx%ni;
			if ( _b_desc(i,j,k)==1 )//in flow
			{
				_un(i,j,k) = _U_in;//5.0f;
				_un(i+1,j,k) = _U_in;//5.0f;
				_vn(i,j,k) = 0;
				_vn(i,j+1,k) = 0;
				_wn(i,j,k) = 0;
				_wn(i,j,k+1) = 0;
			}
		});


		//set up solver
		compute_num = ni*nj*nk;
		slice = ni*nj;
		tbb::parallel_for(0,compute_num,1,[&](int thread_idx)
		{
			int k = thread_idx/slice;
			int j = (thread_idx%slice)/ni;
			int i = thread_idx%ni;
			if(i>=1 && i<ni-1 && j>=1 && j<nj-1 && k>=1 && k<nk-1)
			{
				int index = i + ni*j + ni*nj*k;

				rhs[index] = 0;
				pressure[index] = 0;

				//if( _b_desc(i,j,k)==0 )//a fluid cell 
				{

					//right neighbour
					//if( _b_desc(i+1,j,k)==0 ) {//a fluid cell
					//	matrix.add_to_element(index, index, 1.0/_hx/_hx);
					//	matrix.add_to_element(index, index + 1, -1.0/_hx/_hx);
					//}
					//else if( _b_desc(i+1,j,k)==1 )//an empty cell
					//{
					//	matrix.add_to_element(index, index, 1.0/_hx/_hx);
					//}
					//rhs[index] -= _un(i+1,j,k) / _hx;

					float term = _u_weights(i+1,j,k)  / sqr(_hx);
					matrix.add_to_element(index, index, term);
					if(i<ni-2) matrix.add_to_element(index,index+1, -term);
					rhs[index] -= _u_weights(i+1,j,k)*_un(i+1,j,k)/_hx;


					//left neighbour
					//if( _b_desc(i-1,j,k)==0 ) {
					//	matrix.add_to_element(index, index, 1.0/_hx/_hx);
					//	matrix.add_to_element(index, index - 1, -1.0/_hx/_hx);
					//}
					//else if( _b_desc(i-1,j,k)==1 ){

					//	matrix.add_to_element(index, index, 1.0/_hx/_hx);
					//}
					//rhs[index] += _un(i,j,k) / _hx;
					term = _u_weights(i,j,k)  /sqr(_hx);
					matrix.add_to_element(index,index,term);
					if(i>=2) matrix.add_to_element(index,index-1,-term);
					rhs[index] += _u_weights(i,j,k)*_un(i,j,k)/_hx;


					//top neighbour
					//if( _b_desc(i,j+1,k)==0 ) {//a fluid cell
					//	matrix.add_to_element(index, index, 1.0/_hx/_hx);
					//	matrix.add_to_element(index, index + ni, -1.0/_hx/_hx);
					//}
					//else if( _b_desc(i,j+1,k)==1 )//an empty cell
					//{
					//	matrix.add_to_element(index, index, 1.0/_hx/_hx);
					//}
					//rhs[index] -= _vn(i,j+1,k) / _hx;
					term = _v_weights(i,j+1,k)  /sqr(_hx);
					matrix.add_to_element(index, index,term);
					if(j<nj-2) matrix.add_to_element(index,index+ni,-term);
					rhs[index] -= _v_weights(i,j+1,k)*_vn(i,j+1,k)/_hx;

					//bottom neighbour
					//if( _b_desc(i,j-1,k)==0 ) {
					//	matrix.add_to_element(index, index, 1.0/_hx/_hx);
					//	matrix.add_to_element(index, index - ni, -1.0/_hx/_hx);
					//}
					//else if( _b_desc(i,j-1,k)==1 ){

					//	matrix.add_to_element(index, index, 1.0/_hx/_hx);
					//}
					//rhs[index] += _vn(i,j,k) / _hx;
					//rhs[index] += _burn_div(i,j,k);

					term = _v_weights(i,j,k)  /sqr(_hx);
					matrix.add_to_element(index,index,term);
					if(j>=2) matrix.add_to_element(index,index-ni,-term);
					rhs[index] += _v_weights(i,j,k)*_vn(i,j,k)/_hx;



					//back neighbour
					//if( _b_desc(i,j,k+1)==0 ) {//a fluid cell
					//	matrix.add_to_element(index, index, 1.0/_hx/_hx);
					//	matrix.add_to_element(index, index + ni*nj, -1.0/_hx/_hx);
					//}
					//else if( _b_desc(i,j,k+1)==1 )//an empty cell
					//{
					//	matrix.add_to_element(index, index, 1.0/_hx/_hx);
					//}
					//rhs[index] -= _wn(i,j,k+1) / _hx;
					term = _w_weights(i,j,k+1)  /sqr(_hx);
					matrix.add_to_element(index,index,term);
					if(k<nk-2) matrix.add_to_element(index,index+ni*nj,-term);
					rhs[index] -= _w_weights(i,j,k+1)*_wn(i,j,k+1)/_hx;

					//front neighbour
					/*if( _b_desc(i,j,k-1)==0 ) {
						matrix.add_to_element(index, index, 1.0/_hx/_hx);
						matrix.add_to_element(index, index - ni*nj, -1.0/_hx/_hx);
					}
					else if( _b_desc(i,j,k-1)==1 ){

						matrix.add_to_element(index, index, 1.0/_hx/_hx);
					}
					rhs[index] += _wn(i,j,k) / _hx;*/
					term = _w_weights(i,j,k)  /sqr(_hx);
					matrix.add_to_element(index,index,term);
					if(k>=2) matrix.add_to_element(index,index-ni*nj,-term);
					rhs[index] += _w_weights(i,j,k)*_wn(i,j,k)/_hx;

					//rhs[index] += _burn_div(i,j,k);
					if(matrix(index,index)>0)
					{

						rhs[index] += 1.0/_hx*((_u_weights(i+1,j,k)-1.0f)*_u_solid(i+1,j,k)
							-(_u_weights(i,j,k)  -1.0f)*_u_solid(i,j,k)
							+(_v_weights(i,j+1,k)-1.0f)*_v_solid(i,j+1,k)
							-(_v_weights(i,j,k)  -1.0f)*_v_solid(i,j,k)
							+(_w_weights(i,j,k+1)-1.0f)*_w_solid(i,j,k+1)
							-(_w_weights(i,j,k)  -1.0f)*_w_solid(i,j,k));
					}
				}
			}
		});

		//Solve the system using a AMGPCG solver

		double tolerance;
		int iterations;
		//solver.set_solver_parameters(1e-6, 1000);
		//bool success = solver.solve(matrix, rhs, pressure, tolerance, iterations);
		
		bool success = AMGPCGSolve(matrix,rhs,pressure,1e-6,50,tolerance,iterations,_nx,_ny,_nz);
		printf("Solver took %d iterations and had residual %e\n", iterations, tolerance);
		if(!success) {
			printf("WARNING: Pressure solve failed!************************************************\n");
		}

		//apply grad
		u_valid.assign(0);
		compute_num = _un._nx*_un._ny*_un._nz;
		slice = _un._nx*_un._ny;
		tbb::parallel_for(0, compute_num, 1, [&](int thread_idx)
		{
			int k = thread_idx/slice;
			int j = (thread_idx%slice)/_un._nx;
			int i = thread_idx%_un._nx;
			if(k<_un._nz-1&&k>0 && j<_un._ny-1 &&j>0 && i<_un._nx-2 && i>1)
			{
				int index = i + j*ni + k*ni*nj;
				//if(_b_desc(i,j,k) == 0 || _b_desc(i-1,j,k) == 0) {
				if(_u_weights(i,j,k) >0){
					_un(i,j,k) -=  (float)(pressure[index] - pressure[index-1]) / _hx ; 
					u_valid(i,j,k) = 1;
				}
				else
				{
					_un(i,j,k) = _u_solid(i,j,k);
				}

			}
		});

		v_valid.assign(0);
		compute_num = _vn._nx*_vn._ny*_vn._nz;
		slice = _vn._nx*_vn._ny;
		tbb::parallel_for(0, compute_num, 1, [&](int thread_idx)
		{
			int k = thread_idx/slice;
			int j = (thread_idx%slice)/_vn._nx;
			int i = thread_idx%_vn._nx;
			if(k<_vn._nz-1 &&k>0 && j>1 && j<_vn._ny-2 && i<_vn._nx-1 && i>0 )
			{
				int index = i + j*ni + k*ni*nj;
				//if(_b_desc(i,j,k) == 0 || _b_desc(i,j-1,k) == 0)
				if(_v_weights(i,j,k)>0)
				{

					_vn(i,j,k) -=  (float)(pressure[index] - pressure[index-ni]) / _hx ; 
					v_valid(i,j,k) = 1;
				}
				else
				{
					_vn(i,j,k) =  _v_solid(i,j,k);
				}

			}
		});

		w_valid.assign(0);
		compute_num = _wn._nx*_wn._ny*_wn._nz;
		slice = _wn._nx*_wn._ny;
		tbb::parallel_for(0, compute_num, 1, [&](int thread_idx)
		{
			int k = thread_idx/slice;
			int j = (thread_idx%slice)/_wn._nx;
			int i = thread_idx%_wn._nx;
			if(k>1 && k<_wn._nz-2 && j<_wn._ny-1&&j>0 && i<_wn._nx-1&&i>0 )
			{
				int index = i + j*ni + k*ni*nj;
				//if(_b_desc(i,j,k) == 0 || _b_desc(i,j,k-1) == 0) 
				if(_w_weights(i,j,k)>0)
				{

					_wn(i,j,k) -=  (float)(pressure[index] - pressure[index-ni*nj]) / _hx ; 
					w_valid(i,j,k) = 1;
				}
				else
				{
					_wn(i,j,k) = _w_solid(i,j,k);
				}
			}
		});
		
	}

	void SmokeSolver3D::extrapolate(Array3f & grid, buffer3Df & u, Array3c & valid)
	{
		grid.assign(0);
		int compute_num = u._nx*u._ny*u._nz;
		int slice = u._nx*u._ny;
		tbb::parallel_for(0, compute_num, 1, [&](int thread_idx)
		{
			int k = thread_idx/slice;
			int j = (thread_idx%slice)/u._nx;
			int i = thread_idx%u._nx;
			if(k<u._nz && j<u._ny && i<u._nx )
			{
				grid(i,j,k) = u(i,j,k);
			}
		});



		Array3f temp_grid = grid;
		Array3c old_valid(valid.ni,valid.nj,valid.nk);
		for(int layers = 0; layers < 10; ++layers) {
			old_valid = valid;
			int num = grid.ni*grid.nj*grid.nk;
			int s = grid.ni*grid.nj;
			tbb::parallel_for(0,num,1,[&](int thread_idx)
				//for(int k = 1; k < grid.nk-1; ++k) for(int j = 1; j < grid.nj-1; ++j) for(int i = 1; i < grid.ni-1; ++i) 
			{
				int k = thread_idx/s;
				int j = (thread_idx%s)/grid.ni;
				int i = thread_idx%grid.ni;
				if(k>=1 && k<grid.nk-1 && j>=1 && j<grid.nj-1 && i>=1 && i<grid.ni-1){

					float sum = 0;
					int count = 0;

					if(!old_valid(i,j,k)) 
					{

						if(old_valid(i+1,j,k)) {
							sum += grid(i+1,j,k);
							++count;
						}
						if(old_valid(i-1,j,k)) {
							sum += grid(i-1,j,k);
							++count;
						}
						if(old_valid(i,j+1,k)) {
							sum += grid(i,j+1,k);
							++count;
						}
						if(old_valid(i,j-1,k)) {
							sum += grid(i,j-1,k);
							++count;
						}
						if(old_valid(i,j,k+1)) {
							sum += grid(i,j,k+1);
							++count;
						}
						if(old_valid(i,j,k-1)) {
							sum += grid(i,j,k-1);
							++count;
						}

						//If any of neighbour cells were valid, 
						//assign the cell their average value and tag it as valid
						if(count > 0) {
							temp_grid(i,j,k) = sum /(float)count;
							valid(i,j,k) = 1;
						}

					}
				}
			});
			grid = temp_grid;

		}



		compute_num = u._nx*u._ny*u._nz;
		slice = u._nx*u._ny;
		tbb::parallel_for(0, compute_num, 1, [&](int thread_idx)
		{
			int k = thread_idx/slice;
			int j = (thread_idx%slice)/u._nx;
			int i = thread_idx%u._nx;
			if(k<u._nz && j<u._ny && i<u._nx )
			{
				u(i,j,k) = grid(i,j,k);
			}
		});


	}
	void SmokeSolver3D::getDelta(buffer3Df &field_old, buffer3Df &field_new, buffer3Df &fieldtemp)
	{
		int compute_elements = fieldtemp._blockx*fieldtemp._blocky*fieldtemp._blockz;

		int slice = fieldtemp._blockx*fieldtemp._blocky;

		tbb::parallel_for(0, compute_elements, 1, [&](int thread_idx) {

			uint bk = thread_idx/slice;
			uint bj = (thread_idx%slice)/fieldtemp._blockx;
			uint bi = thread_idx%(fieldtemp._blockx);

			for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
			{
				uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
				if( i<fieldtemp._nx && j<fieldtemp._ny && k<fieldtemp._nz )
				{
					fieldtemp(i,j,k) = field_new(i,j,k) - field_old(i,j,k);
				}
			}
		});
	}
	void SmokeSolver3D::compute_delta_u()
	{
		getDelta(_utemp,_un,_utemp);
		getDelta(_vtemp,_vn,_vtemp);
		getDelta(_wtemp,_wn,_wtemp);
	}

	void SmokeSolver3D::hash_FLIP()
	{
		particle_hash.reserve(_nx*_ny*_nz);
		for (int i=0;i<_nx*_ny*_nz;i++)
		{
			particle_hash[i].resize(0);
			//particle_hash[i].reserve(20);
		}
		for (int i=0;i<FLIP.size();i++)
		{
			Vec3f xx = FLIP[i].pos;
			xx += Vec3f(0.5,0.5,0.5)*_hx;
			Vec3i ijk(xx/_hx);
			if (ijk[0]>=0&&ijk[1]>=0&&ijk[2]>=0&&ijk[0]<_nx&&ijk[1]<_ny&&ijk[2]<_nz)
			{
				int idx = ijk[2]*_nx*_ny + ijk[1]*_nx + ijk[0];
				particle_hash[idx].push_back(i);
			}
		}
	}
	void SmokeSolver3D::reorder_FLIP()
	{
		std::vector<particle> new_FLIP;
		new_FLIP.resize(0);
		//new_FLIP.reserve(FLIP.size());
		for (int k=0;k<_nz;k++)for(int j=0;j<_ny;j++)for(int i=0;i<_nx;i++)
		{
			int idx = k*_nx*_ny + j*_nx + i;
			if (particle_hash[idx].size()<8)
			{
				for (int p=0;p<particle_hash[idx].size();p++)
				{

					new_FLIP.push_back(FLIP[particle_hash[idx][p]]);
				}
				for (int p=0;p<8-particle_hash[idx].size();p++)
				{
					Vec3f x = Vec3f(i,j,k)*_hx + Vec3f(frand(-1,1), frand(-1,1),frand(-1,1))*0.49f*_hx;
					Vec3f vel = get_velocity_self(x);
					new_FLIP.push_back(particle(x[0],x[1],x[2],vel[0],vel[1],vel[2],0,0));
				}
			}
			else if(particle_hash[idx].size()<10)
			{
				for (int p=0;p<particle_hash[idx].size();p++)
				{
					
					new_FLIP.push_back(FLIP[particle_hash[idx][p]]);
				}
			}
			else  
			{
				for (int p=0;p<8;p++)
				{

					new_FLIP.push_back(FLIP[particle_hash[idx][p]]);
				}
			}
		}
		//FLIP.swap(new_FLIP);
		FLIP.resize(0);
		//FLIP.reserve(new_FLIP.size());
		for (int p=0;p<new_FLIP.size();p++)
		{
			if (get_solid_phi(new_FLIP[p].pos)>=0)
			{
				FLIP.push_back(new_FLIP[p]);
			}
		}
	}

	void SmokeSolver3D::particle_from_grid()
	{
		int num = FLIP.size();
		tbb::parallel_for(0,num,1,[&](int thread_idx){
			FLIP[thread_idx].vel += get_velocity_diff(FLIP[thread_idx].pos);
			//FLIP[thread_idx].vel = get_velocity(FLIP[thread_idx].pos);
		});

	}
	void SmokeSolver3D::advect_FLIP(float dt)
	{
		tbb::parallel_for((size_t)0,
				(size_t)FLIP.size(),
				(size_t)1,
				[&](size_t i)
			{
				Vec3f pos = Vec3f(FLIP[i].pos[0],FLIP[i].pos[1],FLIP[i].pos[2]);
				
				//if(get_solid_phi(pos)>0)
				{

					pos = traceRK3(pos,dt);
					/*Vec3f vel = get_velocity(pos);
					Vec3f mpos = pos + 0.5f * dt * vel;
					pos = pos + dt*get_velocity(mpos);*/
					//pos += dt * get_velocity(pos);
					float phi_value = get_solid_phi(pos);
					if (phi_value<0)
					{
						Vec3f normal = get_solid_normal(pos);
						pos -= phi_value*normal;

					}
					FLIP[i].pos = Vec3f(pos[0],pos[1],pos[2]);
				}
				//else
				//{
				//	pos = traceRK3(pos,dt);
				//	/*Vec3f vel = get_velocity(pos);
				//	Vec3f mpos = pos + 0.5f * dt * vel;
				//	pos = pos + dt*get_velocity(mpos);*/
				//	float phi_value = get_solid_phi(pos);
				//	if (phi_value>0)
				//	{
				//		Vec3f normal =get_solid_normal(pos);
				//		pos -= phi_value*normal;

				//	}
				//	FLIP[i].pos = Vec3f(pos[0],pos[1],pos[2]);
				//}
			});
	}

	void SmokeSolver3D::particle_to_grid(buffer3Df & field,buffer3Df &coef, int component)
	{
		field.setZero();
		coef.setZero();
		int compute_elements = field._blockx*field._blocky*field._blockz;

		int slice = field._blockx*field._blocky;

		tbb::parallel_for(0, compute_elements, 1, [&](int thread_idx){
			uint bk = thread_idx/slice;
			uint bj = (thread_idx%slice)/field._blockx;
			uint bi = thread_idx%(field._blockx);

			

			for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
			{
				uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;

				if(i<field._nx&&j<field._ny&&k<field._nz)
				{
					float world_x = ((float)i-field._ox)*_hx;
					float world_y = ((float)j-field._oy)*_hx;
					float world_z = ((float)k-field._oz)*_hx;
					float total_w = 1e-5;
					float total_val = 0;
					for (int hkk=-1; hkk<=1;hkk++) for(int hjj=-1;hjj<=1;hjj++) for(int hii=-1;hii<=1;hii++)
					{
						int iii = (int)i+hii, jjj = (int)j+hjj, kkk = (int)k+hkk;
						if(iii>=0&&iii<_nx &&jjj>=0 && jjj<_ny && kkk>=0 && kkk<_nz)
						{

							int hash_idx = kkk*_nx*_ny + jjj*_nx + iii;
							for (int p=0; p<particle_hash[hash_idx].size();p++)
							{
								int pp = particle_hash[hash_idx][p];
								//if(get_solid_phi(FLIP[pp].pos)>0)
								{

									float weight = compute_weight(world_x,world_y,world_z,
										FLIP[pp].pos[0],FLIP[pp].pos[1],FLIP[pp].pos[2]);

									total_w += weight;

									total_val += weight*(FLIP[pp].vel[component]);
								}


							}
							
						}
					}
					coef(i,j,k) = total_w;
					field(i,j,k) = total_val/total_w;
				}
			
			}
			
		});
	}

	void SmokeSolver3D::seedHalo()
	{
		for(int k=0;k<_nz;k++)for(int j=0;j<_ny;j++)for(int i=0;i<_nx;i++)
		{
			if (sub_domain.is_halo(Vec3f(i,j,k)*_hx, 2.0f*_hx)&&
				(!sub_domain.is_in2(Vec3f(i,j,k)*_hx)))
			{
				for(int r=0;r<512;r++)
				{

					Vec3f x = Vec3f(i,j,k)*_hx + Vec3f(frand(-1,1), frand(-1,1),frand(-1,1))*0.49f*_hx;
					Vec3f vel = get_velocity_self(x);
					FLIP.push_back(particle(x[0],x[1],x[2],vel[0],vel[1],vel[2],0,0));
				}
			}
		}
	}
	void SmokeSolver3D::seed_ghost()
	{
		for(int k=0;k<sub_domain._nz;k++)for(int j=0;j<sub_domain._ny;j++)for(int i=0;i<sub_domain._nx;i++)
		{
			if(get_solid_phi(Vec3f(i,j,k)*sub_domain._hx + sub_domain.bmin)<3.0f*sub_domain._hx)
			{

				for(int r=0;r<8;r++)
				{

					Vec3f x = Vec3f(i,j,k)*sub_domain._hx + Vec3f(frand(-1,1), frand(-1,1),frand(-1,1))*0.49f*sub_domain._hx + sub_domain.bmin;
					float phi = sub_domain.get_solid_phi(x);
					
					if(phi<=2.0*sub_domain._hx&&phi>-2.0*sub_domain._hx)
					{
						//Vec3f normal = _nodal_solid_phi.sample_grad(x[0],x[1],x[2]);
						
						//normalize(normal);
						//Vec3f normal = get_solid_normal(x);
						//x = x - phi*normal;
						//Vec3f solid_vel = get_solid_vel(x);
						//Vec3f vel = solid_vel;//- get_solid_vel(x-2.0f*phi*normal);;
						Vec3f vel = get_velocity(x);
						if(phi<=0) vel = get_solid_vel(x);
						sub_domain.FLIP.push_back(particle(x[0],x[1],x[2],vel[0],vel[1],vel[2],0,0));
					}
					
				}
			}
		}
	}
	void SmokeSolver3D::diffuse(float dt)
	{
		int num = sub_domain.FLIP.size();
		printf("after merge subdomain: %d\n",num);
		seed_ghost();
		printf("after seed ghost:%d\n",sub_domain.FLIP.size());
		sub_domain.hash_FLIP();
		int compute_num = sub_domain.FLIP.size();
		vector<Vec3f> du;
		du.resize(compute_num);
		//tbb::parallel_for(0, compute_num, 1, [&](int i){
		int count = 0;
		for(int i = 0; i < compute_num; i++) {
			Vec3f xi = sub_domain.FLIP[i].pos;
			Vec3f ui = sub_domain.FLIP[i].vel;
			Vec3i ijk = Vec3i((xi-sub_domain.bmin)/sub_domain._hx+Vec3f(0.5,0.5,0.5));
			du[i] = Vec3f(0.0f,0.0f,0.0f);
			if ( get_solid_phi(xi)<=2.0*_hx
				&& ijk[0]>=0&&ijk[0]<sub_domain._nx
				&&ijk[1]>=0&&ijk[1]<sub_domain._ny
				&&ijk[2]>=0&&ijk[2]<sub_domain._nz)
			{
				count++;
				float coef = 1e-9;
				//coef += PSE_kernel(xi,xi,dt*nu+1e-10);
				Vec3f u_exchange(0);
				for (int kk = max(ijk[2]-1,0);kk<=min(ijk[2]+1,(int)sub_domain._nz-1);kk++)
				{

					for (int jj = max(ijk[1]-1,0);jj<=min(ijk[1]+1,(int)sub_domain._ny-1);jj++)
					{

						for (int ii = max(ijk[0]-1,0);ii<=min(ijk[0]+1,(int)sub_domain._nx-1);ii++)
						{
							int hash_idx = kk*sub_domain._nx*sub_domain._ny + jj*sub_domain._nx + ii;
							for (int p=0; p<(int)(sub_domain.particle_hash[hash_idx].size());p++)
							{
								//int pp = rand()%(particle_hash[hash_idx].size());
								
								particle pj;
								pj.pos = sub_domain.FLIP[sub_domain.particle_hash[hash_idx][p]].pos;
								pj.vel = sub_domain.FLIP[sub_domain.particle_hash[hash_idx][p]].vel;
								u_exchange += (pj.vel - ui) * (float) PSE_kernel(xi,pj.pos,dt*nu+1e-10);
								coef += PSE_kernel(xi,pj.pos,dt*nu+1e-10);
							}
						}
					}
				}
				du[i] = u_exchange/coef;
			}

		}//);
		printf("*************%d\n", count);
		tbb::parallel_for(0, compute_num, 1, [&](int i){
			sub_domain.FLIP[i].vel += du[i];
			/*if (get_solid_phi(sub_domain.FLIP[i].pos)<=0)
			{
				sub_domain.FLIP[i].vel = get_solid_vel(sub_domain.FLIP[i].pos);
			}*/
		});
		sub_domain.FLIP.resize(num);
		printf("after diffusion:%d\n",sub_domain.FLIP.size());
	}

	void SmokeSolver3D::gen_heat(float dt)
	{

		int compute_num = _rho._nx*_rho._ny*_rho._nz;
		int slice = _rho._nx*_rho._ny;
		tbb::parallel_for(0, compute_num, 1, [&](int thread_idx)
		{
			int k = thread_idx/slice;
			int j = (thread_idx%slice)/_rho._nx;
			int i = thread_idx%_rho._nx;
			if(k<_rho._nz && j<_rho._ny && i<_rho._nx )
			{
				_rho(i,j,k) =_rho(i,j,k)/(1.0f + _dens_decay * dt);
				_Temperature(i,j,k) = _Temperature(i,j,k)/(1.0f + _temp_decay * dt);
				
				
				Vec3f pos = Vec3f(i,j,k)*(float)(_rho._hx);
				
				

				if(dist(Vec3f(pos[0],0,pos[2]),Vec3f(_smoke_center[0],0,_smoke_center[2]))<=_smoke_radius && fabs(pos[1]-_smoke_center[1])<=0.25*_smoke_radius)
				{
					_rho(i,j,k) = _smoke_dens;
					_Temperature(i,j,k) = _smoke_heat;
				}
			}
		});
	}
	void SmokeSolver3D::add_smokeforce(float dt)
	{
		tbb::parallel_for((size_t)0,
				(size_t)FLIP.size(),
				(size_t)1,
				[&](size_t i)
			{
				Vec3f pos = Vec3f(FLIP[i].pos[0],FLIP[i].pos[1],FLIP[i].pos[2]);
				
				float rho = _rho.sample_linear(pos[0],pos[1],pos[2]);
				float T   = _Temperature.sample_linear(pos[0],pos[1],pos[2]);

				float force_y = dt *(-_alpha*rho + _beta*T);
			
				FLIP[i].vel += Vec3f(0,force_y,0);
			});
	}
}