#ifndef _POST_REFINEMENT_
#define _POST_REFINEMENT_

#include "Smoke_solver3D.h"

using namespace std;
namespace ivock{
	class PostResimulation
	{
	public:
		PostResimulation(){}
		~PostResimulation(){ fine_solver.Finalize(); coarse_solver.Finalize();}


		SmokeSolver3D fine_solver;
		SmokeSolver3D coarse_solver;
		buffer3Df U,V,W, Un, Vn, Wn;
		Vec3f getVelocity(Vec3f & pos1, Vec3f & pos2)
		{
			float u = U.sample_linear(pos1[0],pos1[1],pos1[2]) 
				+ fine_solver._un.sample_linear(pos2[0],pos2[1],pos2[2]);

			float v = V.sample_linear(pos1[0],pos1[1],pos1[2]) 
				+ fine_solver._vn.sample_linear(pos2[0],pos2[1],pos2[2]);

			float w = W.sample_linear(pos1[0],pos1[1],pos1[2]) 
				+ fine_solver._wn.sample_linear(pos2[0],pos2[1],pos2[2]);

			return Vec3f(u,v,w);

		}
		void advect_field(float dt, buffer3Df & field, buffer3Df &field_new)
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

						float world_x = ((float)i-field._ox)*fine_solver._hx + fine_solver.origin[0];
						float world_y = ((float)j-field._oy)*fine_solver._hy + fine_solver.origin[1];
						float world_z = ((float)k-field._oz)*fine_solver._hz + fine_solver.origin[2];
						Vec3f pos1(world_x,world_y,world_z);
						Vec3f pos2 = pos1 - fine_solver.origin;
						//Vec3f trace_pos = trace(dt, pos);

						Vec3f vel = getVelocity(pos1,pos2);
						

						float px = world_x - dt * vel[0], py = world_y - dt *vel[1], pz = world_z - dt*vel[2];

						float SLv = field.sample_linear(px - fine_solver.origin[0],
							py - fine_solver.origin[1],
							pz - fine_solver.origin[2]);
						if (fine_solver._b_desc.sample_linear(px - fine_solver.origin[0],
							py - fine_solver.origin[1],
							pz - fine_solver.origin[2])==2)
						{
							SLv = 0;
						}
						field_new(i,j,k) = SLv;

					}
				}
			});
		}
		void advect(float dt)
		{
			advect_field(dt, fine_solver._un, fine_solver._unp1);
			advect_field(dt, fine_solver._vn, fine_solver._vnp1);
			advect_field(dt, fine_solver._wn, fine_solver._wnp1);
		}
		void lerp_buffer(float c, buffer3Df &un, buffer3Df &unp1, buffer3Df &u)
		{
			int compute_elements = u._blockx*u._blocky*u._blockz;

			int slice = u._blockx*u._blocky;

			tbb::parallel_for(0, compute_elements, 1, [&](int thread_idx) {

				uint bk = thread_idx/slice;
				uint bj = (thread_idx%slice)/u._blockx;
				uint bi = thread_idx%(u._blockx);

				for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
				{
					uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
					if(i<u._nx && j<u._ny && k<u._nz)
					{

						u(i,j,k) = (1-c)*un(i,j,k) + c*unp1(i,j,k);

					}
				}
			});
		}
		void coarse_fine_addition(buffer3Df & coarse, buffer3Df &fine, float alpha)
		{
			int compute_elements = fine._blockx*fine._blocky*fine._blockz;

			int slice = fine._blockx*fine._blocky;

			tbb::parallel_for(0, compute_elements, 1, [&](int thread_idx) {

				uint bk = thread_idx/slice;
				uint bj = (thread_idx%slice)/fine._blockx;
				uint bi = thread_idx%(fine._blockx);

				for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
				{
					uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
					if(i<fine._nx && j<fine._ny && k<fine._nz)
					{

						float world_x = ((float)i-fine._ox)*fine_solver._hx + fine_solver.origin[0];
						float world_y = ((float)j-fine._oy)*fine_solver._hy + fine_solver.origin[1];
						float world_z = ((float)k-fine._oz)*fine_solver._hz + fine_solver.origin[2];
						


						fine(i,j,k) += alpha*coarse.sample_linear(world_x,world_y,world_z);

					}
				}
			});
		}
		void TimeStep(float dt)
		{
			//Un = coarse_solver.U;
			Un.copy(coarse_solver._un);
			Vn.copy(coarse_solver._vn);
			Wn.copy(coarse_solver._wn);

			coarse_solver.time_step(dt,0);

			//Unp1 = coarse_solver.Un;

			{
				//advect fine solver,
				lerp_buffer(0.25*(float)i,Un,coarse_solver._un,U);
				lerp_buffer(0.25*(float)i,Vn,coarse_solver._vn,V);
				lerp_buffer(0.25*(float)i,Wn,coarse_solver._wn,W);
				// u_np1 = advect u_n by (U + u_n);
				advect(0.25*dt);
				//diffuse velocity
				//fine_solver.diffuse_buffer(0.01,fine_solver._unp1);
				//fine_solver.diffuse_buffer(0.01,fine_solver._vnp1);
				//fine_solver.diffuse_buffer(0.01,fine_solver._wnp1);
				fine_solver._un.copy(fine_solver._unp1);
				fine_solver._vn.copy(fine_solver._vnp1);
				fine_solver._wn.copy(fine_solver._wnp1);

				lerp_buffer(0.25*(float)i + 0.25,Un,coarse_solver._un,U);
				lerp_buffer(0.25*(float)i + 0.25,Vn,coarse_solver._vn,V);
				lerp_buffer(0.25*(float)i + 0.25,Wn,coarse_solver._wn,W);


				//fine_solver.u += U_np1;
				coarse_fine_addition(U,fine_solver._un,1);
				coarse_fine_addition(V,fine_solver._vn,1);
				coarse_fine_addition(W,fine_solver._wn,1);
				//fine_solver.diffuse_buffer(0.01,fine_solver._un);
				//fine_solver.diffuse_buffer(0.01,fine_solver._vn);
				//fine_solver.diffuse_buffer(0.01,fine_solver._wn);
				//fine_solver.project;
				fine_solver.pcg_projection();
				//fine_solver.u -= U_np1;

				fine_solver.emit_tracers();
				fine_solver.advect_tracers(0.25*dt);


				coarse_fine_addition(U,fine_solver._un,-1);
				coarse_fine_addition(V,fine_solver._vn,-1);
				coarse_fine_addition(W,fine_solver._wn,-1);
				fine_solver.clearBoundary();
			}
			


		}

	};
}





#endif